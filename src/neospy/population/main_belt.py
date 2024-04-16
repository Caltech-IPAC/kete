from __future__ import annotations
import numpy as np
import pandas as pd  # type: ignore
import logging
from scipy.optimize import least_squares  # type: ignore

from ..conversion import compute_H
from ..time import Time
from ..vector import CometElements
from ..flux import NeatmParams


from .utils import (
    sample_albedos,
    sample_beaming,
    sample_orbits,
    visibility_test,
)
from .diameters import CumulativePowerLaw
from .definitions import (
    mba_inner_complete,
    mba_middle_complete,
    mba_outer_complete,
)
from ..mpc import num_to_base62
from ..pds import pds_data_filtered

logger = logging.getLogger(__name__)


mba_inner_diameters = CumulativePowerLaw(
    slopes=[2.87492, 1.0, 2.88164],
    cutoffs=[13.4677, 82.9235],
    num_1km_obj=257038,
    max_possible=196.371,
)
"""
Diameter fit for the inner main belt, these were fit using the
`mba_fit_diameter_distribution` function defined below, and the parameters are saved
here as reference.
"""

mba_middle_diameters = CumulativePowerLaw(
    slopes=[2.89788, 1.03716, 5.0],
    cutoffs=[13.06984, 115.4161],
    num_1km_obj=928025,
    max_possible=231.689,
)
"""
Diameter fit for the middle main belt, these were fit using the
`mba_fit_diameter_distribution` function defined below, and the parameters are saved
here as reference.
"""

mba_outer_diameters = CumulativePowerLaw(
    slopes=[2.54969, 1.21393, 3.2284],
    cutoffs=[19.58096, 87.5279],
    num_1km_obj=1816513,
    max_possible=453.239,
)
"""
Diameter fit for the outer main belt, these were fit using the
`mba_fit_diameter_distribution` function defined below, and the parameters are saved
here as reference.
"""


def mba_fit_diameter_distribution(filt):
    """
    Creates cumulative object count functions which define the number of MBAs greater
    than a specified diameter in km.

    This returns a CumulativePowerLaw class.

    This class enables estimating the number of objects larger than a specified
    diameter. They are constructed using the PDS dataset for the main belt, and the
    data is linearly extrapolated to down to arbitrarily small diameters.

    Parameters
    ----------
    filt :
        The `neospy.population.definition` function which selects an population group.
        This should be `mba_inner_complete` or similar.
    """

    # Plotting the cumulative total number of objects as a function of size shows that
    # the diameters seem to be split into 3 distinct regimes, along with a falloff of
    # small object counts due to them being too small to observe. Each of these 3
    # zones follow a linear slope, along with a cutoff size where no objects are seen
    # below that size.

    diams = pds_data_filtered("neowise_mainbelt.xml", filt, "D").Diameter

    # Finding the minimum size which exists in the dataset, this uses an approximation
    # of the hessian to determine when the population stops growing as the diameters
    # get smaller.
    cutoff = []
    for i in range(50, 200):
        bins = np.logspace(np.log10(min(diams)), np.log10(max(diams)), i)
        counts, bins = np.histogram(diams, bins=bins)
        bin_centers = (bins[1:] + bins[:-1]) / 2
        cutoff.append(bin_centers[np.argmin(np.diff(counts))])
    cutoff = np.median(cutoff)

    # Now that the observability cutoff has been estimated, we can find the slope of the
    # 3 distinct regions.  Select and histogram the data after this cutoff
    bins = np.linspace(cutoff, max(diams), 50)
    counts, bins = np.histogram(diams, bins=bins)

    bounds = ([1.0, 50, 0, 1, 1, 1], [50, 500, 1e12, 5, 5, 5])

    def _diff(params, diams, freq):
        cutoff_0, cutoff_1, offset_0, slope_0, slope_1, slope_2 = params
        t = CumulativePowerLaw(
            slopes=[slope_0, slope_1, slope_2],
            cutoffs=[cutoff_0, cutoff_1],
            num_1km_obj=offset_0,
            max_possible=np.inf,
        )
        return np.array(freq) - t.num_larger_than(diams)

    fit = least_squares(
        _diff,
        [15.0, 100.0, 300000.0, 2.6, 1.2, 4.7],
        args=(bins[:-1], np.cumsum(counts[::-1])[::-1]),
        bounds=bounds,
    ).x
    cutoff_0, cutoff_1, offset_0, slope_0, slope_1, slope_2 = fit

    return CumulativePowerLaw(
        slopes=[slope_0, slope_1, slope_2],
        cutoffs=[cutoff_0, cutoff_1],
        num_1km_obj=offset_0,
        max_possible=np.max(diams),
    )


class MBAGroup:
    def __init__(self, name, group_id, filter_func, diameter_sampler):
        self.name = name
        self.group_id = group_id
        self.diameter_sampler = diameter_sampler
        self.filter_func = filter_func

    def sample_orbits(self, n_samples, seed=42):
        return sample_orbits(self.filter_func, n_samples, seed=seed)

    def sample_beaming(self, n_samples, seed=43):
        return sample_beaming(
            "neowise_mainbelt.xml", self.filter_func, n_samples, seed=seed
        )

    def sample_albedos(self, n_samples, seed=44):
        return sample_albedos(
            "neowise_mainbelt.xml", self.filter_func, n_samples, seed=seed
        )

    def sample_diameters(self, min_size, max_size=np.inf, n_objects=None, seed=42):
        scale = None if n_objects is not None else 1.0
        return self.diameter_sampler.sample(
            min_size, max_size=max_size, seed=seed, n_objects=n_objects, scale=scale
        )

    def sample_objects(
        self,
        min_size,
        max_size=np.inf,
        n_objects=None,
        batch_size=100_000,
        seed=42,
        id_start=0,
        epoch=None,
        keep_only_visible=True,
    ):
        """
        Construct a generator of an MBA population down to the specified diameter size.

        To reduce memory pressure, this is structured as a generator function.

        Parameters
        ----------
        min_size : float
            The minimum diameter of objects to make, this is in units of km.
        max_size : float
            The maximum diameter of objects to make, units of km.
        n_objects:
            If this is None, then the expected number of objects are made for the
            specified size range, if this is set to a value, then exactly that number
            of objects are made.
        batch_size : int
            The number of objects inside of each generated pandas dataframe.
        seed : int
            The random seed used for all generation of random numbers. This can be used
            to accurately reconstruct the same objects.
        id_start : int
            Each object has a object id assigned to it, this is the offset value to add
            to all objects created in this call.
        keep_only_visible : bool
            If this is `True`, only objects which are potentially visible are kept.
        """

        if epoch is None:
            epoch = Time.from_current_time().jd
        prefix = "S" + num_to_base62(seed, 2) + self.group_id

        diams = self.sample_diameters(
            min_size, seed=seed, n_objects=n_objects, max_size=max_size
        )
        n_obj = len(diams)
        logger.info("Generated % diameters", n_obj)

        n_batches = int(np.ceil(len(diams) / batch_size))

        for i in range(n_batches):
            batch = diams[i * batch_size : (i + 1) * batch_size]
            batch_len = len(batch)
            ids = np.arange(i * batch_size, (i + 1) * batch_size)[:batch_len] + id_start

            albedos = self.sample_albedos(batch_len, seed=seed + 1 + i * 6)
            beaming = self.sample_beaming(batch_len, seed=seed + 2 + i * 6)
            orbits = self.sample_orbits(batch_len, seed=seed + 3 + i * 6)

            objects = pd.DataFrame(
                np.transpose(
                    [
                        beaming.Beaming_param,
                        albedos.V_albedo,
                        albedos.IR_albedo,
                        batch,
                        orbits.peri_dist,
                        orbits.ecc,
                        orbits.incl,
                        orbits.lon_node,
                        orbits.peri_arg,
                        orbits.peri_time,
                    ]
                ),
                columns=[
                    "beaming",
                    "vis_albedo",
                    "ir_albedo",
                    "diameter",
                    "peri_dist",
                    "ecc",
                    "incl",
                    "lon_node",
                    "peri_arg",
                    "peri_time",
                ],
            )

            objects["h_mag"] = compute_H(
                albedo=objects.vis_albedo, diameter=objects.diameter
            )
            objects["g_phase"] = 0.15
            objects["desig"] = [prefix + num_to_base62(n, 6) for n in ids]
            objects["epoch"] = epoch
            objects["efrho"] = 0.0

            if keep_only_visible:
                vis = visibility_test(
                    objects.peri_dist,
                    objects.vis_albedo,
                    objects.beaming,
                    objects.diameter,
                )
                objects = objects[vis]

            yield objects


MBA_Inner = MBAGroup(
    name="Inner Main Belt",
    group_id="10",
    filter_func=mba_inner_complete,
    diameter_sampler=mba_inner_diameters,
)


MBA_Middle = MBAGroup(
    name="Middle Main Belt",
    group_id="11",
    filter_func=mba_middle_complete,
    diameter_sampler=mba_middle_diameters,
)


MBA_Outer = MBAGroup(
    name="Outer Main Belt",
    group_id="12",
    filter_func=mba_outer_complete,
    diameter_sampler=mba_outer_diameters,
)


def create_mba_population(
    db_name,
    min_sizes=None,
    seed=43,
    epoch=None,
    batch_size=1_000_000,
    keep_only_visible=True,
):
    """
    Create a new databases using the samplers above.

    There will be a number of databases for each MBA group (inner, middle, outer).

    Objects will not be placed into the database if they are never visible from the NEO
    Surveyor, this results in significant reduction of the total number of objects
    saved.

    These databases will at most contain ``batch_size`` number of objects. This is to
    reduce the memory pressure during simulation, 100_000 will result in around 50 gigs
    of memory used during a regular simulation run. Something around 25k is appropriate
    for most laptops.

    Databases are constructed with the provided name plus `_0001.db` added on the end.

    .. testcode::
        :skipif: True

        # Create all objects for the 3 sections of the main belt, larger than 3, 4, and
        # 7 km in size for the 3 groups.
        create_mba_population("mbas", min_size=[3, 4, 7])

    Parameters
    ----------
    db_name : str
        The filename for the new database.
    min_sizes : float
        A list of 3 minimum diameters of objects for the inner/middle/outer MBAs, this
        is in units of km. For example: [1, 1.5, 3] for 1km+ IMB, 1.5km+ MMB, 3km+ OMB
    seed : int
        The random seed used for all generation of random numbers. The same seed and min
        diameter will result in the same database being built.
    epoch : None or float
        The epoch of observation for the database. If no epoch is specified, the
        database will be built using the first day of the survey. This is the ideal time
        to minimize simulation time.
    batch_size : int
        The largest number of objects to place inside of a database before creating a
        new database.
    keep_only_visible : bool
        If this is `True`, only objects which are potentially visible are kept.
    """
    min_sizes = [0.4, 0.5, 0.75] if min_sizes is None else min_sizes
    logger.info("Beginning construction of database...")

    if epoch is None:
        epoch = Time.from_current_time().jd

    diam_samplers = [MBA_Inner, MBA_Middle, MBA_Outer]

    db_idx = 0
    for min_size, sampler in zip(min_sizes, diam_samplers):
        logger.info("Working on %s", sampler.name)
        id_start = 0
        for vals in sampler.sample_objects(
            min_size=min_size,
            batch_size=batch_size,
            seed=seed + db_idx,
            id_start=id_start,
            keep_only_visible=keep_only_visible,
        ):
            group_name = sampler.name.split(" ")[0].lower()
            full_name = db_name + f"_{group_name}_{str(db_idx).zfill(4)}"
            id_start += len(vals)
            db_idx += 1
            state = CometElements(
                epoch,
                list(vals.desig),
                vals.ecc,
                vals.incl,
                vals.peri_dist,
                vals.peri_arg,
                vals.peri_time,
                vals.lon_node,
            ).as_state

            # pylint: disable=too-many-function-args
            properties_list = [
                NeatmParams.new_neos(
                    o.desig,
                    [o.ir_albedo, o.ir_albedo],
                    o.h_mag,
                    o.diameter,
                    o.v_albedo,
                    o.beaming,
                    o.g_phase,
                )
                for o in vals.itertuples()
            ]

            state.save(f"{full_name}_state.bin")
            NeatmParams.save_list(properties_list, f"{full_name}_props.bin")

            logger.info(
                "Completed database %s of %d objects %0.2f+ km.",
                full_name,
                len(vals),
                min(vals.diameter),
            )
        logger.info("Finished %s, constructed %d objects", sampler.name, id_start)
