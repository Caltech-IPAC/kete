from __future__ import annotations
import numpy as np
import pandas as pd  # type: ignore
import logging

from . import utils
from .definitions import (
    neo_amor_complete,
    neo_apollo_complete,
    neo_aten_complete,
    neo_complete,
    nearly_neo_complete,
    neo_atira,
)
from .power_law import CumulativePowerLaw
from ..conversion import compute_H
from ..time import Time
from ..vector import CometElements, State
from ..mpc import num_to_base62
from ..flux import NeatmParams


logger = logging.getLogger(__name__)


neo_sdt_diameters = CumulativePowerLaw(
    slopes=[3.2, 1.6217413993016, 2.75],
    cutoffs=[0.07, 1.5],
    num_1km_obj=993,
    max_possible=50,
)
"""
Diameter fit from the 2017 SDT Report
The middle slope in the paper is 1.6, here it is set to 1.64 to make the total
number of objects larger than 140m be 25k objects. This is within the error bars of
the paper."""


class NEOGroup:
    """
    Parameters
    ----------
    name :
        The name of the group.
    group_id:
        An id string to be inserted into the object name.
    frac_dark : float
        The percentage of objects which come from the darker of the two Rayleigh
        distributions of visible albedo.
    filter_func:
        The orbital filter function for this group as defined in
        `neospy.population.definitions`
    pop_fraction : float
        The relative fraction of the total NEO population which is made up of the
        specified `filt` function, for example, the Apollos are about `0.551`.

        Population fractions are from Granvik 2018
        Debiased orbit and absolute-magnitude distributions for near-Earth objects
        https://www.boulder.swri.edu/~bottke/Reprints/Granvik_2018_Icarus_312_181_Debiased_Orbit_Mag_NEO.pdf

        That papers relative fractions dont add up to 100%, so these values have been
        normalized to make them 100%, and Vatiras count was added to the Atiras
    """

    def __init__(self, name, group_id, frac_dark, filter_func, pop_fraction):
        self.name = name
        self.group_id = group_id
        self.frac_dark = frac_dark
        self.filter_func = filter_func
        self.pop_fraction = pop_fraction

    def sample_orbits(self, n_samples, seed=42):
        orbits = utils.sample_orbits(
            nearly_neo_complete, n_samples * 10, seed=seed, bandwidth=0.2
        )
        keep = self.filter_func(orbits.peri_dist, orbits.ecc, np.zeros(len(orbits)))
        orbits = orbits[keep][:n_samples]

        # The filter specified is too small of a percentage of the full dataset to do
        # sampling and selection. So here we sample directly from this set. This is
        # primarily for the Atiras, since they are such a tiny fraction.
        if len(orbits) < n_samples:
            orbits = utils.sample_orbits(
                self.filter_func, n_samples, seed=seed, bandwidth=0.2
            )
        return orbits

    def sample_beaming(self, n_samples, seed=42):
        return utils.sample_beaming(
            "neowise_neos.xml",
            neo_complete,
            n_samples,
            seed=seed + 2,
            bandwidth=0.1,
        )

    def sample_vis_albedo(self, n_samples, seed=42):
        """
        Sample visible albedos using the double Rayleigh distribution described in the
        2016 Wright paper: https://arxiv.org/pdf/1606.07421.pdf

        The double Rayleigh distribution is fully described by 3 values, the sigma
        values for each distribution, and the relative fraction of one distribution vs
        the other.

        The sigmas are hard coded from the paper above, but its expected that each NEO
        group will have a different relative share of "dark" vs "light" albedos.

        Parameters
        ----------
        n_samples : int
            The number of visible and ir albedos to compute.
        seed : int
            The random seed used for all generation of random numbers. This can be used
            to accurately reconstruct the same objects.
        """
        dark_sigma = 0.030
        bright_sigma = 0.168
        rng = np.random.default_rng(seed)

        # This breaks the sampling problem into 2 bins, if the object came from the
        # light or dark distributions.
        # Each Rayleigh distribution is invertible on its own, but together it doesn't
        # have a nice analytical form. By splitting it up it makes it trivial to sample.
        in_pop_a = rng.random(n_samples) <= self.frac_dark

        # Using CDF inverse transformation sampling on each population
        sample_a = dark_sigma * np.sqrt(-2 * np.log(1 - rng.random(n_samples)))
        sample_b = bright_sigma * np.sqrt(-2 * np.log(1 - rng.random(n_samples)))

        # Selecting the samples based off if it came from the first or second
        # distribution
        sample = sample_a.copy()
        sample[~in_pop_a] = sample_b[~in_pop_a]

        return sample

    def sample_albedos(self, n_samples, bandwidth=0.1, seed=42):
        """
        Sample NEO albedos using existing data, returns a pandas dataframe of 2 columns.

        The visible albedos are sampled using the `neo_sample_vis_albedo` function,
        which accepts a single parameter, `dark_frac`. This is the relative fraction of
        objects which come from the darker of the two Rayleigh distributions, other than
        that, the visible albedos are fully defined.

        The `ir_albedo` is created by taking the full dataset of all observationally
        complete NEOs and computing the `IR_Albedo / V_Albedo` distribution, then
        sampling from that. These values are multiplied by the `V_albedo` distribution
        from the analytically defined function `neo_sample_vis_albedo` to estimate an
        `IR_albedo`.

        Then there is a rejection process where the albedos are required to be within a
        physically reasonably bound, see the code for details.

        Parameters
        ----------
        n_samples : int
            The number of visible and ir albedos to compute.
        bandwidth : float
            The number bandwidth of the Guassian KDE used to sample the IR/Vis albedo
            ratio.
        seed : int
            The random seed used for all generation of random numbers. This can be used
            to accurately reconstruct the same objects.
        """
        albedo_ratio_data = utils.pds_data_filtered(
            "neowise_neos.xml", neo_complete, fit_code="VBI"
        )

        for scale in [1.5, 3, 10]:
            extra_samples = int(np.ceil(n_samples * scale))
            v_albedo = self.sample_vis_albedo(extra_samples, seed=seed)
            ir_ratios = utils.sample_values(
                np.transpose(albedo_ratio_data[["V_albedo", "IR_albedo"]]),
                extra_samples,
                bandwidth=bandwidth,
                seed=seed + 1,
            )
            ir_ratios = np.array(ir_ratios[1] / ir_ratios[0]).ravel()
            ir_albedo = v_albedo * ir_ratios

            keep = (v_albedo > 0.015) & (v_albedo < 0.7)
            keep = keep & (ir_albedo > 0.015) & (ir_albedo < 0.85)
            v_albedo = v_albedo[keep][:n_samples]
            ir_albedo = ir_albedo[keep][:n_samples]
            if len(v_albedo) == n_samples:
                return pd.DataFrame(
                    np.transpose([v_albedo, ir_albedo]),
                    columns=["V_albedo", "IR_albedo"],
                )
        raise ValueError(
            "Filter was too aggressive with provided data, failed to find valid sample."
        )

    def sample_diameters(self, min_size, max_size=np.inf, n_objects=None, seed=42):
        """
        Generate a sorted fair sample of diameters for this group.

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
        seed : int
            The random seed used for all generation of random numbers. This can be used
            to accurately reconstruct the same objects.
        """
        scale = self.pop_fraction if n_objects is None else None
        return neo_sdt_diameters.sample(
            min_size, max_size=max_size, scale=scale, seed=seed, n_objects=n_objects
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
    ):
        """
        Construct a generator of an NEO population down to the specified diameter size.

        To reduce memory pressure, this is structured as a generator function, and all
        objects can be collected using something like:

        .. testcode::
            :skipif: True

            for obj_batch in neo_sample_objects(filt, min_size, 0.254):
                # Do something with the batch of object properties

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
        """

        if epoch is None:
            epoch = Time.from_current_time().jd
        prefix = "S" + num_to_base62(seed, 2) + self.group_id

        diams = self.sample_diameters(
            min_size, seed=seed, n_objects=n_objects, max_size=max_size
        )
        logger.info("Generated % diameters", len(diams))

        n_batches = int(np.ceil(len(diams) / batch_size))

        for i in range(n_batches):
            batch = diams[i * batch_size : (i + 1) * batch_size]
            batch_len = len(batch)
            ids = np.arange(i * batch_size, (i + 1) * batch_size)[:batch_len] + id_start

            albedos = self.sample_albedos(batch_len, seed=seed + 1)
            beaming = self.sample_beaming(batch_len, seed=seed + 2)
            orbits = self.sample_orbits(batch_len, seed=seed + 3)

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
                        orbits.epoch,
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
                    "epoch",
                ],
            )

            objects["h_mag"] = compute_H(
                albedo=objects.vis_albedo, diameter=objects.diameter
            )
            objects["g_phase"] = 0.15
            objects["emissivity"] = 0.9
            objects["desig"] = [prefix + num_to_base62(n, 6) for n in ids]
            objects["efrho"] = 0.0

            yield objects


Amor = NEOGroup(
    name="Amor",
    group_id="03",
    frac_dark=0.358,
    filter_func=neo_amor_complete,
    pop_fraction=0.400,
)
"""Sampler for the Amor Group"""

Aten = NEOGroup(
    name="Aten",
    group_id="01",
    frac_dark=0.104,
    filter_func=neo_aten_complete,
    pop_fraction=0.035,
)
"""Sampler for the Aten Group"""

Apollo = NEOGroup(
    name="Apollo",
    group_id="02",
    frac_dark=0.308,
    filter_func=neo_apollo_complete,
    pop_fraction=0.551,
)
"""Sampler for the Apollo Group"""

Atira = NEOGroup(
    name="Atira",
    group_id="00",
    frac_dark=0.104,
    filter_func=neo_atira,
    pop_fraction=0.014,
)
"""Sampler for the Atria Group"""


def create_neo_population(
    filename: str,
    min_size: float = 0.14,
    max_size: float = np.inf,
    seed: int = 42,
    epoch=None,
):
    """
    Create a new NEO population containing objects down to the specified minimum
    diameter in km.

    .. testcode::
        :skipif: True

        create_neo_population("neos")

    Parameters
    ----------
    filename :
        The filename for the new database.
    min_size :
        The minimum diameter of objects to make, this is in units of km.
    max_size :
        The maximum diameter of objects to make, this is in units of km.
    seed :
        The random seed used for all generation of random numbers. The same seed and min
        diameter will result in the same database being built.
    epoch : None or float
        The epoch of observation for the database. If no epoch is specified, the
        database will be built using the first day of the survey. This is the ideal time
        to minimize simulation time.
    """
    logger.info("Beginning construction of database...")

    if epoch is None:
        epoch = Time.from_current_time().jd

    id_start = 0

    names = []
    pos = []
    vel = []
    properties: list = []
    for group in [Atira, Aten, Apollo, Amor]:
        logger.info("Working on %s", group.name)
        for vals in group.sample_objects(
            min_size=min_size,
            max_size=max_size,
            batch_size=1_000_000,
            seed=seed,
            id_start=id_start,
        ):
            vals = vals.copy()
            logger.info(
                "%d  %s  %0.1f+ m diameter",
                len(vals),
                group.name,
                min(vals.diameter * 1000),
            )
            state = CometElements(
                epoch,
                vals.desig,
                vals.ecc,
                vals.incl,
                vals.peri_dist,
                vals.peri_arg,
                vals.peri_time,
                vals.lon_node,
            ).as_state
            names.extend(vals.desig)
            pos.append(state.pos)
            vel.append(state.vel)
            properties.extend(
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
            )

            id_start += len(vals)
    pos_arr = np.vstack(pos)
    vel_arr = np.vstack(vel)

    state = State(epoch, names, pos_arr, vel_arr, [None for _ in names])
    state.save(f"{filename}_state.bin")
    NeatmParams.save_list(properties, f"{filename}_props.bin")
