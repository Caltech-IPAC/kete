from __future__ import annotations
from ..mpc import mpc_known_orbit_filtered
from ..pds import pds_data_filtered
from ..vector import Vector
from ..flux import neatm
import numpy as np
from scipy.stats import gaussian_kde  # type: ignore
import pandas as pd  # type: ignore
from numpy.typing import NDArray
from typing import Optional


def sample_values(
    data: NDArray[np.floating],
    n_samples: float,
    bandwidth: Optional[float] = None,
    seed: int = 42,
) -> np.ndarray:
    """
    Create new samples which approximate the distribution of the provided data using
    gaussian kernel density estimation (KDE).

    Data is normalized against its standard deviations, making the bandwidth parameter
    be meaningful across a large scale difference in inputs.

    Parameters
    ----------
    data:
        An array of values which can be cast to a numpy array from which to generate
        samples.
    n_samples:
        How many samples to generate from the provided dataset.
    bandwidth:
        The bandwidth of the gaussian KDEs in fractions of standard deviation of the
        data.
    seed:
        A random seed to use to generate the values, the same seed will result in the
        same draws.
    """
    std = np.array(np.std(data, axis=1))[:, np.newaxis]
    data /= std
    kde = gaussian_kde(data, bw_method=bandwidth)
    return kde.resample(n_samples, seed=seed) * std


def sample_orbits(filt, n_samples, bandwidth=0.05, seed=42):
    """
    Sample the MPC orbital database using the specified filter function and generate
    new orbital elements which also obey the filter.

    .. testcode::
        :skipif: True

        mpc_filter = neospy.population.definitions.complete_mba_inner_filter
        new = neospy.population.utils.sample_orbits(mpc_filter, len(mpc_data))

    Parameters
    ----------
    filter:
        Filter function which defines which group to select. The filter function must
        accept 3 parameters, `peri_dist, eccentricity, h_mag`, and return a bool.
        See `neospy.population.definitions` for a collection of filter functions which
        are used to generation model populations.
    n_samples:
        How many samples to generate from the provided dataset.
    bandwidth:
        The bandwidth of the gaussian KDEs in fractions of standard deviation of the
        orbital parameters.
    seed:
        A random seed to use to generate the values, the same seed will result in the
        same draws.
    """
    columns = ["incl", "peri_arg", "peri_dist", "peri_time", "ecc", "lon_node"]
    mpc_orbits = mpc_known_orbit_filtered(filt)
    epoch = mpc_orbits.epoch.median()
    for scale in [1.25, 3, 10]:
        extra_samples = int(np.ceil(n_samples * scale))
        data = pd.DataFrame(
            sample_values(
                np.transpose(mpc_orbits[columns]),
                extra_samples,
                bandwidth,
                seed=seed,
            ).T,
            columns=columns,
        )
        h_mag = np.zeros(len(data))
        filt_bool = filt(data.peri_dist, data.ecc, h_mag)
        filt_bool = filt_bool & (data.incl >= 0.0)
        data['epoch'] = epoch
        if len(data[filt_bool]) >= n_samples:
            data["peri_arg"] %= 360
            data["lon_node"] %= 360
            return data[filt_bool][:n_samples]
    raise ValueError(
        "Filter was too aggressive with provided data, failed to find valid sample."
    )


def sample_beaming(collection, orbit_filter, n_samples, bandwidth=0.05, seed=42):
    """
    Create a Beaming Parameter sample.

    This samples the PDS dataset with the specified orbital elements filter, returning
    new beaming parameters. This assumes that beaming is uncorrelated with other values
    within the specified orbit filter.

    Parameters
    ----------
    collection:
        A PDS data filename, such as ``"neowise_neos.xml"``
    orbit_filter:
        Filter function which defines which group to select. The filter function must
        accept 3 parameters, `peri_dist, eccentricity, h_mag`, and return a bool.
        See `neospy.population.definitions` for a collection of filter functions which
        are used to generation model populations.
    n_samples:
        How many samples to generate from the provided dataset.
    bandwidth:
        The bandwidth of the gaussian KDEs in fractions of standard deviation of the
        orbital parameters.
    seed:
        A random seed to use to generate the values, the same seed will result in the
        same draws.
    """
    beaming_data = pds_data_filtered(collection, orbit_filter, fit_code="B")
    for scale in [1.5, 3, 10]:
        extra_samples = int(np.ceil(n_samples * scale))

        beaming = sample_values(
            np.array(beaming_data["Beaming_param"])[np.newaxis, :],
            extra_samples,
            bandwidth=bandwidth,
            seed=seed,
        ).ravel()
        keep = (0 < beaming) & (beaming < np.pi)
        beaming = beaming[keep][:n_samples]
        if len(beaming) == n_samples:
            return pd.DataFrame(np.transpose([beaming]), columns=["Beaming_param"])
    raise ValueError(
        "Filter was too aggressive with provided data, failed to find valid sample."
    )


def sample_albedos(collection, orbit_filter, n_samples, bandwidth=0.05, seed=42):
    """
    Create a V and IR albedo sample using the PDS dataset.

    This samples the PDS dataset with the specified orbital elements filter, returning
    new albedo parameters. This assumes that albedo is uncorrelated with other values
    within the specified orbit filter.

    This assumes that the objects with `V` fits in the PDS data are I.I.D. (in other
    words a representative sample of the group specified by the filter).

    Objects which have an `I` fit do not appear to be an IID sample, but those fits
    are used to estimate the ratio of V_Albedo / IR_albedo, and the V_albedo calculated
    from the strictly `V` fits are then used to create an IR albedo.

    Parameters
    ----------
    collection:
        A PDS data filename, such as ``"neowise_neos.xml"``
    orbit_filter:
        Filter function which defines which group to select. The filter function must
        accept 3 parameters, `peri_dist, eccentricity, h_mag`, and return a bool.
        See `neospy.population.definitions` for a collection of filter functions which
        are used to generation model populations.
    n_samples:
        How many samples to generate from the provided dataset.
    bandwidth:
        The bandwidth of the gaussian KDEs in fractions of standard deviation of the
        orbital parameters.
    seed:
        A random seed to use to generate the values, the same seed will result in the
        same draws.
    """
    albedo_ratio_data = pds_data_filtered(collection, orbit_filter, fit_code="VBI")
    vis_albedo_data = pds_data_filtered(collection, orbit_filter, fit_code="V")
    for scale in [1.5, 3, 10]:
        extra_samples = int(np.ceil(n_samples * scale))
        v_albedo = sample_values(
            np.array(vis_albedo_data["V_albedo"])[np.newaxis, :],
            extra_samples,
            bandwidth=bandwidth,
            seed=seed,
        ).ravel()
        ir_ratios = sample_values(
            np.transpose(albedo_ratio_data[["V_albedo", "IR_albedo"]]),
            extra_samples,
            bandwidth=bandwidth,
            seed=seed,
        )
        ir_ratios = np.array(ir_ratios[1] / ir_ratios[0]).ravel()
        ir_albedo = v_albedo * ir_ratios

        keep = (v_albedo > 0.015) & (v_albedo < 0.7)
        keep = keep & (ir_albedo > 0.015) & (ir_albedo < 0.7)
        v_albedo = v_albedo[keep][:n_samples]
        ir_albedo = ir_albedo[keep][:n_samples]
        if len(v_albedo) == n_samples:
            return pd.DataFrame(
                np.transpose([v_albedo, ir_albedo]), columns=["V_albedo", "IR_albedo"]
            )
    raise ValueError(
        "Filter was too aggressive with provided data, failed to find valid sample."
    )


def visibility_test(
    peri_dist: NDArray[np.floating],
    vis_albedo: NDArray[np.floating],
    beaming: NDArray[np.floating],
    diameter: NDArray[np.floating],
    wavelength: float = 8e-06,
    solar_elong: float = 120,
    min_flux: float = 50e-6,
) -> np.ndarray:
    """
    Given lists of physical parameters, determine if the object would ever be visible
    if it is located at the optimal geometry for observation

    This places the object at the specified solar elongation angle from the observer, at
    the perihelion distance from the sun.

    Parameters
    ----------
    peri_dist:
        Perihelion distance of the object in AU.
    vis_albedo:
        The visible albedo of the object.
    beaming:
        The objects beaming parameter.
    diameter:
        The diameter in km.
    wavelength:
        The wavelength at which to calculate the flux from the object using NEATM.
    solar_elong:
        The solar elongation at which to place the object for the test.
    min_flux:
        The minimum flux to determine if the object could possibly be visible.
    """

    visible = []
    lon = (
        np.pi
        - np.arcsin(np.sin(np.radians(solar_elong)) / peri_dist)
        - np.radians(solar_elong)
    )
    obj2sun = -Vector(
        np.transpose(
            [peri_dist * np.cos(lon), peri_dist * np.sin(lon), np.zeros_like(lon)]
        )
    )

    obj2sc = Vector(
        np.transpose(
            [1 - peri_dist * np.cos(lon), -peri_dist * np.sin(lon), np.zeros_like(lon)]
        )
    )
    for o2s, o2o, alb, beam, r in zip(  # type: ignore
        obj2sun,
        obj2sc,
        vis_albedo,
        beaming,
        diameter,
    ):
        flux_ujy = neatm(
            o2s,
            o2o,
            geom_albedo=alb,
            G=0.15,
            beaming=beam,
            emissivity=0.9,
            diameter=r,
            wavelength=wavelength,
        )
        visible.append(flux_ujy > min_flux)

    return np.array(visible)
