"""
Useful conversion functions between various physical values or representations.
"""

from __future__ import annotations
import numpy as np
import logging
from typing import Union
from numpy.typing import NDArray
from . import _core  # pylint: disable=no-name-in-module
from . import constants


logger = logging.getLogger(__name__)


def compute_H(diameter: NDArray, albedo: NDArray, c_hg=constants.C_V) -> np.ndarray:
    """
    Compute the H magnitude of the object given the diameter in km and the albedo.

    Parameters
    ----------
    diameter:
        The diameter of the object in Km.
    albedo:
        The albedo of the object.
    """
    diameter = np.array(diameter, dtype=float, copy=False)
    return -5.0 * np.log10(diameter * np.sqrt(albedo) / c_hg)


def compute_albedo(diameter: NDArray, H: NDArray, c_hg=constants.C_V) -> np.ndarray:
    """
    Compute the albedo of the object given the diameter in km and H magnitude.

    Parameters
    ----------
    diameter:
        The diameter of the object in Km.
    H:
        The H magnitude of the object.
    """
    return np.clip((c_hg * 10.0 ** (-0.2 * H) / diameter) ** 2.0, 0, 1)


def compute_diameter(albedo: NDArray, H: NDArray, c_hg=constants.C_V) -> np.ndarray:
    """
    Compute the diameter of the object in km given the albedo and H magnitude.

    Parameters
    ----------
    albedo:
        The albedo of the object.
    H:
        The H magnitude of the object.
    """
    return (c_hg / np.sqrt(albedo)) * 10.0 ** (-0.2 * H)


def compute_semi_major(peri_dist: NDArray, ecc: NDArray) -> np.ndarray:
    """
    Calculate semi major axis, returning nan when it cannot be computed.

    Parameters
    ----------
    peri_dist:
        Perihelion distance in units of AU.
    ecc:
        Eccentricity of the orbit.
    """
    return np.divide(
        peri_dist,
        (1.0 - ecc),
        out=np.full_like(ecc, np.nan, dtype=float),
        where=ecc < 1.0,
    )


def compute_aphelion(peri_dist: NDArray, ecc: NDArray) -> np.ndarray:
    """
    Calculate aphelion, returning nan when it cannot be computed.

    Parameters
    ----------
    peri_dist :
        Perihelion distance in units of AU.
    ecc :
        Eccentricity of the orbit.
    """
    return np.divide(
        peri_dist * (1.0 + ecc),
        (1.0 - ecc),
        out=np.full_like(ecc, np.nan, dtype=float),
        where=ecc < 1.0,
    )


def compute_earth_radius(geodetic_latitude: float) -> float:
    """
    Compute the effective Earth's radii at the specified geodetic latitude assuming that
    the Earth is an oblate spheroid.

    Returns the radii in units of AU.

    Parameters
    ----------
    geodetic_latitude:
        The geodetic latitude in degrees.
    """
    # https://en.wikipedia.org/wiki/Earth_radius#Geocentric_radius
    geodetic_latitude = np.radians(geodetic_latitude)
    a = constants.EARTH_MAJOR_AXIS_M
    b = constants.EARTH_MINOR_AXIS_M
    a_cos = a * np.cos(geodetic_latitude)
    b_sin = b * np.sin(geodetic_latitude)
    return (
        np.sqrt(((a * a_cos) ** 2 + (b * b_sin) ** 2) / (a_cos**2 + b_sin**2))
        / constants.AU_M
    )


def compute_eccentric_anomaly(
    eccentricity: NDArray, mean_anomaly: NDArray, peri_dist: NDArray
) -> np.ndarray:
    """
    Solve Kepler's equation for the eccentric anomaly.

    Parameters
    ----------
    eccentricity:
        The eccentricity of the orbit, greater than or equal to 0.
    mean_anomaly:
        The mean anomaly of the orbit in degrees.
    peri_dist:
        The perihelion distance, only required for parabolic objects. (Units of AU).
    """
    eccentricity = np.atleast_1d(eccentricity)
    mean_anomaly = np.radians(np.atleast_1d(mean_anomaly))
    peri_dist = np.atleast_1d(peri_dist)
    return np.degrees(
        _core.compute_eccentric_anomaly(eccentricity, mean_anomaly, peri_dist),
        dtype=float,
    )


def dec_degrees_to_dms(dec: NDArray) -> Union[str, list[str]]:
    """
    Convert a declination in degrees to a "degrees arcminutes arcseconds" string.

    Parameters
    ----------
    dec:
        Declination in decimal degrees.
    """
    if np.any(np.abs(dec) > 90):
        raise ValueError("Dec must be between -90 and 90")

    dec = np.array(dec)

    def _val(dec):
        degree = np.fix(dec)
        arc_min = int((dec - degree) * 60)
        arc_sec = (((dec - degree) * 60) - arc_min) * 60
        return f"{degree:+03.0f} {abs(arc_min):02d} {abs(arc_sec):05.2f}"

    res = np.vectorize(_val)(dec)
    if res.ndim == 0:
        return str(res)
    if res.ndim == 1:
        return list(res)
    return res


def dec_dms_to_degrees(dec: str) -> float:
    """
    Convert a declination from "degrees arcminutes arcseconds" string to degrees.

    This must be formatted with a space between the terms.

    Parameters
    ----------
    dec:
        Declination in degrees-arcminutes-arcseconds.
    """
    dec = dec.rstrip()

    if dec[0] not in "+- ":
        raise ValueError("Dec must have a sign as the first character")

    sign = -1 if dec[0] == "-" else 1

    dms = dec.split()
    if len(dms) != 3:
        raise ValueError("Dec must be formatted with space-separation")

    return sign * (abs(int(dms[0])) + int(dms[1]) / 60.0 + float(dms[2]) / 3600.0)


def ra_degrees_to_hms(ra: NDArray) -> Union[str, list[str]]:
    """
    Convert a Right Ascension in decimal degrees to an "hours minutes seconds" string.

    Parameters
    ----------
    ra:
        Right Ascension in decimal degrees.
    """
    ra = np.array(ra) % 360

    def _val(ra):
        ra_time = ra / 15.0
        hours = int(ra_time)
        minutes = int((ra_time - hours) * 60)
        seconds = (((ra_time - hours) * 60) - minutes) * 60
        return f"{hours:02d} {minutes:02d} {seconds:06.3f}"

    res = np.vectorize(_val)(ra)
    if res.ndim == 0:
        return str(res)
    if res.ndim == 1:
        return list(res)
    return res


def ra_hms_to_degrees(ra: str) -> float:
    """
    Convert a right ascension from "hours minutes seconds" string to degrees.

    This must be formatted with a space between the terms.

    Parameters
    ----------
    ra:
        Right ascension in hours-minutes-seconds.
    """

    hms = ra.split()
    if len(hms) != 3:
        raise ValueError("RA must be formatted with space-separation")

    return (int(hms[0]) + int(hms[1]) / 60.0 + float(hms[2]) / 3600.0) * 15


def flux_to_mag(flux: float, zero_point=3631) -> float:
    """
    Convert flux in Jy to AB Magnitude, assuming it is a single frequency source.
    Note that this assumes that the band is infinitely narrow.

    See: https://en.wikipedia.org/wiki/AB_magnitude

    Parameters
    ----------
    flux:
        Flux in Jy.
    zero_point:
        Flux in Jy where the magnitude is zero.
    """
    if flux < 1e-14:
        return np.inf
    return -2.5 * np.log10(flux / zero_point)


def mag_to_flux(mag: float, zero_point=3631) -> float:
    """
    Convert AB Magnitude to flux in Jy, assuming it is a single frequency source.
    Note that this assumes that the band is infinitely narrow.

    See: https://en.wikipedia.org/wiki/AB_magnitude

    Parameters
    ----------
    flux:
        Flux in Jy.
    zero_point:
        Flux in Jy where the magnitude is zero.
    """
    return 10 ** (mag / -2.5) * zero_point
