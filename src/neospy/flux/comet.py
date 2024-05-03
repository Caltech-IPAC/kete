from ..vector import Vector

# pylint: disable=no-name-in-module
from .. import _core  # type: ignore


__all__ = [
    "comet_apparent_mags",
]


def comet_apparent_mags(
    sun2obj: Vector,
    sun2obs: Vector,
    m1: float,
    k1: float,
    m2: float,
    k2: float,
    phase_mag_slope=0.035,
) -> list[float]:
    """
    Given the M1/K1 and M2/K2 values, compute the apparent Comet visible magnitudes.

    The model for apparent magnitudes are:

    m1 + k1 * log10(sun2obj.r) + 5.0 * log10(obj2obs.r) + phase_mag_slope * phase
    m2 + k2 * log10(sun2obj.r) + 5.0 * log10(obj2obs.r) + phase_mag_slope * phase

    Where m1/k1 are related to total magnitudes and m2/k2 are nucleus magnitudes.

    This model is based off of these:
    https://ssd.jpl.nasa.gov/horizons/manual.html#obsquan  (see section 9)
    https://en.wikipedia.org/wiki/Absolute_magnitude#Cometary_magnitudes

    Note that the above model does not include a 2.5x term attached to the K1/2 terms
    which are present in the wikipedia definition, this matches the definitions used by
    JPL Horizons.

    This does a best effort to compute both magnitudes, if any values are missing this
    will return None in the respective calculation.

    Parameters
    ----------
    sun2obj :
        A vector-like object containing the X/Y/Z coordinates pointing to the object
        from the sun in units of AU.
    sun2obs :
        A vector-like object containing the X/Y/Z coordinates pointing from the sun
        to the observer in units of AU.
    m1 :
        Total Absolute magnitude of the comet.
    k1 :
        Total Absolute magnitude slope as a function of helio distance of the comet.
    m2 :
        Absolute magnitude of the nucleus of the comet.
    k2 :
        Absolute magnitude slope of the nucleus as a function of helio distance of the
        comet.
    phase_mag_slope :
        Magnitude variation of the comet as a function of observing phase, units are
        Mag/Deg of phase, this defaults to 0.035 Mag/Deg which is a reasonable
        assumption.

    Returns
    -------
    float
        (Total apparent magnitude, Magnitude of the nucleus)
    """
    if m1 is not None and k1 is not None:
        mk1 = [m1, k1]
    else:
        mk1 = None
    if m2 is not None and k2 is not None:
        mk2 = [m2, k2]
    else:
        mk2 = None

    return _core.comet_apparent_mags(sun2obj, sun2obs, mk1, mk2, phase_mag_slope)
