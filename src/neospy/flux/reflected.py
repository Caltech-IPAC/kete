import numpy as np
from ..vector import Vector

# pylint: disable=no-name-in-module
from .. import _core  # type: ignore


def hg_phase_curve_correction(G: float, phase: float) -> float:
    """
    This computes the phase curve correction in the IAU format.

    Parameters
    ----------
    G :
        The G parameter of the object.
    phase :
        The angular separation between the sun and observer as viewed from the object.

    Returns
    -------
    float
        The phase curve correction.
    """

    def f(A, B, C, phase):
        phase = np.radians(phase)
        s_phase = np.sin(phase)
        theta_L = np.exp(-A * np.tan(0.5 * phase) ** B)
        theta_S = 1.0 - C * s_phase / (0.119 + 1.341 * s_phase - 0.754 * s_phase**2.0)
        W = np.exp(-90.56 * np.tan(0.5 * phase) ** 2.0)
        return W * theta_S + (1.0 - W) * theta_L

    return (1.0 - G) * f(3.332, 0.631, 0.986, phase) + G * f(1.862, 1.218, 0.238, phase)


def hg_apparent_flux(
    sun2obj: Vector,
    sun2obs: Vector,
    H: float,
    G: float,
    c_hg: float,
    diameter: float,
    wavelength: float,
    geom_albedo: float,
) -> float:
    """
    Calculate the reflected flux from an object using the HG asteroid model.

    This assumes that the object is an ideal disk facing the sun and applies the IAU
    correction curve to the reflected light, returning units of Jy per unit frequency.

    This treats the sun as a flat black body disk, which is a good approximation as long
    as the object is several solar radii away.

    Parameters
    ----------
    obj2sun :
        A vector-like object containing the X/Y/Z coordinates pointing from the sun
        to the object in units of AU.
    sun2obs :
        A vector-like object containing the X/Y/Z coordinates pointing from the sun
        to the observer in units of AU.
    H :
        H magnitude of the object in V band.
    G :
        The G parameter of the object.
    c_hg :
        The relationship constant between H, D, and pV for the bandpass
    diameter :
        The diameter of the object in km.
    wavelength :
        The wavelength of the object in nanometers.
    geom_albedo :
        The albedo of the object at the specified wavelength.

    Returns
    -------
    float
        Flux in Jy per unit frequency.
    """
    sun2obj = Vector(sun2obj)
    sun2obs = Vector(sun2obs)
    # Jy / steradian per unit freq
    return _core.hg_apparent_flux(
        sun2obj, sun2obs, H, G, c_hg, diameter, wavelength, geom_albedo
    )


def hg_absolute_to_apparent_mag(
    sun2obj: Vector, sun2obs: Vector, G: float, H: float
) -> float:
    """
    Compute the apparent magnitude of an object using the absolute magnitude H, G, and
    positional information.

    The HG IAU model is not technically defined above 120 degrees phase, however this
    will continue to return values fit to the model until 160 degrees. Phases larger
    than 160 degrees will return an apparent magnitude of infinity.

    Parameters
    ----------
    obj2sun :
        A vector-like object containing the X/Y/Z coordinates pointing from the sun
        to the object in units of AU.
    sun2obs :
        A vector-like object containing the X/Y/Z coordinates pointing from the sun
        to the observer in units of AU.
    G :
        The G parameter of the object.
    H :
        The absolute magnitude of the object.

    Returns
    -------
    float
        The apparent magnitude of the object.
    """
    sun2obj = Vector(sun2obj)
    sun2obs = Vector(sun2obs)

    return _core.hg_apparent_mag(sun2obj, sun2obs, H, G)
