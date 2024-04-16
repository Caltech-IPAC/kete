"""Thermal calculations."""

from __future__ import annotations
import numpy as np
from ..vector import Vector

# pylint: disable=no-name-in-module
from .. import _rust  # type: ignore

# pylint: disable=import-error
from .._rust import (  # type: ignore
    neatm_facet_temps,
    frm_facet_temps,
    FrmParams,
    NeatmParams,
    ModelResults,
)


__all__ = [
    "subsolar_temp",
    "neatm_facet_temps",
    "frm_facet_temps",
    "black_body_radiation",
    "neatm",
    "frm",
    "FrmParams",
    "NeatmParams",
    "ModelResults",
]


def black_body_radiation(temperature: Vector, wavelength: float) -> np.ndarray:
    """
    Compute the black body flux at the specified temperature and wavelength.

    Flux is in units Janskys / steradian.
    This can be converted to being per unit wavelength by multiplying the result of this
    by `c / wavelength**2`.

    This is limited by the size of float 32s, numbers exceeding 1e-300 will be rounded
    to 0. This should not raise any warnings, and should be relatively numerically
    stable and correct over most domains of interest.

    Temperature may be a vector.

    Parameters
    ----------
    temperature :
        In units of Kelvin.
    wavelength :
        In units of nanometers.

    Returns
    -------
    Iterable
        Janskys per steradian
    """
    temperature = np.array(temperature, copy=False)
    if temperature.ndim == 0:
        return _rust.black_body_flux([temperature], wavelength)[0]
    return _rust.black_body_flux(list(temperature), wavelength)


def subsolar_temp(
    obj2sun: Vector,
    geom_albedo: float,
    G: float,
    emissivity: float,
    beaming: float = np.pi,
) -> float:
    """
    Compute the temperature at the sub-solar point on the object.

    Parameters
    ----------
    obj2sun :
        The vector from the object to the sun, in units of AU.
    geom_albedo :
        The geometric albedo of the object.
    G :
        The `G` parameter of the object.
    emissivity :
        The emissivity of the object.
    beaming :
        The beaming parameter of the object, for the FRM model this is set to pi.

    Returns
    -------
    float
        Temperature at the sub-solar point on the object in Kelvin.
    """
    obj2sun = Vector(obj2sun)
    return _rust.sub_solar_temperature(obj2sun, geom_albedo, G, emissivity, beaming)


def neatm(
    sun2obj: Vector,
    sun2obs: Vector,
    geom_albedo: float,
    G: float,
    beaming: float,
    emissivity: float,
    diameter: float,
    wavelength: float = 23.68e3,
) -> float:
    """
    Calculate the flux in Janskys from an object using the NEATM thermal model.

    See :doc:`../auto_examples/plot_thermal_model`

    Parameters
    ----------
    sun2obj :
        A vector-like object containing the X/Y/Z coordinates pointing from the sun
        to the object in units of AU.
    sun2obs :
        A vector-like object containing the X/Y/Z coordinates pointing from the sun
        to the observer in units of AU.
    geom_albedo :
        The albedo of the object.
    G :
        The G parameter of the object
    beaming :
        The beaming parameter.
    emissivity :
        The emissivity of the object.
    diameter :
        The diameter of the object in km.
    wavelength :
        The wavelength of the object in nanometers.
    geometry :
        The geometry which defines the surface normals of the object.

    Returns
    -------
    float
        Flux in units of Jy.
    """
    sun2obj = Vector(sun2obj)
    sun2obs = Vector(sun2obs)
    return _rust.neatm_thermal(
        sun2obj,
        sun2obs,
        geom_albedo,
        G,
        beaming,
        diameter,
        wavelength,
        emissivity,
    )


def frm(
    sun2obj: Vector,
    sun2obs: Vector,
    geom_albedo: float,
    G: float,
    emissivity: float,
    diameter: float,
    wavelength: float = 23.68e3,
) -> float:
    """
    Calculate the flux from an object using the FRM thermal model in Janskys.

    See :doc:`../auto_examples/plot_thermal_model`

    Parameters
    ----------
    sun2obj :
        A vector-like object containing the X/Y/Z coordinates pointing from the sun
        to the object in units of AU.
    sun2obs :
        A vector-like object containing the X/Y/Z coordinates pointing from the sun
        to the observer in units of AU.
    geom_albedo :
        The albedo of the object.
    G :
        The G parameter of the object
    emissivity :
        The emissivity of the object.
    diameter :
        The diameter of the object in km.
    wavelength :
        The wavelength of the object in nanometers.
    geometry :
        The geometry which defines the surface normals of the object.

    Returns
    -------
    float
        Flux in units of Jy.
    """
    sun2obj = Vector(sun2obj)
    sun2obs = Vector(sun2obs)
    return _rust.frm_thermal(
        sun2obj,
        sun2obs,
        geom_albedo,
        G,
        diameter,
        wavelength,
        emissivity,
    )
