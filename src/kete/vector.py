"""
Representation of States, Vectors, and coordinate Frames.
"""

from __future__ import annotations
from . import conversion

from ._core import (
    Vector,
    Frames,
    State,
    CometElements,
    SimultaneousStates,
    wgs_lat_lon_to_ecef,
    ecef_to_wgs_lat_lon,
)


__all__ = [
    "Frames",
    "Vector",
    "State",
    "CometElements",
    "SimultaneousStates",
    "wgs_lat_lon_to_ecef",
    "ecef_to_wgs_lat_lon",
]

Vector.dec_dms = property(
    fget=lambda self: conversion.dec_degrees_to_dms(self.dec),
    doc="The declination, in degrees-arcminutes-arcseconds string format.",
)

Vector.ra_hms = property(
    fget=lambda self: conversion.ra_degrees_to_hms(self.ra),
    doc="The right ascension, in hours-minutes-seconds string format.",
)
