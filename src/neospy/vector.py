"""
Representation of States and vectors.
"""

from __future__ import annotations
from . import conversion

# pylint: disable=import-error
from ._rust import Vector, Frames, State, CometElements  # type: ignore


__all__ = ["Frames", "Vector", "State", "CometElements"]

Vector.dec_dms = property(
    fget=lambda self: conversion.dec_degrees_to_dms(self.dec),
    doc="The declination, in degrees-arcminutes-arcseconds string format.",
)

Vector.ra_hms = property(
    fget=lambda self: conversion.ra_degrees_to_hms(self.ra),
    doc="The right ascension, in hours-minutes-seconds string format.",
)
