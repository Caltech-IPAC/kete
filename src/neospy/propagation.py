"""
Propagation of objects using orbital mechanics, this includes a simplified 2 body model
as well as a N body model which includes some general relativistic effects.
"""

from __future__ import annotations
from typing import Optional
from scipy import optimize
import numpy as np

from .vector import State
from . import spice
from ._core import (
    NonGravModel,
    propagate_n_body,
    propagate_n_body_long,
    propagate_two_body,
)


__all__ = [
    "propagate_n_body",
    "propagate_n_body_long",
    "propagate_two_body",
    "NonGravModel",
    "moid",
]


def _moid_single(obj0: State, other: State):
    """
    Given the state of 2 objects, compute the MOID between them. This is used by the
    moid function below and is not intended to be used directly.
    """
    obj0_elem = obj0.elements
    obj1_elem = other.elements
    self_center = obj0_elem.peri_time
    self_period = obj0_elem.orbital_period
    other_center = obj1_elem.peri_time
    other_period = obj1_elem.orbital_period

    def _err(x):
        jd0, jd1 = x
        jd0 = jd0 * self_period / 4 + self_center
        jd1 = jd1 * other_period / 4 + other_center
        pos0 = propagate_two_body([obj0], jd0)[0].pos
        pos1 = propagate_two_body([other], jd1)[0].pos
        return np.linalg.norm(pos0 - pos1)

    soln = []
    soln.append(optimize.minimize(_err, [1, 1]).fun)
    soln.append(optimize.minimize(_err, [-1, -1]).fun)
    soln.append(optimize.minimize(_err, [-1, 1]).fun)
    soln.append(optimize.minimize(_err, [1, -1]).fun)
    return min(soln)


def moid(state: State, other: Optional[State] = None):
    """
    Compute the MOID between two objects assuming 2 body mechanics.

    If other is not provided, it is assumed to be Earth.

    Parameters
    ----------
    state:
        The state describing an object.
    other:
        The state of the object to calculate the MOID for, if this is not provided,
        then Earth is fetched from :mod:`~neospy.spice` and is used in the
        calculation.
    """
    if other is None:
        other = spice.get_state("Earth", state.jd)
    return _moid_single(state, other)
