"""
Propagation of objects using orbital mechanics, this includes a simplified 2 body model
as well as a N body model which includes some general relativistic effects.
"""

from __future__ import annotations
import logging
from typing import Optional
from scipy import optimize  # type: ignore
import numpy as np

from .spice import SpiceKernels
from .vector import Vector, State
from .fov import FOVList
from . import _core  # pylint: disable=no-name-in-module

logger = logging.getLogger(__name__)


def propagate_n_body(
    states: list[State],
    jd: float,
    include_asteroids: bool = False,
    a_terms: Optional[list[Optional[tuple[float, float, float, bool]]]] = None,
    suppress_errors: bool = True,
) -> list[State]:
    """
    Propagate the provided :class:`~neospy.State` using N body mechanics to the
    specified times, no approximations are made, this can be very CPU intensive.

    This does not compute light delay, however it does include corrections for general
    relativity due to the Sun.

    Parameters
    ----------
    states:
        The initial states, this is a list of multiple State objects.
    jd:
        A JD to propagate the initial states to.
    include_asteroids:
        If this is true, the computation will include the largest 5 asteroids.
        The asteroids are: Ceres, Pallas, Interamnia, Hygiea, and Vesta.
    a_terms:
        A list of non-gravitational terms for each object. If provided, then every
        object must have a defined tuple containing (A_1, A_2, A_3, bool), where
        A_1 is the radial 1/r^2 correction, A_2 is the correction along the
        tangential direction of motion, and A_3 is the normal term. The bool defines
        if this object obeys the cometary force fall-off or simple 1/r^2.
        (True for comets, False for 1/r^2).
        A_1 is equivalent to the Beta term used in Cometary dust. These values are
        what are available on the JPL Horizons website for some objects.
    suppress_errors:
        If True, errors during propagation will return NaN for the relevant state
        vectors, but propagation will continue.

    Returns
    -------
    Iterable
        A :class:`~neospy.State` at the new time.
    """
    if a_terms is None:
        a_terms = [None for _ in range(len(states))]
    return _core.propagate_n_body_spk(
        states, jd, include_asteroids, a_terms, suppress_errors
    )


def propagate_two_body(
    states: list[State],
    jd: float,
    observer_pos: Optional[Vector] = None,
) -> list[State]:
    """
    Propagate the :class:`~neospy.State` for all the objects to the specified time.
    This assumes 2 body interactions.

    Parameters
    ----------
    states:
        The input vector state to propagate.
    jd:
        The desired time at which to estimate the objects' state.
    observer_pos:
        A vector of length 3 describing the position of an observer. If this is
        provided then the estimated states will be returned as a result of light
        propagation delay.

    Returns
    -------
    State
        Final state after propagating to the target time.
    """
    return _core.propagate_two_body(states, jd, observer_pos)


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
        The state of the object to calculate the MOID for, if this is not provided, then
        Earth is fetched from :class:`~neospy.spice.SpiceKernels` and is used in
        the calculation.
    """
    if other is None:
        other = SpiceKernels.state("Earth", state.jd)
    return _moid_single(state, other)


def state_visible(
    states: list[State], fovs: FOVList, dt: float = 3
) -> list[list[State]]:
    """
    Given states and field of view, return only the objects which are visible to the
    observer, adding a correction for optical light delay.

    Objects are propagated using 2 body physics to the time of the FOV if time steps are
    less than the specified `dt`.

    parameters
    ----------
    states:
        States which do not already have a specified FOV.
    fov:
        A field of view from which to subselect objects which are visible.
    dt:
        Length of time in days where 2-body mechanics is a good approximation.
    """
    return _core.fov_checks(states, fovs, dt)


def spice_visible(desigs: list[str], fovs) -> list[State]:
    """
    Given a list of object names and field of views, return only the objects which are
    visible to the observer, adding a correction for optical light delay.

    Objects are queried from the loaded SPK files. This does a best effort lookup and
    may silently not return states if an object is not loaded or doesn't have data for
    the specified epochs.

    parameters
    ----------
    desigs:
        Designations to lookup.
    fov:
        A list of field of views from which to subselect objects which are visible.
    """
    obj_ids = [SpiceKernels.name_lookup(n)[1] for n in desigs]
    return _core.fov_spk_checks(obj_ids, fovs)
