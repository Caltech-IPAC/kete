"""
Field of view definitions, along with tests for checks to see if objects are within the
FOVs.
"""

from ._core import (
    NeosCmos,
    NeosVisit,
    WiseCmos,
    ZtfCcdQuad,
    ZtfField,
    RectangleFOV,
    FOVList,
    fov_static_check,
    fov_state_check,
    fov_spk_check as _fov_spk_check,
)
from .vector import State
from .spice import SpiceKernels


__all__ = [
    "NeosCmos",
    "NeosVisit",
    "WiseCmos",
    "ZtfCcdQuad",
    "ZtfField",
    "RectangleFOV",
    "FOVList",
    "fov_static_check",
    "fov_state_check",
    "fov_spice_check",
]


def fov_spice_check(desigs: list[str], fovs) -> list[State]:
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
    return _fov_spk_check(obj_ids, fovs)
