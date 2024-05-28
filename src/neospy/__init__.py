from .vector import (
    Vector,
    Frames,
    State,
    CometElements,
    SimultaneousStates,
)
from . import (
    constants,
    covariance,
    wise,
    neos,
    cache,
    population,
    flux,
    mpc,
    irsa,
    ztf,
    fov,
)
from .spice import SpiceKernels
from .propagation import (
    propagate_n_body,
    propagate_two_body,
    moid,
    spice_visible,
    state_visible,
)
from .time import Time
from .conversion import (
    compute_H,
    compute_albedo,
    compute_diameter,
    compute_semi_major,
    compute_aphelion,
    mag_to_flux,
    flux_to_mag,
)
from .horizons import HorizonsProperties

import logging


__all__ = [
    "cache",
    "constants",
    "CometElements",
    "covariance",
    "irsa",
    "Frames",
    "moid",
    "Vector",
    "State",
    "population",
    "Time",
    "set_logging",
    "SpiceKernels",
    "flux",
    "fov",
    "wise",
    "neos",
    "mpc",
    "SimultaneousStates",
    "propagate_n_body",
    "propagate_two_body",
    "compute_H",
    "compute_albedo",
    "compute_diameter",
    "compute_semi_major",
    "compute_aphelion",
    "spice_visible",
    "state_visible",
    "mag_to_flux",
    "flux_to_mag",
    "Vector",
    "HorizonsProperties",
    "ztf",
]


def set_logging(level=logging.INFO, fmt="%(asctime)s - %(message)s"):
    """
    Output logging information to the console.

    Parameters
    ----------
    level:
        The logging level to output, if this is set to 0 logging is disabled.
    fmt:
        Format of the logging messages, see the ``logging`` package for format string
        details. Here is a more verbose output example:
        "%(asctime)s %(name)s:%(lineno)s - %(message)s"
    """
    logger = logging.getLogger(__name__)
    logger.setLevel(level)

    # If there is already a handler in the logger, dont add another
    logger.handlers.clear()

    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(logging.Formatter(fmt))
    logger.addHandler(ch)
    return logger


set_logging()
