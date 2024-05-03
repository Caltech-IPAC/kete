"""
Definitions for population groups.

Group definitions are computed strictly from perihelion distance and eccentricity.
"""

from __future__ import annotations
import matplotlib.pyplot as plt  # type: ignore
import numpy as np
from ..conversion import compute_aphelion, compute_semi_major
from numpy.typing import NDArray

MAX_ECCENTRICITY = 0.995
"""Maximum allowed eccentricity to be considered non-hyperbolic."""


def mba_inner(
    peri_dist: NDArray[np.floating], ecc: NDArray[np.floating], *_
) -> np.ndarray:
    """
    Orbital element filter to select inner main belt objects.

    Returns `True` if the object is in the inner main belt.

    Parameters
    ----------
    peri_dist:
        Perihelion distance in units of AU.
    ecc:
        Eccentricity of the orbit.
    """
    a = compute_semi_major(peri_dist, ecc)
    return (peri_dist > 1.3) & (a <= 2.5) & (ecc < MAX_ECCENTRICITY) & (ecc >= 0.0)


def mba_middle(
    peri_dist: NDArray[np.floating], ecc: NDArray[np.floating], *_
) -> np.ndarray:
    """
    Orbital element filter to select middle main belt objects.

    Returns `True` if the object is in the middle main belt.

    Parameters
    ----------
    peri_dist:
        Perihelion distance in units of AU.
    ecc:
        Eccentricity of the orbit.
    """
    a = compute_semi_major(peri_dist, ecc)
    return (
        (peri_dist > 1.3)
        & (a > 2.5)
        & (a <= 2.82)
        & (ecc < MAX_ECCENTRICITY)
        & (ecc >= 0.0)
    )


def mba_outer(
    peri_dist: NDArray[np.floating], ecc: NDArray[np.floating], *_
) -> np.ndarray:
    """
    Orbital element filter to select outer main belt objects.

    Returns `True` if the object is in the outer main belt.

    Parameters
    ----------
    peri_dist:
        Perihelion distance in units of AU.
    ecc:
        Eccentricity of the orbit.
    """
    a = compute_semi_major(peri_dist, ecc)
    return (
        (peri_dist > 1.3)
        & (a > 2.82)
        & (a <= 3.6)
        & (ecc < MAX_ECCENTRICITY)
        & (ecc >= 0.0)
    )


def neo_atira(
    peri_dist: NDArray[np.floating], ecc: NDArray[np.floating], *_
) -> np.ndarray:
    """
    Orbital element filter to select Atira NEO objects.

    Returns `True` if the object is in the Atiras.

    Parameters
    ----------
    peri_dist:
        Perihelion distance in units of AU.
    ecc:
        Eccentricity of the orbit.
    """
    a = compute_semi_major(peri_dist, ecc)
    aphelion = compute_aphelion(peri_dist, ecc)
    return (
        (peri_dist <= 1.3)
        & (peri_dist >= 0.0)
        & (a < 1.0)
        & (aphelion < 0.983)
        & (ecc < MAX_ECCENTRICITY)
        & (ecc >= 0.0)
    )


def neo_aten(
    peri_dist: NDArray[np.floating], ecc: NDArray[np.floating], *_
) -> np.ndarray:
    """
    Orbital element filter to select Aten NEO objects.

    Returns `True` if the object is in the Atens.

    Parameters
    ----------
    peri_dist:
        Perihelion distance in units of AU.
    ecc:
        Eccentricity of the orbit.
    """
    a = compute_semi_major(peri_dist, ecc)
    aphelion = compute_aphelion(peri_dist, ecc)
    return (
        (peri_dist <= 1.3)
        & (peri_dist >= 0.0)
        & (a < 1.0)
        & (aphelion >= 0.983)
        & (ecc < MAX_ECCENTRICITY)
        & (ecc >= 0.0)
    )


def neo_apollo(
    peri_dist: NDArray[np.floating], ecc: NDArray[np.floating], *_
) -> np.ndarray:
    """
    Orbital element filter to select Apollo NEO objects.

    Returns `True` if the object is in the Apollo.

    Parameters
    ----------
    peri_dist:
        Perihelion distance in units of AU.
    ecc:
        Eccentricity of the orbit.
    """
    a = compute_semi_major(peri_dist, ecc)
    return (peri_dist <= 1.017) & (a >= 1.0) & (ecc < MAX_ECCENTRICITY) & (ecc >= 0.0)


def neo_amor(
    peri_dist: NDArray[np.floating], ecc: NDArray[np.floating], *_
) -> np.ndarray:
    """
    Orbital element filter to select Amor NEO objects.

    Returns `True` if the object is in the Amors.

    Parameters
    ----------
    peri_dist:
        Perihelion distance in units of AU.
    ecc:
        Eccentricity of the orbit.
    """
    a = compute_semi_major(peri_dist, ecc)
    return (
        (peri_dist <= 1.3)
        & (a > 1.0)
        & (peri_dist > 1.017)
        & (ecc < MAX_ECCENTRICITY)
        & (ecc >= 0.0)
    )


def jup_trojan(
    peri_dist: NDArray[np.floating], ecc: NDArray[np.floating], *_
) -> np.ndarray:
    """
    Orbital element filter to select objects which appear to be jupiter trojans.

    Returns `True` if the object is not one of the above.

    Parameters
    ----------
    peri_dist:
        Perihelion distance in units of AU.
    ecc:
        Eccentricity of the orbit.
    """
    a = compute_semi_major(peri_dist, ecc)
    return (a < 5.5) & (a > 4.6) & (ecc < 0.3)


def distant(
    peri_dist: NDArray[np.floating], ecc: NDArray[np.floating], *_
) -> np.ndarray:
    """
    Orbital element filter to select distant objects outside the main belt.

    Returns `True` if the object is not one of the above.

    Parameters
    ----------
    peri_dist:
        Perihelion distance in units of AU.
    ecc:
        Eccentricity of the orbit.
    """
    a = compute_semi_major(peri_dist, ecc)
    return (
        (peri_dist > 1.3)
        & (a > 3.6)
        & (ecc < MAX_ECCENTRICITY)
        & (ecc >= 0.0)
        & ~jup_trojan(peri_dist, ecc)
    )


# Functions below use the orbital elements definitions above.


def which_group(
    peri_dist: NDArray[np.floating], ecc: NDArray[np.floating], *_
) -> list[str]:
    """
    Return a label for the object orbital elements provided.

    Parameters
    ----------
    peri_dist:
        Perihelion distance in units of AU.
    ecc:
        Eccentricity of the orbit.
    """
    groups = []
    for peri, e in zip(np.atleast_1d(peri_dist), np.atleast_1d(ecc)):
        for func in [
            mba_inner,
            mba_middle,
            mba_outer,
            neo_apollo,
            neo_aten,
            neo_amor,
            neo_atira,
            distant,
            jup_trojan,
        ]:
            if func(peri, e):
                groups.append(func.__name__)
                break
        else:
            groups.append("Other")
    return groups


def neo_amor_complete(peri_dist, ecc, h_mag):
    """
    Observationally complete set of Amors.

    Returns `True` if the object is in the Amors and is a part of what is expected to
    be observationally complete.

    Parameters
    ----------
    peri_dist:
        Perihelion distance in units of AU.
    ecc:
        Eccentricity of the orbit.
    h_mag:
        The H magnitude of the object.
    """
    return neo_amor(peri_dist, ecc) & (h_mag < 20.0)


def neo_apollo_complete(peri_dist, ecc, h_mag):
    """
    Observationally complete set of Apollos.

    Returns `True` if the object is in the Apollos and is a part of what is expected to
    be observationally complete.

    Parameters
    ----------
    peri_dist:
        Perihelion distance in units of AU.
    ecc:
        Eccentricity of the orbit.
    h_mag:
        The H magnitude of the object.
    """
    return neo_apollo(peri_dist, ecc) & (h_mag < 20.0)


def neo_aten_complete(peri_dist, ecc, h_mag):
    """
    Observationally complete set of Atens.

    Returns `True` if the object is in the Atens and is a part of what is expected to
    be observationally complete.

    Parameters
    ----------
    peri_dist:
        Perihelion distance in units of AU.
    ecc:
        Eccentricity of the orbit.
    h_mag:
        The H magnitude of the object.
    """
    return neo_aten(peri_dist, ecc) & (h_mag < 20.0)


def mba(peri_dist, ecc, *_):
    """
    Orbital element filter to select all main belt objects.

    Returns `True` if the object is in the main belt.

    Parameters
    ----------
    peri_dist:
        Perihelion distance in units of AU.
    ecc:
        Eccentricity of the orbit.
    """
    return (
        mba_inner(peri_dist, ecc)
        | mba_middle(peri_dist, ecc)
        | mba_outer(peri_dist, ecc)
    )


def neo(peri_dist, ecc, *_):
    """
    Orbital element filter to select all NEO objects.

    Returns `True` if the object is an NEO.

    Parameters
    ----------
    peri_dist:
        Perihelion distance in units of AU.
    ecc:
        Eccentricity of the orbit.
    """
    return (
        neo_atira(peri_dist, ecc)
        | neo_aten(peri_dist, ecc)
        | neo_apollo(peri_dist, ecc)
        | neo_amor(peri_dist, ecc)
    )


def neo_complete(peri_dist, ecc, h_mag):
    """
    Orbital element filter to select all NEO objects.

    Returns `True` if the object is an NEO and is a part of what is expected to be
    observationally complete.

    Parameters
    ----------
    peri_dist:
        Perihelion distance in units of AU.
    ecc:
        Eccentricity of the orbit.
    h_mag:
        The H magnitude of the object.
    """
    return (
        neo_atira(peri_dist, ecc)
        | neo_aten(peri_dist, ecc)
        | neo_apollo(peri_dist, ecc)
        | neo_amor(peri_dist, ecc)
    ) & (h_mag < 20.0)


def nearly_neo_complete(peri_dist, ecc, h_mag):
    """
    Orbital element filter to select all objects which are either NEO or nearly NEO.

    Returns `True` if the object is ~NEO and is a part of what is expected to be
    observationally complete.

    This function exists to solve a KDE orbit sampling problem:
    A histogram of the perihelion distance of NEOs shows an increase in object count up
    to the NEO/Mar-crosser boundary, leading to a hard cut at around perihelion of 1.3
    The KDE estimator has a commonly known issue with dataset with hard cuts like this,
    and has a tendency to "round" the corners of the histogram near this cutoff. To
    solve this issue, we sample from a collection of data beyond the limit of our group
    of interest, and select only the objects which match the properties we want. What
    this in effect does is over-sample, then subselect, this moves the "rounding" issue
    to an area where we will always reject the data, thus resulting in clean cuts where
    we expect.

    Parameters
    ----------
    peri_dist:
        Perihelion distance in units of AU.
    ecc:
        Eccentricity of the orbit.
    h_mag:
        The H magnitude of the object.
    """
    return (peri_dist <= 1.4) & (ecc < MAX_ECCENTRICITY) & (ecc >= 0.0) & (h_mag < 20.0)


def plot_groups():
    """
    Plot orbital element groups defined above.


    .. plot::
        :context: close-figs

        import neospy
        neospy.population.definitions.plot_groups();


    """

    ecc = np.linspace(0.01, 1.1, 501)
    peri = np.linspace(0.01, 5, 501)

    eccs, peris = np.meshgrid(ecc, peri)

    for func in [
        mba_inner,
        mba_middle,
        mba_outer,
        neo_apollo,
        neo_aten,
        neo_amor,
        neo_atira,
        distant,
        jup_trojan,
    ]:
        plt.scatter(
            peris[func(peris, eccs)], eccs[func(peris, eccs)], s=1, label=func.__name__
        )
    plt.xlabel("Perihelion Dist (AU)")
    plt.ylabel("Eccentricity")
    plt.legend(loc=1)
