from __future__ import annotations
import os
import pandas as pd  # type: ignore
import logging
from functools import lru_cache
import numpy as np

from .mpc import mpc_known_orbit_filtered, normalize_names
from .time import Time
from .data import cached_zip_download

logger = logging.getLogger(__name__)


@lru_cache()
def fetch_pds_data(filepath, url=None):
    """
    Fetch PDS data from the specified URL, if no URL is provided, default to the neowise
    diameter albedo v2.0 dataset.
    """
    # PDS4_tools appears to have some custom module finder which fails unless the neospy
    # is fully initialized.
    import pds4_tools  # type: ignore

    if url is None:
        url = (
            "https://sbnarchive.psi.edu/pds4/non_mission"
            "/neowise_diameters_albedos_V2_0.zip"
        )
    path = cached_zip_download(url)
    data = pds4_tools.read(os.path.join(path, filepath), quiet=True)
    return pd.DataFrame(data[0].data)


def pds_data_filtered(
    collection, orbit_filter, fit_code="DVBI", max_jd=Time.from_ymd(2010, 8, 6.5).jd
):
    """
    Given a file in the PDS dataset, an orbital parameter filter, and optionally a
    fit code (or partial fit code), return a DataFrame of all of the valid data
    collected inside of the PDS file specified.

    .. testcode::
        :skipif: True

        pds_data_filtered("neowise_mainbelt.xml",
                          complete_mba_center_filter,
                          fit_code='B')

    Parameters
    ----------
    collection:
        A PDS data filename, such as ``"neowise_neos.xml"``
    filter:
        Filter function which defines which group to select. The filter function must
        accept 3 parameters, `peri_dist, eccentricity, h_mag`, and return a bool.
        See `neospy.population.definitions` for a collection of filter functions which
        are used to generation model populations.
    fit_code:
        A string which will be used to subselect the fit code in the PDS dataset. This
        can be "DVBI" or any subset.
    max_jd:
        The maximum JD time to collect from the PDS data, defaults to the end of the
        cryo mission.
    """
    fit_code = "" if fit_code is None else fit_code

    if orbit_filter is not None:
        # Names of the objects in the MPC database which obey this filter
        names = set(mpc_known_orbit_filtered(orbit_filter).desig)

    dataset = normalize_names(fetch_pds_data("data/" + collection))
    dataset.fillna(value=np.nan)

    filt = dataset.Mean_JD < max_jd
    filt &= [fit_code in x for x in dataset.Fit_code]
    filt &= (0.015 <= dataset.V_albedo) & (dataset.V_albedo <= 0.6)
    filt &= (0.015 <= dataset.IR_albedo) & (dataset.IR_albedo <= 0.9)
    filt &= (0.6 <= dataset.Beaming_param) & (dataset.Beaming_param <= np.pi)
    dataset = dataset[filt]
    dataset = dataset.sort_values("Diameter_err", ascending=True).drop_duplicates(
        ["MPC_packed_name"]
    )
    if orbit_filter is not None:
        keep = [row in names for row in dataset["MPC_packed_name"]]
        return dataset[keep]
    return dataset
