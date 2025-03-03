"""
PTF Related Functions and Data.
"""

import os
import logging
from functools import lru_cache
from collections import defaultdict
from astropy.io import fits

from .cache import download_file, cache_path
from .fov import PtfCcd, PtfField, FOVList
from .time import Time
from .tap import query_tap
from .mpc import find_obs_code
from .vector import Vector, State
from . import spice


__all__ = ["fetch_fovs", "fetch_frame"]


logger = logging.getLogger(__name__)


@lru_cache(maxsize=2)
def fetch_fovs(year: int):
    """
    Load all FOVs taken during the specified mission year of PTF.

    This will download and cache all FOV information for the given year from IRSA.

    This can take about 20 minutes per year of survey, each year is 2-3 GB of data.

    Parameters
    ----------
    year :
        Which year of PTF.
    """
    year = int(year)
    cache_dir = cache_path()
    dir_path = os.path.join(cache_dir, "fovs")
    filename = os.path.join(dir_path, f"ptf_fields_{year}.bin")

    if not os.path.isdir(dir_path):
        os.makedirs(dir_path)
    if os.path.isfile(filename):
        return FOVList.load(filename)

    table = "ptf.ptf_procimg"
    cols = [
        "fieldid",
        "obsjd",
        "ccdid",
        "filter",
        "pfilename",
        "infobits",
        "seeing",
        "ra1",
        "dec1",
        "ra2",
        "dec2",
        "ra3",
        "dec3",
        "ra4",
        "dec4",
    ]
    jd_start = Time.from_ymd(year, 1, 1).jd
    jd_end = Time.from_ymd(year + 1, 1, 1).jd

    irsa_query = query_tap(
        f"SELECT {', '.join(cols)} FROM {table} "
        f"WHERE obsjd between {jd_start} and {jd_end}",
        verbose=True,
    )

    # Exposures are 30 seconds
    jds = [Time(x, scaling="utc") for x in irsa_query["obsjd"]]
    obs_info = find_obs_code("ZTF")

    # PTF fields are made up of up to 11 individual CCDs, here we first construct
    # the individual CCD information.
    fovs = []
    for jd, row in zip(jds, irsa_query.itertuples()):
        corners = []
        for i in range(4):
            ra = getattr(row, f"ra{i + 1}")
            dec = getattr(row, f"dec{i + 1}")
            corners.append(Vector.from_ra_dec(ra, dec))
        observer = spice.earth_pos_to_ecliptic(jd, *obs_info[:-1])
        observer = State("PTF", observer.jd, observer.pos, observer.vel)

        try:
            fov = PtfCcd(
                corners,
                observer,
                row.fieldid,
                row.ccdid,
                row.filter,
                row.pfilename,
                row.infobits,
                row.seeing,
            )
        except Exception:
            print(
                corners,
                observer,
                row.fieldid,
                row.ccdid,
                row.filter,
                row.pfilename,
                row.infobits,
                row.seeing,
            )
        fovs.append(fov)

    # Now group the quad information into full 64 size Fields
    grouped = defaultdict(list)
    for fov in fovs:
        grouped[fov.field].append(fov)

    # Sort the fovs by ccdid and make PTF Fields
    final_fovs = []
    for value in grouped.values():
        value = sorted(value, key=lambda x: (x.ccdid))
        fov = PtfField(value)
        final_fovs.append(fov)

    # finally save and return the result
    fov_list = FOVList(final_fovs)
    fov_list.save(filename)
    return fov_list


def fetch_frame(
    fov: PtfCcd,
    force_download=False,
    retry=2,
):
    """
    Given a PTF FOV, return the FITs file associated with it.

    This downloads the fits file into the cache.

    Parameters
    ----------
    fov :
        A single CCD FOV.
    force_download :
        Optionally force a re-download if the file already exists in the cache.
    """
    ptf_base = (
        f"https://irsa.ipac.caltech.edu/ibe/data/ptf/images/level1/{fov.filename}"
    )
    file = download_file(
        ptf_base, force_download=force_download, auto_zip=True, subfolder="ptf_frames"
    )

    try:
        return fits.open(file)[0]
    except OSError as exc:
        if retry == 0:
            raise ValueError("Failed to fetch PTF frame.") from exc
        logger.info("PTF file appears corrupted, attempting to fetch again.")
        os.remove(file)
        return fetch_frame(fov, force_download, retry - 1)
