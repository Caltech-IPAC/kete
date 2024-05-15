from functools import lru_cache
import os
import numpy as np
from collections import defaultdict

from .data import cached_file_download, cache_path
from .fov import ZtfCcdQuad, ZtfField, FOVList
from .time import Time
from .irsa import query_irsa_tap
from .mpc import find_obs_code
from .vector import Vector
from .spice import SpiceKernels


__all__ = ["fetch_ZTF_file", "fetch_ZTF_fovs"]

SURVEY_START_JD = Time.from_ymd(2018, 3, 20).jd
"""First image in ZTF dataset."""

ZTF_IRSA_TABLES = {
    "ztf.ztf_current_meta_cal": "ZTF Calibration Metadata Table",
    "ztf.ztf_current_meta_deep": "ZTF Deep Reference Images",
    "ztf.ztf_current_meta_raw": "ZTF Raw Metadata Table",
    "ztf.ztf_current_meta_ref": "ZTF Reference (coadd) Images",
    "ztf.ztf_current_meta_sci": "ZTF Science Exposure Images",
    "ztf.ztf_current_path_cal": "ZTF Calibration Product Paths",
    "ztf.ztf_current_path_deep": "ZTF Deep Reference Product Paths",
    "ztf.ztf_current_path_raw": "ZTF Raw Product Paths",
    "ztf.ztf_current_path_ref": "ZTF Reference Product Paths",
    "ztf.ztf_current_path_sci": "ZTF Science Product Paths",
    "ztf_objects": "ZTF Objects",
    "ztf_objects_dr16": "ZTF Data Release 16 Objects",
    "ztf_objects_dr17": "ZTF Data Release 17 Objects",
    "ztf_objects_dr18": "ZTF Data Release 18 Objects",
    "ztf_objects_dr19": "ZTF Data Release 19 Objects",
    "ztf_objects_dr20": "ZTF Data Release 20 Objects",
}


@lru_cache(maxsize=3)
def fetch_ZTF_fovs(year: int):
    """
    Load all FOVs taken during the specified mission year of ZTF.

    This will download and cache all FOV information for the given year from IRSA.

    This can take about 20 minutes per year of survey, each year is 2-3 GB of data.

    Parameters
    ----------
    year :
        Which year of ZTF, 2018 through 2024.
    """
    year = int(year)
    if year not in [2018, 2019, 2020, 2021, 2022, 2023, 2024]:
        raise ValueError("Year must only be in the range 2018-2024")
    cache_dir = cache_path()
    dir_path = os.path.join(cache_dir, "fovs")
    filename = os.path.join(dir_path, f"ztf_fields_{year}.bin")

    if not os.path.isdir(dir_path):
        os.makedirs(dir_path)
    if os.path.isfile(filename):
        return FOVList.load(filename)

    table = "ztf.ztf_current_meta_sci"
    cols = [
        "field",
        "filefracday",
        "ccdid",
        "filtercode",
        "imgtypecode",
        "qid",
        "obsdate",
        "maglimit",
        "fid",
        "ra",
        "dec",
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

    irsa_query = query_irsa_tap(
        f"SELECT {', '.join(cols)} FROM {table} "
        f"WHERE obsjd between {jd_start} and {jd_end}",
        verbose=True,
    )

    # Exposures are 30 seconds, add 15 to select the midpoint of the observations
    jds_str = [x.split("+")[0] for x in irsa_query["obsdate"]]
    jds = np.array(Time(jds_str, "iso", "utc").jd) + 15 / 60 / 60 / 24

    obs_info = find_obs_code("ZTF")

    # ZTF fields are made up of up to 64 individual CCD quads, here we first construct
    # the individual CCD quad information.
    fovs = []
    for jd, row in zip(jds, irsa_query.itertuples()):
        corners = []
        for i in range(4):
            ra = row.__getattribute__(f"ra{i+1}")
            dec = row.__getattribute__(f"dec{i+1}")
            corners.append(Vector.from_ra_dec(ra, dec))
        observer = SpiceKernels.earth_pos_to_ecliptic(jd, *obs_info[:-1])

        fov = ZtfCcdQuad(
            corners,
            observer,
            row.field,
            row.filefracday,
            row.ccdid,
            row.filtercode,
            row.imgtypecode,
            row.qid,
            row.maglimit,
            row.fid,
        )
        fovs.append(fov)

    # Now group the quad information into full 64 size Fields
    grouped = defaultdict(list)
    for fov in fovs:
        key = (fov.filefracday, fov.fid, fov.filtercode)
        grouped[key].append(fov)

    # Sort the quads by ccdid and qid and make ZTF Fields
    final_fovs = []
    for value in grouped.values():
        value = sorted(value, key=lambda x: (x.ccdid, x.qid))
        fov = ZtfField(value)
        final_fovs.append(fov)

    # finally save and return the result
    fov_list = FOVList(final_fovs)
    fov_list.save(filename)
    return fov_list


def file_frac_day_split(filefracday):
    """
    Given the file frac day value from a field of view, return the year, month, and
    fraction of day.

    Note that the fraction of a day is in approximately JD UTC time, however there is
    about a 0.4 second offset in the ZTF time conversions for JD time.
    """
    filefracday = str(filefracday)
    year = int(filefracday[:4])
    month = int(filefracday[4:6])
    day = int(filefracday[6:8])
    frac_day = int(filefracday[8:])
    return (year, month, day, frac_day)


def fetch_ZTF_file(
    field,
    filefracday,
    filter_code,
    ccdid,
    image_type_code,
    qid,
    products="sci",
    im_type="sciimg.fits",
    force_download=False,
):
    """
    Fetch a ZTF file directly from the IPAC server, returning the path to where it was
    saved.
    """

    ztf_base = f"https://irsa.ipac.caltech.edu/ibe/data/ztf/products/{products}/"
    filefracday = str(filefracday)
    year = filefracday[:4]
    month_day = filefracday[4:8]
    frac_day = filefracday[8:]
    field = str(field).zfill(6)
    ccdid = str(ccdid).zfill(2)

    path = f"{year}/{month_day}/{frac_day}/"
    file = (
        f"ztf_{filefracday}_{field}_"
        f"{filter_code}_c{ccdid}_{image_type_code}_q{qid}_{im_type}"
    )

    url = ztf_base + path + file

    return cached_file_download(
        url, force_download=force_download, subfolder="ztf_frames"
    )
