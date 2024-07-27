"""
ZTF Related Functions and Data.
"""

from functools import lru_cache
import os
from collections import defaultdict
import matplotlib.pyplot as plt
import warnings

from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import (
    PowerDistStretch,
    AsymmetricPercentileInterval,
    ImageNormalize,
)


from .cache import cached_file_download, cache_path
from .fov import ZtfCcdQuad, ZtfField, FOVList
from .time import Time
from .irsa import query_irsa_tap
from .mpc import find_obs_code
from .vector import Vector, State
from . import spice


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

    # Exposures are 30 seconds
    jds = [Time.from_iso(x + ":00").jd for x in irsa_query["obsdate"]]
    obs_info = find_obs_code("ZTF")

    # ZTF fields are made up of up to 64 individual CCD quads, here we first construct
    # the individual CCD quad information.
    fovs = []
    for jd, row in zip(jds, irsa_query.itertuples()):
        corners = []
        for i in range(4):
            ra = getattr(row, f"ra{i + 1}")
            dec = getattr(row, f"dec{i + 1}")
            corners.append(Vector.from_ra_dec(ra, dec))
        observer = spice.earth_pos_to_ecliptic(jd, *obs_info[:-1])
        observer = State("ZTF", observer.jd, observer.pos, observer.vel)

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


def plot_ztf_fov(
    fov: ZtfCcdQuad,
    cmap="grey",
    products="sci",
    im_type="sciimg.fits",
    force_download=False,
):
    """
    Given a ztf FOV, plot the associated frame.

    This returns the associate WCS which is constructed.

    Parameters
    ----------
    fov :
        A single CCD Quad FOV.
    cmap :
        Colormap of the plot.
    products :
        Which data product to fetch.
    im_type :
        Image extension, this must match the products variable.
    force_download :
        Optionally force a re-download if the file already exists in the cache.
    """
    file = fetch_ZTF_file(
        fov.field,
        fov.filefracday,
        fov.filtercode,
        fov.ccdid,
        fov.imgtypecode,
        fov.qid,
        products,
        im_type,
        force_download,
    )

    frame = fetch_frame_from_fov(
        fov, products=products, im_type=im_type, force_download=force_download
    )

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        wcs = WCS(frame.header)

    if not plt.get_fignums():
        plt.figure(dpi=300, figsize=(6, 6), facecolor="w")

    ax = plt.subplot(projection=wcs)

    norm = ImageNormalize(
        frame.data,
        interval=AsymmetricPercentileInterval(10, 99.5),
        stretch=PowerDistStretch(0.25),
    )

    ax.imshow(frame.data, origin="lower", norm=norm, cmap=cmap)
    ax.set_xlabel("RA")
    ax.set_ylabel("DEC")
    ax.set_aspect("equal", adjustable="box")
    return wcs


def fetch_frame_from_fov(
    fov: ZtfCcdQuad,
    products="sci",
    im_type="sciimg.fits",
    force_download=False,
):
    """
    Given a ztf FOV, return the FITs file associated with it.

    This downloads the fits file into the cache.

    Parameters
    ----------
    fov :
        A single CCD Quad FOV.
    products :
        Which data product to fetch.
    im_type :
        Image extension, this must match the products variable.
    force_download :
        Optionally force a re-download if the file already exists in the cache.
    """
    file = fetch_ZTF_file(
        fov.field,
        fov.filefracday,
        fov.filtercode,
        fov.ccdid,
        fov.imgtypecode,
        fov.qid,
        products,
        im_type,
        force_download,
    )
    return fits.open(file)[0]


def fetch_ZTF_file(
    field,
    filefracday,
    filtercode,
    ccdid,
    imgtypecode,
    qid,
    products="sci",
    im_type="sciimg.fits",
    force_download=False,
):
    """
    Fetch a ZTF file directly from the IPAC server, returning the path to where it was
    saved.

    Parameters
    ----------
    field :
        Field identifier number, integer between 1 and ~2200.
    filefracday :
        String describing the fraction of a day, a record keeping string based on time.
    filtercode :
        Which filter was used for the exposure.
    ccdid :
        The CCD identified. (1-16)
    imgtypecode:
        Image type code.
    qid :
        Which quad of the ccd was used. (1-4)
    products :
        Exposure products, "sci" for science images.
    im_type :
        Image extension, this must match the products variable.
    force_download :
        Optionally force a re-download if the file already exists in the cache.

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
        f"{filtercode}_c{ccdid}_{imgtypecode}_q{qid}_{im_type}"
    )

    url = ztf_base + path + file

    return cached_file_download(
        url, force_download=force_download, subfolder="ztf_frames"
    )
