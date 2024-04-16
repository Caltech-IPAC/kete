from __future__ import annotations


import os
import shutil
from collections import namedtuple
from functools import lru_cache
import tempfile
import requests
import warnings
import matplotlib.pyplot as plt  # type: ignore
import numpy as np
from typing import Optional

from astropy.io import fits  # type: ignore
from astropy.wcs import WCS  # type: ignore
from astropy import units as u  # type: ignore
from astropy.coordinates import SkyCoord  # type: ignore


from .data import cache_path
from .time import Time
from .spice import SpiceKernels
from .vector import Vector
from .irsa import IRSA_URL, query_irsa_tap

# pylint: disable-next=import-error
from ._rust import WiseCmos, FOVList  # type: ignore


# This table contains calibrated flux correction values for the 4 different WISE bands.
COLOR_CORR = np.array(
    [  # Tbb  K_W1     K_W2    K_W3     K_W4
        [100, 17.2062, 3.9096, 2.6588, 1.0032],
        [110, 9.3213, 3.1120, 2.1424, 0.9955],
        [120, 6.5668, 2.5973, 1.8071, 0.9905],
        [130, 5.0963, 2.2458, 1.5776, 0.9873],
        [140, 4.1630, 1.9949, 1.4139, 0.9853],
        [150, 3.5203, 1.8096, 1.2935, 0.9841],
        [160, 3.0554, 1.6687, 1.2027, 0.9834],
        [170, 2.7070, 1.5591, 1.1328, 0.9830],
        [180, 2.4387, 1.4722, 1.0782, 0.9830],
        [190, 2.2273, 1.4021, 1.0350, 0.9831],
        [200, 2.0577, 1.3448, 1.0006, 0.9833],
        [210, 1.9194, 1.2975, 0.9728, 0.9836],
        [220, 1.8051, 1.2579, 0.9503, 0.9839],
        [230, 1.7095, 1.2245, 0.9321, 0.9843],
        [240, 1.6286, 1.1961, 0.9172, 0.9847],
        [250, 1.5596, 1.1717, 0.9051, 0.9851],
        [260, 1.5003, 1.1508, 0.8952, 0.9855],
        [270, 1.4489, 1.1326, 0.8872, 0.9860],
        [280, 1.4040, 1.1167, 0.8808, 0.9864],
        [290, 1.3647, 1.1029, 0.8755, 0.9868],
        [300, 1.3300, 1.0907, 0.8714, 0.9871],
        [310, 1.2993, 1.0799, 0.8682, 0.9875],
        [320, 1.2719, 1.0704, 0.8657, 0.9879],
        [330, 1.2474, 1.0619, 0.8639, 0.9882],
        [340, 1.2255, 1.0543, 0.8626, 0.9886],
        [350, 1.2057, 1.0476, 0.8617, 0.9889],
        [360, 1.1879, 1.0416, 0.8612, 0.9892],
        [370, 1.1718, 1.0361, 0.8611, 0.9895],
        [380, 1.1571, 1.0313, 0.8613, 0.9898],
        [390, 1.1438, 1.0269, 0.8616, 0.9901],
        [400, 1.1316, 1.0229, 0.8622, 0.9903],
    ]
)
# All lists below are indexed as follows W1 = [0], W2 = [1]. W3 = [2], W4 = [3]

# Flux in the reflection model should be scaled by these values.
SUN_COLOR_CORRECTION: list[float] = [1.0049, 1.0193, 1.0024, 1.0012]

# Non-color corrected values for zero mag corrections in Janskys.
# magnitude can then be computed via -2.5 log10(flux Jy / zero_point)
ZERO_MAGS: list[float] = [306.681, 170.663, 29.0448, 8.2839]

# Color Corrected Zero Mags in units of Jy
ZERO_MAGS_COLOR_CORRECTED: list[float] = [
    ZERO_MAGS[0],
    ZERO_MAGS[1],
    1.08 * ZERO_MAGS[2],
    0.96 * ZERO_MAGS[3],
]

# Non-color corrected values for the effective central wavelength of the bands (nm)
BAND_WAVELENGTHS: list[float] = [3352.6, 4602.8, 11560.8, 22088.3]

# Color Corrected Band Wavelengths in units of nm
BAND_WAVELENGTHS_COLOR_CORRECTED: list[float] = [
    BAND_WAVELENGTHS[0],
    BAND_WAVELENGTHS[1],
    0.96 * BAND_WAVELENGTHS[2],
    1.025 * BAND_WAVELENGTHS[3],
]

# The data in the color correction table above is relatively well fit to a 1 / f(x)
# where f is a 4th order polynomial. It is the least accurate for the final band, but
# W4 is almost a constant. These fits perform much better than linear except near 100k.
# These functions should not be referenced directly, as they also exist in the rust
# code and are only left here for reference to how the rust was constructed.

# from numpy.polynomial import Polynomial
# _COLOR_FITS = [
#     Polynomial.fit(WISE_COLOR_CORR[:, 0], 1 / WISE_COLOR_CORR[:, 1], 4),
#     Polynomial.fit(WISE_COLOR_CORR[:, 0], 1 / WISE_COLOR_CORR[:, 2], 4),
#     Polynomial.fit(WISE_COLOR_CORR[:, 0], 1 / WISE_COLOR_CORR[:, 3], 4),
#     Polynomial.fit(WISE_COLOR_CORR[:, 0], 1 / WISE_COLOR_CORR[:, 4], 4),
# ]

# WISE Field of view is ~47 arc-minutes square.
FOV_WIDTH: float = 47 / 60

# Convert directly from DN to Jy using this Jy/DN conversion factor
# from https://wise2.ipac.caltech.edu/docs/release/prelim/expsup/sec2_3f.html
DN_TO_JY = [1.9350e-06, 2.7048e-06, 2.9045e-06, 5.2269e-05]


MissionPhase = namedtuple(
    "MissionPhase",
    [
        "name",
        "jd_start",
        "jd_end",
        "bands",
        "frame_url",
        "frame_meta_table",
        "source_table",
    ],
)

MISSION_PHASES = {
    "Cryo": MissionPhase(
        name="Cryo",
        jd_start=Time.from_ymd(2009, 12, 14).jd,
        jd_end=2455414.941783008,
        bands=(1, 2, 3, 4),
        frame_url=IRSA_URL + "/ibe/data/wise/allsky/4band_p1bm_frm/",
        frame_meta_table="allsky_4band_p1bs_frm",
        source_table="allsky_4band_p1bs_psd",
    ),
    "3-Band": MissionPhase(
        name="3-Band",
        jd_start=2455414.9417830084,
        jd_end=2455469.278276,
        bands=(1, 2, 3),
        frame_url=IRSA_URL + "/ibe/data/wise/cryo_3band/3band_p1bm_frm/",
        frame_meta_table="allsky_3band_p1bs_frm",
        source_table="allsky_3band_p1bs_psd",
    ),
    "Post-Cryo": MissionPhase(
        name="Post-Cryo",
        jd_start=2455469.278277,
        jd_end=2455593.96119803,
        bands=(1, 2),
        frame_url=IRSA_URL + "/ibe/data/wise/postcryo/2band_p1bm_frm/",
        frame_meta_table="allsky_2band_p1bs_frm",
        source_table="allsky_2band_p1bs_psd",
    ),
    "Reactivation_2014": MissionPhase(
        name="Reactivation_2014",
        jd_start=Time.from_ymd(2013, 12, 13).jd,
        jd_end=Time.from_ymd(2015, 1, 1).utc.jd,
        bands=(1, 2),
        frame_url=IRSA_URL + "/ibe/data/wise/neowiser/p1bm_frm/",
        frame_meta_table="neowiser_p1bs_frm",
        source_table="neowiser_p1bs_psd",
    ),
    "Reactivation_2023": MissionPhase(
        name="Reactivation_2023",
        jd_start=Time.from_ymd(2023, 1, 1).utc.jd,
        jd_end=Time.from_ymd(2023, 12, 13.121151733).jd,
        bands=(1, 2),
        frame_url=IRSA_URL + "/ibe/data/wise/neowiser/p1bm_frm/",
        frame_meta_table="neowiser_p1bs_frm",
        source_table="neowiser_p1bs_psd",
    ),
}

for year in range(2015, 2023):
    MISSION_PHASES[f"Reactivation_{year}"] = MissionPhase(
        name=f"Reactivation_{year}",
        jd_start=Time.from_ymd(year, 1, 1).utc.jd,
        jd_end=Time.from_ymd(year + 1, 1, 1).utc.jd,
        bands=(1, 2),
        frame_url=IRSA_URL + "/ibe/data/wise/neowiser/p1bm_frm/",
        frame_meta_table="neowiser_p1bs_frm",
        source_table="neowiser_p1bs_psd",
    )
MISSION_PHASES = dict(sorted(MISSION_PHASES.items(), key=lambda x: x[1].jd_start))


def WISE_phase_from_jd(jd: float):
    """
    Return which mission phase is associated with the provided JD. Returns None if no
    matching mission phase.
    """
    for mission in MISSION_PHASES.values():
        if mission.jd_start <= jd and jd < mission.jd_end:
            return mission
    return None


def WISE_phase_from_scan(scan_id: str) -> Optional[MissionPhase]:
    """
    The WISE mission has different phases based on the cryogen on board and the
    hibernation period, resulting in different useable bands and different locations in
    the IRSA archive.

    Parameters
    ----------
    scan_id :
        The scan id to be converted, this is fully sufficient to decide which mission
        phase.
    """
    scan_num = int(scan_id[0:5])
    if scan_id[-1] in ["r", "s", "t"]:
        # late reactivation
        if scan_num <= 45686:  # final scan number of year 9
            return MISSION_PHASES["Reactivation_2022"]
        return None
    else:
        # first pass through the 5 digit scan number cycle
        if scan_num > 40000:
            # early reactivation
            return MISSION_PHASES["Reactivation_2014"]
        elif scan_num > 8745:
            return MISSION_PHASES["Post-Cryo"]
        elif scan_num <= 7101 and scan_id[-1] in "a":
            return MISSION_PHASES["Cryo"]
        elif scan_num <= 7097 and scan_id[-1] == "b":
            return MISSION_PHASES["Cryo"]
        elif scan_num <= 6217 and scan_id[-1] == "c":
            return MISSION_PHASES["Cryo"]
        elif scan_num <= 831 and scan_id[-1] == "j":
            return MISSION_PHASES["Cryo"]
        else:
            return MISSION_PHASES["3-Band"]


def cache_WISE_frame(scan_id, frame_num, band=3, im_type="int"):
    """
    Download and save a WISE frame to the cache, unless the file already exists.

    Returns the cached filepath of the file.
    """
    scan_id = str(scan_id)
    frame_num = f"{frame_num:03d}"
    scan_grp = str(scan_id[-2:])
    band = f"w{band:1d}"
    ext = "fits" if im_type == "int" else "fits.gz"
    filename = f"{scan_id}{frame_num}-{band}-{im_type}-1b.{ext}"

    dir_path = os.path.join(cache_path(), "wise_frames")
    full_path = os.path.join(dir_path, filename)

    if not os.path.isfile(full_path):
        if not os.path.isdir(dir_path):
            os.makedirs(dir_path)

        phase = WISE_phase_from_scan(scan_id)
        url = f"{phase.frame_url}{scan_grp}/{scan_id}/{frame_num}/{filename}"
        res = requests.get(url, timeout=30)
        if res.status_code == 404:
            raise ValueError(
                f"scan {scan_id}  frame {frame_num}  {band}  {im_type} was not found."
            )
        if res.status_code != 200:
            raise ValueError(res.content.decode())

        with tempfile.NamedTemporaryFile("wb") as f:
            f.write(res.content)
            shutil.copyfile(f.name, full_path)

    return full_path


def fetch_WISE_frame(scan_id, frame_num, band=3, im_type="int", cache=True):
    """
    Fetch a WISE image directly from the IPAC server, returning a FITS image.

    WISE Frames are stored in 3 `im_type` flavors:
    - 'int' - Intensity frames.
    - 'unc' - Uncertainty frames (these are saved as `.fits.gz` files)
    - 'msk' - Mask frames (these are saved as `.fits.gz` files)

    If `cache` is `True`, then the image is stored in the image cache directory.
    """

    ext = "fits" if im_type == "int" else "fits.gz"

    if cache:
        path = cache_WISE_frame(scan_id, frame_num, band=band, im_type="int")
        return fits.open(path)[0]

    scan_id = str(scan_id)
    frame_num = f"{frame_num:03d}"
    scan_grp = str(scan_id[-2:])
    band = f"w{band:1d}"

    filename = f"{scan_id}{frame_num}-{band}-{im_type}-1b.{ext}"

    phase = WISE_phase_from_scan(scan_id)
    url = f"{phase.frame_url}{scan_grp}/{scan_id}/{frame_num}/{filename}"
    res = requests.get(url, timeout=30)

    with tempfile.NamedTemporaryFile("wb") as f:
        f.write(res.content)
        try:
            fit = fits.open(f.name)[0]
        except OSError as exc:
            raise ValueError(
                "Failed to load frame at url ", url, phase, res.content
            ) from exc
    return fit


def plot_4band(
    scan_id,
    frame_num,
    ra=None,
    dec=None,
    zoom=True,
    bands=None,
):
    """
    Plot a 2x2 grid of images showing the W1, W2, W3, and W4 data for the scan/frame
    pair, centered on the ra,dec provided
    """
    ecolor = [None, "blue", "green", "#ff7700", "red"]

    if bands is None:
        bands = WISE_phase_from_scan(scan_id).bands

    plt.figure(dpi=120, figsize=(8, 8), facecolor="w")
    plt.suptitle(f"Scan: {scan_id}  Frame: {frame_num}")
    plt.subplots_adjust(
        wspace=0.3, hspace=0.3, left=0.1, right=0.9, top=0.9, bottom=0.1
    )
    for band in bands:
        try:
            fit = fetch_WISE_frame(
                scan_id,
                frame_num,
                band=band,
            )
        except OSError:
            continue
        data = np.nan_to_num(fit.data)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            wcs = WCS(fit.header)
            ax = plt.subplot(2, 2, band, projection=wcs)

        data_no_bkg = data - np.median(data)
        # np.std below is doing a full frame std, which grabs the flux
        # from stars and so is not a great estimate for W1 and W2
        data_perc = np.percentile(data_no_bkg, [16, 84])
        good_std = (data_perc[1] - data_perc[0]) / 2.0

        ax.pcolormesh(np.clip(data_no_bkg / good_std, -3, 7), cmap="bone_r")
        plt.gca().set_aspect("equal", adjustable="box")

        if ra is not None and dec is not None:
            loc = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame="icrs")
            plt.scatter(
                loc.ra,
                loc.dec,
                transform=plt.gca().get_transform("world"),
                edgecolor=ecolor[band],
                s=100,
                linewidth=1,
                facecolor="none",
            )
            if zoom and len(np.atleast_1d(ra).ravel()) > 1:
                warnings.warn("More than one object found in file, cannot zoom.")
            elif zoom:
                pixloc = wcs.world_to_pixel(loc)
                if band == 4:
                    # band 4's pixels are twice as big as the other three, so zoom more
                    # to ensure the scene is the same
                    pixspan = 37 * float(zoom)
                else:
                    pixspan = 75 * float(zoom)
                plt.xlim(pixloc[0] - pixspan, pixloc[0] + pixspan)
                plt.ylim(pixloc[1] - pixspan, pixloc[1] + pixspan)
        plt.xlabel("RA")
        plt.ylabel("Dec")
        plt.title(f"W{band:1d}")


@lru_cache()
def fetch_WISE_fovs(phase):
    """
    Load all FOVs taken during the specified mission phase of WISE.

    This will download and cache all FOV information for the given mission phase from
    IRSA.

    Parameters
    ----------
    phase :
        A mission phase object.
    """
    cache_dir = cache_path()
    filename = os.path.join(cache_dir, f"wise_{phase.name}_fovs.bin")

    if os.path.isfile(filename):
        return FOVList.load(filename)

    table = phase.frame_meta_table
    cols = ["scan_id", "frame_num", "mjd", "ra", "dec"]
    mjd_start = Time(phase.jd_start).utc.mjd
    mjd_end = Time(phase.jd_end).utc.mjd

    res = query_irsa_tap(
        f"SELECT {', '.join(cols)} FROM {table} "
        f"WHERE mjd >= {mjd_start} and mjd < {mjd_end}"
    )

    # Adding 4.4 seconds for the offset due to W3/4 exposures taking 8.8 seconds
    res["jd"] = np.array(Time(list(res.mjd), "mjd", "utc").jd) + 5.092592592592592e-05

    fovs = []
    for row in res.itertuples():
        state = SpiceKernels.state("WISE", row.jd)

        pointing = Vector.from_ra_dec(row.ra, row.dec).as_ecliptic
        fov = WiseCmos(
            pointing,
            0.0,
            state,
            row.frame_num,
            row.scan_id,
        )

        fovs.append(fov)
    fovs = FOVList(fovs)
    fovs.sort()
    fovs.save(filename)
    return fovs
