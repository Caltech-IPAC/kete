from .data import cached_file_download

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

    return cached_file_download(url, force_download=force_download, subfolder="ztf")
