from __future__ import annotations
import io
from functools import lru_cache
import time
import logging
import warnings
from xml.etree import ElementTree
import matplotlib.pyplot as plt
import requests
import numpy as np
import pandas as pd


from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.visualization import (
    PowerDistStretch,
    AsymmetricPercentileInterval,
    ImageNormalize,
    LinearStretch,
    ZScaleInterval,
)

from .vector import Vector

IRSA_URL = "https://irsa.ipac.caltech.edu"


logger = logging.getLogger(__name__)


@lru_cache()
def IRSA_column_data(table_name, base_url=IRSA_URL, auth=None):
    """
    Retrieve the column data for a specified IRSA table.

    This will return a dataframe containing the column properties of the target table.

    Parameters
    ----------
    table_name :
        The name of the table to query the columns of.
    base_url :
        The URL of the TAPS service to query, this defaults to IRSA.
    auth :
        An optional (username, password), this may be used to access restricted data.
    """
    return query_irsa_tap(
        f"""SELECT * FROM TAP_SCHEMA.columns WHERE table_name='{table_name}'""",
        base_url=base_url,
        auth=auth,
    )


def query_irsa_tap(
    query, upload_table=None, base_url=IRSA_URL, auth=None, timeout=None, verbose=False
):
    """
    Query IRSA's TAP service, optionally upload a table which will be included in the
    query data. The pandas dataframe table will be labeled as `my_table` and columns in
    the query can be used like so:

    .. testcode::
        :skipif: True

        import kete
        import pandas as pd

        # Column names cannot match the column names in the IRSA table you are querying
        # 0 has been added to the end of these column names to satisfy this constraint.

        data = pd.DataFrame([['foo', 56823.933738, 186.249070833, 22.8977],
                            ['bar', 55232.963786, 49.14175, 21.63811111]],
                            columns=['name', 'mjd0', 'ra0', 'dec0'])

        jd = kete.Time.from_mjd(56823.933738).jd

        # This time corresponds to this phase:
        phase = kete.wise.mission_phase_from_jd(jd)

        # Single source table on IRSA is then: phase.source_table

        # The columns of data available in this source table are
        column_information = kete.irsa.IRSA_column_data(phase.source_table)

        # Note that lots of information is available in column_information

        # Now select which columns we want IRSA to return.
        # Using TAP_UPLOAD.my_table.name we can get back the column of data we sent
        columns_to_fetch = "TAP_UPLOAD.my_table.name, mjd, ra, dec"

        query = (f"select {columns_to_fetch} from {phase.source_table} where " +
                "CONTAINS(POINT('J2000',ra,dec)," +
                "         CIRCLE('J2000'," +
                "                TAP_UPLOAD.my_table.ra0," +
                "                TAP_UPLOAD.my_table.dec0," +
                "                0.01)" +
                "         )=1 " +
                " and (((mjd - mjd0) < 0.0001) " +
                " and ((mjd0 - mjd) < 0.0001))")

        result = kete.irsa.query_irsa_tap(query, upload_table=data)

    This is a blocking operation using TAP Async queries. This submits the query,
    receives a response from IRSA containing a URL, then queries that URL for job
    status. This continues until the job either completes or errors.

    Parameters
    ----------
    query :
        An SQL text query.
    upload_table :
        An optional pandas dataframe.
    base_url :
        The URL of the TAPS service to query, this defaults to IRSA.
    auth :
        An optional (username, password), this may be used to access restricted data.
    timeout :
        Timeout for web queries. This will raise an exception if the servers to not
        respond within this time.
    verbose :
        Print status responses as they are fetched from IRSA.
    """
    data = dict(FORMAT="CSV", QUERY=query)
    files = None
    if upload_table is not None:
        data["UPLOAD"] = "my_table,param:table.tbl"

        csv_output = io.StringIO()
        pd.DataFrame(upload_table).to_csv(csv_output, index=False)
        csv_output.seek(0)
        files = {"table.tbl": csv_output.read().encode()}

    submit = requests.post(
        base_url + "/TAP/async", data=data, files=files, auth=auth, timeout=timeout
    )
    submit.raise_for_status()

    tree = ElementTree.fromstring(submit.content.decode())
    element = tree.find("{*}results")

    urls = [v for k, v in element[0].attrib.items() if "href" in k]
    if len(urls) != 1:
        raise ValueError("Unexpected results: ", submit.content.decode())
    url = urls[0]

    phase_url = url.replace("results/result", "phase")

    status = requests.get(phase_url, timeout=timeout, auth=auth)
    status.raise_for_status()

    # Status results can have one of 4 outcomes:
    # QUEUED, EXECUTING, ERROR, COMPLETED

    start = time.time()
    time.sleep(0.15)
    delay = 0.85
    while status.content.decode().upper() in ["QUEUED", "EXECUTING"]:
        elapsed = time.time() - start
        # job is not complete
        if verbose:
            logger.info(
                f"IRSA response ({elapsed:0.1f} sec elapsed): %s",
                status.content.decode(),
            )
        time.sleep(delay)
        status = requests.get(phase_url, timeout=timeout, auth=auth)
        status.raise_for_status()

        # Increase time between queries until there is 30 seconds between.
        # Then continue forever.
        if delay < 30:
            delay += 1

    if status.content.decode().upper() != "COMPLETED":
        raise ValueError("Job Failed: ", status.content.decode())

    result = requests.get(url, timeout=timeout, auth=auth)
    result.raise_for_status()
    return pd.read_csv(io.StringIO(result.text))


def plot_fits_image(fit, percentiles=(0.1, 99.95), power_stretch=0.5, cmap="gray"):
    """
    Plot a FITS image, returning a WCS object which may be used to plot future points
    correctly onto the current image.

    This estimates the standard deviation, subtracts the median, and scales the
    displayed image by number of standard deviations from the median value.

    This returns the WCS which is constructed during the plotting process.

    This will use the existing matplotlib plotting axis if available.

    Parameters
    ----------
    fit:
        Fits file from Astropy.
    percentiles :
        Statistical percentile limit for which data to plot. By default this is set
        to 0.1% and 99.95%. If this is set to `None`, then this uses Astropy's
        `ZScaleInterval`.
    power_stretch :
        The scaling of the intensity of the plot is a power law, this defines the power
        of that power law. By default plots are sqrt scaled. If this is set to 1, then
        this becomes a linear scaling.
    cmap :
        Color map to use for the plot.
    """

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        wcs = WCS(fit.header)

    if not plt.get_fignums():
        plt.figure(dpi=200, figsize=(6, 6), facecolor="w")

    # This is a little fancy footwork to get this to play nice with subplots
    # This gets the current subplot geometry, pulls up the current axis and
    # nukes it. Then in its place it inserts a new axis with the correct
    # projection.
    fig = plt.gcf()
    rows, cols, start, _ = plt.gca().get_subplotspec().get_geometry()
    fig.axes[start].remove()
    ax = fig.add_subplot(rows, cols, start + 1, projection=wcs)
    ax.set_aspect("equal", adjustable="box")
    fig.axes[start] = ax

    if power_stretch == 1:
        stretch = LinearStretch()
    else:
        stretch = PowerDistStretch(power_stretch)

    if percentiles is None:
        interval = ZScaleInterval()
    else:
        interval = AsymmetricPercentileInterval(*percentiles)

    norm = ImageNormalize(fit.data, interval=interval, stretch=stretch)
    data = np.nan_to_num(fit.data, nan=np.nanpercentile(fit.data, 5))

    ax.imshow(data, origin="lower", norm=norm, cmap=cmap)
    ax.set_xlabel("RA")
    ax.set_ylabel("DEC")
    ax.set_aspect("equal", adjustable="box")
    return wcs


def _ra_dec(ra, dec=None):
    """
    Given either a RA/Dec pair, or Vector, return an RA/Dec pair.
    """
    if isinstance(ra, Vector):
        vec = ra.as_equatorial
        ra = vec.ra
        dec = vec.dec
    return ra, dec


def zoom_plot(wcs, ra, dec=None, zoom=100):
    """
    Given a WCS, zoom the current plot to the specified RA/Dec.

    Parameters
    ----------
    wcs :
        An Astropy World Coordinate system from the image.
    ra :
        The RA in degrees, can be a `Vector`, if so then dec is ignored.
    dec :
        The DEC in degrees.
    zoom :
        Optional zoom region in pixels
    """
    ra, dec = _ra_dec(ra, dec)
    pix = wcs.world_to_pixel(SkyCoord(ra, dec, unit="deg"))
    plt.gca().set_xlim(pix[0] - zoom, pix[0] + zoom)
    plt.gca().set_ylim(pix[1] - zoom, pix[1] + zoom)


def annotate_plot(
    wcs,
    ra,
    dec=None,
    text=None,
    px_gap=70,
    length=50,
    lw=1,
    c="red",
    text_color="White",
    style="+",
    text_dx=0,
    text_dy=0,
    text_fs=None,
):
    """
    Add an annotation for a point in a FITS plot, this requires a world coordinate
    system (wcs) as returned by the plotting function above.

    Parameters
    ----------
    wcs :
        An Astropy World Coordinate system from the image.
    ra :
        The RA in degrees, can be a `Vector`, if so then dec is ignored.
    dec :
        The DEC in degrees.
    text :
        Optional text to display.
    px_gap :
        How many pixels should the annotation be offset from the specified RA/DEC.
    length :
        Length of the bars in pixels.
    lw :
        Line width of the marker.
    c :
        Color of the marker, uses matplotlib colors.
    text_color :
        If text is provided, this defines the text color.
    style :
        Style of marker, this can be either "o", "+", or "L".
    text_dx :
        Offset of the text x location in pixels.
    text_dy :
        Offset of the text y location in pixels.
    text_fs :
        Text font size.
    """
    ra, dec = _ra_dec(ra, dec)
    x, y = wcs.world_to_pixel(SkyCoord(ra, dec, unit="deg"))
    total = length + px_gap
    if style == "+":
        plt.plot([x - total, x - px_gap], [y, y], c=c, lw=lw)
        plt.plot([x + px_gap, x + total], [y, y], c=c, lw=lw)
        plt.plot([x, x], [y - px_gap, y - total], c=c, lw=lw)
        plt.plot([x, x], [y + px_gap, y + total], c=c, lw=lw)
    elif style == "L":
        plt.plot([x + px_gap, x + total], [y, y], c=c, lw=lw)
        plt.plot([x, x], [y + px_gap, y + total], c=c, lw=lw)
    elif style == "o":
        plt.scatter(x, y, fc="None", ec=c, s=5 * px_gap, lw=lw)
    else:
        raise ValueError("Style is not recognized, must be one of: o, +, L")

    if text:
        plt.text(x + text_dx, y + text_dy, text, c=text_color, fontsize=text_fs)
