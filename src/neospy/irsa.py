from __future__ import annotations
import warnings
import io
from functools import lru_cache
import requests
import numpy as np
import pandas as pd  # type: ignore
import matplotlib.pyplot as plt  # type: ignore


from astropy.wcs import WCS  # type: ignore
from astropy.coordinates import SkyCoord  # type: ignore

IRSA_URL = "https://irsa.ipac.caltech.edu"


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
    query, upload_table=None, base_url=IRSA_URL, auth=None, timeout=None
):
    """
    Query IRSA's TAP service, optionally upload a table which will be included in the
    query data. The pandas dataframe table will be labeled as `my_table` and columns in
    the query can be used like so:

    .. testcode::
        :skipif: True

        import neospy
        import pandas as pd

        # Column names cannot match the column names in the IRSA table you are querying
        # 0 has been added to the end of these column names to satisfy this constraint.

        data = pd.DataFrame([['foo', 56823.933738, 186.249070833, 22.8977],
                            ['bar', 55232.963786, 49.14175, 21.63811111]],
                            columns=['name', 'mjd0', 'ra0', 'dec0'])

        jd = neospy.Time(56823.933738, 'mjd').jd

        # This time corresponds to this phase:
        phase = neospy.wise.mission_phase_from_jd(jd)

        # Single source table on IRSA is then: phase.source_table

        # The columns of data available in this source table are
        column_information = neospy.irsa.IRSA_column_data(phase.source_table)

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

        result = neospy.irsa.query_irsa_tap(query, upload_table=data)

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
    """
    data = dict(FORMAT="CSV", QUERY=query)
    files = None
    if upload_table is not None:
        data["UPLOAD"] = "my_table,param:table.tbl"

        csv_output = io.StringIO()
        pd.DataFrame(upload_table).to_csv(csv_output, index=False)
        csv_output.seek(0)
        files = {"table.tbl": csv_output.read().encode()}

    res = requests.post(
        base_url + "/TAP/sync", data=data, files=files, auth=auth, timeout=timeout
    )
    if res.text[0] == "<":
        raise ValueError("Query returned non-table results: \n\n" + res.text)
    return pd.read_csv(io.StringIO(res.text))


def plot_fits_image(fit, vmin=-3, vmax=7, cmap="bone_r"):
    """
    Plot a FITS image, returning a WCS object which may be used to plot future points
    correctly onto the current image.
    """
    data = np.nan_to_num(fit.data)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        wcs = WCS(fit.header)
        if not plt.get_fignums():
            plt.figure(dpi=120, figsize=(6, 6), facecolor="w")
        ax = plt.subplot(projection=wcs)

    data_no_bkg = data - np.median(data)
    # np.std below is doing a full frame std, which grabs the flux
    # from stars and so is not a great estimate for W1 and W2
    data_perc = np.percentile(data_no_bkg, [16, 84])
    good_std = (data_perc[1] - data_perc[0]) / 2.0

    ax.pcolormesh(data_no_bkg / good_std, cmap=cmap, vmin=vmin, vmax=vmax)
    plt.gca().set_aspect("equal", adjustable="box")
    return wcs


def annotate_plot(
    wcs, ra, dec, text=None, px_gap=70, length=50, lw=1, c="red", text_color="White"
):
    """
    Add an annotation for a point in a FITS plot, this requires a world coordinate
    system (wcs) as returned by the plotting function above.
    """
    x, y = wcs.world_to_pixel(SkyCoord(ra, dec, unit="deg"))
    total = length + px_gap
    plt.plot([x - total, x - px_gap], [y, y], c=c, lw=lw)
    plt.plot([x + px_gap, x + total], [y, y], c=c, lw=lw)
    plt.plot([x, x], [y - px_gap, y - total], c=c, lw=lw)
    plt.plot([x, x], [y + px_gap, y + total], c=c, lw=lw)
    if text:
        plt.text(x, y + px_gap / 2, "     " + text, c=text_color)
