import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
import numpy as np
from .vector import Vector, Frames

__all__ = ["neos_field_of_regard", "add_rectangle", "celestial_sphere"]


def celestial_sphere(
    fig=None,
    extent=None,
    xticks=range(-180, 181, 45),
    yticks=range(-60, 61, 20),
    projection=ccrs.Mollweide,
    x_inline=False,
    y_inline=False,
    deg_format=".0f",
):
    """
    Plot the celestial sphere, returning the axis so that other things may be added.

    Parameters
    ----------
    fig :
        Optional matplotlib figure, if None, then this checks to see if an existing
        figure is open and uses that one. If None are open at all, this creates a new
        figure.
    extent :
        The extent of the sphere to plot, this is `[x_min, x_max, y_min, y_max]`.
        Defaults to `[-150, 150, -50, 50]`.
    xticks :
        Location of the x ticks.
    yticks :
        Location of the y ticks.
    projection :
        Which spherical projection to use, see `cartopy.crs` for projections.
    x_inline :
        Should the x line labels be inline on the plot.
    y_inline :
        Should the y line labels be inline on the plot.
    deg_format :
        Format string for the degree labels of the ticks.
    """
    if extent is None:
        extent = [-150, 150, -50, 50]
    elif extent == "None":
        extent = None

    if fig is None and len(plt.get_fignums()) > 0:
        fig = plt.gcf()
    elif fig is None:
        fig = plt.figure(dpi=200)

    # set finer resolution when plotting lines. The defaults in cartopy often
    # result in jagged shaped polygons
    proj = projection()
    # pylint: disable=protected-access
    proj._threshold = proj._threshold / 10.0

    ax = plt.axes(projection=proj)
    if extent is not None:
        ax.set_extent(extent, ccrs.PlateCarree())
    else:
        # ugly hack to force the full extent to plot
        ax.scatter(
            [-179.5, 179.5, 0, 0],
            [0, 0, -50, 50],
            c="None",
            transform=ccrs.PlateCarree(),
        )
    gl = ax.gridlines(
        draw_labels=True,
        xlocs=xticks,
        ylocs=yticks,
        x_inline=x_inline,
        y_inline=y_inline,
        transform=ccrs.PlateCarree(),
    )

    gl.xformatter = LongitudeFormatter(
        direction_label=False,
        number_format=deg_format,
        dateline_direction_label=False,
    )
    gl.yformatter = LatitudeFormatter(
        direction_label=False,
        number_format=deg_format,
    )
    gl.xlabel_style = {"size": 6}
    gl.ylabel_style = {"size": 6}
    ax.invert_xaxis()
    return ax


def add_rectangle(xy, width, height, ax=None, **kwargs):
    """
    Add a rectangular patch to an existing celestial sphere plot.

    Parameters
    ----------
    xy :
        Center position of the rectangle.
    width :
        Width of the rectangle.
    height :
        Height of the rectangle.
    ax :
        Optional matplotlib axes, if not found, the current active one is used.
    **kwargs:
        All other kwargs are passed directly to the matplotlib.patched.Rectangle class.
    """
    if ax is None:
        ax = plt.gca()

    xy[0] -= width / 2
    xy[1] -= height / 2

    ax.add_patch(
        mpatches.Rectangle(
            xy=xy,
            width=width,
            height=height,
            transform=ccrs.PlateCarree(),
            **kwargs,
        )
    )


def neos_field_of_regard(
    state,
    include_equatorial=False,
    include_galactic=False,
    elon_range=None,
    lat_range=None,
):
    """
    Plot NEO Surveyors field of regard.

    Parameters
    ----------
    state :
        State defining the position of the observer with respect to the Sun.
    include_equatorial :
        Plot an indicator line for the equatorial plane.
    include_galactic :
        Plot an indicator line for the galactic plane.
    elon_range :
        Range of the solar longitudinal distance where Surveyor will operate.
        Defaults to [45, 120]
    lat_range :
        Range of the Latitude where Surveyor will operate.
        Defaults to [-40, 40].
    """
    elon_range = [45, 120] if elon_range is None else elon_range
    lat_range = [-40, 40] if lat_range is None else lat_range
    sun_pos = -state.pos.as_ecliptic
    ax = celestial_sphere(extent="None")
    ax.scatter(sun_pos.lon, sun_pos.lat, c="Orange", transform=ccrs.PlateCarree())
    width = elon_range[1] - elon_range[0]
    height = lat_range[1] - lat_range[0]
    y_center = (lat_range[1] + lat_range[0]) / 2
    x_center = (elon_range[1] + elon_range[0]) / 2

    add_rectangle([sun_pos.lon - x_center, y_center], width, height, alpha=0.1)
    add_rectangle([sun_pos.lon + x_center, y_center], width, height, alpha=0.1)

    lon_steps = np.linspace(0, 360, 100)
    frames = []
    if include_equatorial:
        frames.append(Frames.Equatorial)
    if include_galactic:
        frames.append(Frames.Galactic)

    for frame in frames:
        vecs = [Vector.from_el_az(0, step, 1, frame=frame) for step in lon_steps]

        vecs = [v.as_ecliptic for v in vecs]
        pos = np.array([[v.az, v.el] for v in vecs])

        # roll the indices so that it plots pretty
        # this is just to make the visualization nice
        pos[:, 0] = (pos[:, 0] + 180) % 360 - 180
        idx = np.argmax(abs(np.diff(pos[:, 0])))
        pos = np.roll(pos, -idx - 1, axis=0)

        ax.plot(*pos.T, label=frame, transform=ccrs.PlateCarree(), lw=0.5)
        vec_0 = Vector([1, 0, 0], frame=frame).as_ecliptic
        ax.scatter(vec_0.az, vec_0.el, transform=ccrs.PlateCarree(), s=10, marker="x")

    return ax
