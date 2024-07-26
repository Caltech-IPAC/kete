"""
Annotate a Comet plot
=====================

Given a specific FITs file containing a comet, annotate the plot with orbital
information, such as diection of motion.
"""

from astropy.wcs import WCS
import neospy
import astropy
import numpy as np
import neospy
import matplotlib.pyplot as plt
import numpy as np


# This is comet NEOWISE as observed by ZTF
state = neospy.HorizonsProperties.fetch("C/2020 F3").state

# Specific frame information for the ZTF frame where neowise was imaged.
# See the tutorials for KONA and Precovery for more information on how to
# find all the times an object is observed.
frame_info = dict(
    field=752,
    filefracday=20200723164595,
    ccdid=2,
    filtercode="zg",
    imgtypecode="o",
    qid=3,
)

# Load the fits file for this ZTF frame.
frame = astropy.io.fits.open(neospy.ztf.fetch_ZTF_file(**frame_info))[0]

# Grab frame information from this file
jd = neospy.Time(frame.header["OBSJD"], scaling="utc").jd
frame_wcs = WCS(frame.header)

corners = []
dx, dy = frame_wcs.pixel_shape
for x, y in zip([dx, dx, 0, 0], [0, dy, dy, 0]):
    coord = frame_wcs.pixel_to_world(x, y).icrs
    corners.append(neospy.Vector.from_ra_dec(coord.ra.deg, coord.dec.deg))
observer_loc = neospy.spice.mpc_code_to_ecliptic("ZTF", jd)

# Build a neospy FOV for this frame
fov = neospy.ZtfCcdQuad(corners, observer_loc, maglimit=np.nan, fid=1, **frame_info)

# Compute the observation information for the comet in this frame
vis = neospy.fov_state_check([state], [fov])[0]


def plot_vector(wcs, vec_a, vec_b, label, x=0.2, y=0.2, c="w", length=0.1, **kwargs):
    """
    Given a world coordinate system, and 2 positions, plot a projected vector from
    the vec_a toward vec_b in the given WCS. This is a fully generic plotting tool
    which is used immediately below to plot the relavent vectors.
    """
    vec_a = vec_a.as_equatorial
    vec_b = vec_b.as_equatorial
    pixels = wcs.world_to_pixel_values([vec_a.ra, vec_b.ra], [vec_a.dec, vec_b.dec])
    diff_dir = np.diff(pixels, axis=1).ravel()
    diff_dir /= np.linalg.norm(diff_dir)

    shape = wcs.array_shape
    length *= max(shape)
    x *= shape[0]
    y *= shape[1]

    kwargs["lw"] = kwargs.get("lw", 0.01)
    kwargs["width"] = kwargs.get("width", 20)
    plt.arrow(x, y, *diff_dir * length, color=c, **kwargs)
    plt.text(*(diff_dir * 1.5 * length + [x, y]), label, c=c, ha="center", va="center")
    plt.scatter(x, y, c="grey", s=5)


def plot_vectors(wcs, state, fov, x=0.2, y=0.2):
    """
    Plot 4 vectors using the WCS, the state of the object, and the FOV.

    The vectors plotted are:
    - Negative Velocity vector of the object, showing its motion.
    - Position vector, showing the direction away from the sun.
    - North vector, showing the Equatorial North.
    - East vector, showing Equatorial East.
    """

    past_vec = (
        neospy.propagate_n_body([state], state.jd - 0.05)[0].pos - fov.observer.pos
    )
    sun_vec = (state.pos * 1.001) - fov.observer.pos
    vec = (state.pos - fov.observer.pos).as_equatorial
    north_vec = neospy.Vector.from_ra_dec(vec.ra, vec.dec + 0.01)
    east_vec = neospy.Vector.from_ra_dec(vec.ra + 0.01, vec.dec)

    plot_vector(wcs, vec, past_vec, r"-$v$", x=x, y=y, c="r")
    plot_vector(wcs, vec, sun_vec, r"r$_\odot$", x=x, y=y, c=(0, 0.5, 1))
    plot_vector(wcs, vec, north_vec, r"N", c="grey", x=x, y=y, ls="--", lw=0.1)
    plot_vector(wcs, vec, east_vec, r"E", c="grey", x=x, y=y, ls="--", lw=0.1)


# Plot the final results
plt.figure(dpi=200)
wcs = neospy.ztf.plot_ztf_fov(vis.fov)
vec = vis.obs_vecs()[0].as_equatorial
neospy.irsa.annotate_plot(wcs, vec.ra, vec.dec, px_gap=0, length=25)
plot_vectors(frame_wcs, vis[0], vis.fov, y=0.85)
# for some reason astropy for this frame is plotting the y-axis inverted
# this just un-inverts it.
plt.gca().invert_yaxis()
plt.show()
