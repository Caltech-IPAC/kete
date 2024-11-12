"""
Field of View
=============

Visualizing Field of Views (FOV).

"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
import kete
import cartopy.crs as ccrs
import matplotlib.ticker as mticker

# Target point
target_point = kete.Vector.from_ra_dec(270.84097, 66.56904)

# Pick your favorite chip, 0, 1, 2, 3
chip_id = 2

# what fraction of the chip do you want the target
x_frac = 0.5
y_frac = 0.5

# what time is it
jd = kete.Time.from_ymd(2028, 4, 1).jd


# mass of the earth in solar masses
earth_mass = 3.00348961546514e-06

# calculating how far the L1 point is from earth
# see https://en.wikipedia.org/wiki/Lagrange_point
mu = earth_mass / (earth_mass + 1)
roots = np.polynomial.Polynomial(
    [1, (mu - 3), (3 - 2 * mu), -mu, +2 * mu, -mu][::-1]
).roots()
l1_ratio = roots[np.argmin(abs(np.abs(roots) - roots))].real

# Calculate an observer at the L1 position
# Approximation of NEOS location
earth = kete.spice.get_state("Earth", jd)
observer = kete.State("NEOS", jd, earth.pos * (1 - l1_ratio), earth.vel)


def place_on_chip(target_ra, target_dec, chip_id, x_frac, y_frac, jd):
    """
    Given a target position, and the chip id, the x/y fraction along the chip, and the
    time, construct a FOV which obeys sunshield restrictions and places the target at
    the desired location.
    """

    def _cost(p):
        pointing = kete.Vector.from_ra_dec(*p)
        rot = kete.neos.sunshield_rotation(observer.pos, pointing)
        fov = kete.neos.NeosVisit(
            kete.neos.FOV_WIDTH,
            kete.neos.FOV_HEIGHT,
            kete.neos.FOV_CHIP_GAP,
            pointing,
            rot,
            observer,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
        )

        chip = fov[chip_id]
        v = (
            (chip.corners[3] - chip.corners[0]) * x_frac
            + (chip.corners[1] - chip.corners[0]) * y_frac
            + chip.corners[0]
        )
        return (target_ra - v.ra) ** 2 + (target_dec - v.dec) ** 2

    f = minimize(_cost, [target_ra, target_dec])
    target_ra, target_dec = f.x
    target_pointing = kete.Vector.from_ra_dec(target_ra, target_dec)

    rot = kete.neos.sunshield_rotation(observer.pos, target_pointing)
    return kete.neos.NeosVisit(
        kete.neos.FOV_WIDTH,
        kete.neos.FOV_HEIGHT,
        kete.neos.FOV_CHIP_GAP,
        target_pointing,
        rot,
        observer,
        0,
        0,
        0,
        0,
        0,
        0,
        1,
    )


fov = place_on_chip(target_point.ra, target_point.dec, chip_id, x_frac, y_frac, jd)

print("RA :", fov.pointing.ra)
print("DEC:", fov.pointing.dec)
print("Rot:", fov.rotation)

chip_id = fov[chip_id]

if fov.pointing.as_ecliptic.lat > 0:
    pole = ccrs.Orthographic(0, 90)
    grid_lines = mticker.FixedLocator([80, 85])
    extent = [0, 359.99, 80, 90]
else:
    pole = ccrs.Orthographic(180, -90)
    grid_lines = mticker.FixedLocator([-80, -85])
    extent = [0, 359.99, -80, -90]
pole.threshold /= 100
plt.figure().add_subplot(projection=pole)
gl = plt.gca().gridlines(draw_labels=True, y_inline=True, dms=True)
gl.ylocator = grid_lines

corners = []
for chip_id in fov:
    c = [c.as_ecliptic for c in chip_id.corners]
    c = [[x.lon, x.lat] for x in c]
    c.append(c[0])
    corners.append(c)

for corner in corners:
    plt.plot(*np.transpose(corner), lw=1, transform=ccrs.Geodetic())
plt.scatter(
    target_point.as_ecliptic.lon,
    target_point.as_ecliptic.lat,
    marker='x',
    c="red", transform=ccrs.Geodetic(), s=15
)

plt.gca().set_extent(extent, crs=ccrs.PlateCarree())
plt.show()
