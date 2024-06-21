"""
Plot visible phases and mags from the MPC
=========================================
"""

import neospy
import numpy as np
import matplotlib.pyplot as plt


# Fetch known orbits from the MPC
orbits = neospy.mpc.fetch_known_orbit_data()

# Subselect NEO population
neos = neospy.population.neo(orbits.peri_dist, orbits.ecc)
neo_subset = orbits[neos]

# Convert the dataframe to States
neos = neospy.mpc.table_to_states(neo_subset)

# Note that some of the objects will have NAN positions, which is result of
# them hitting the earth after discovery.

# Pick your favorite time to observe
jd = neospy.Time("2023-11-04 12:00:00.000", "iso", "utc").jd

# Move the entire population of asteroids to that time using 2-body
# mechanics, this can be directly substituted with propagate_n_body if you
# want more precision.
sun2earth = neospy.SpiceKernels.state("Earth", jd).pos
states = neospy.propagate_two_body(neos, jd)

# Compute the expected V-mags for these objects at this time
mags = []
phases = []
elongs = []
for i, state in enumerate(states):
    sun2obj = state.pos
    obj2earth = sun2earth - sun2obj
    mags.append(
        neospy.flux.hg_apparent_mag(
            -sun2obj,
            obj2earth,
            neo_subset.g_phase.iloc[i],
            neo_subset.h_mag.iloc[i],
        )
    )

    # Do some vector math to get the phase angles.
    phases.append(obj2earth.angle_between(-sun2obj))
    elongs.append((-sun2earth).angle_between(-obj2earth))

mags = np.array(mags)
elongs = np.array(elongs)
phases = np.array(phases)


# Plot the results
plt.subplot(2, 1, 1)
plt.hist(phases[elongs > 90], histtype="step", bins=20)
plt.xlabel("Phase (Deg)")
plt.ylabel("Counts")

plt.subplot(2, 1, 2)
plt.scatter(mags[elongs > 90], phases[elongs > 90], s=1)
plt.xlim(0, 25)
plt.xlabel("Visible Mag")
plt.ylabel("Phase (Deg)")
plt.tight_layout()
plt.show()
