NEO Surveyor Simulation Tools
=============================

NEOSpy is a collection of tools which enable the simulation of ground or spaced based
surveys of minor planets. Included here are:

   - Orbit propagation code capable of accurately calculating the orbits for many
     thousands of minor planets over decades on a laptop.
   - Thermal modeling, including NEATM (Near Earth Asteroid Thermal Model) and the FRM
     (Fast Rotator Model) for asteroids. Including support for non-spherical asteroids.
   - Optical modeling.
   - Interfaces to IPAC's IRSA, JPL's Horizons, and the Minor Planet Centers web tools.
   - Multi-core support for reading SPICE SPK files (SPICE itself is not used).
   - Geometric calculations to check if an object is visible from an observer's
     position.

Many of these tools are independent from one another, and can be combined in many ways
to answer questions. There are a collection of worked examples below which demonstrate
some of the capabilities.

It is recommended that the background reading section be covered before using NEOSpy,
as it is very important to understand the coordinate frames and reference times used.
Conversions such as UTC <-> TDB time can result in over a minute of difference in the
position of minor planets, which for fast moving objects is long enough that they may
move out of frame.

.. toctree::
   :maxdepth: 2

   background
   code_structure
   propagation

.. toctree::
   :maxdepth: 1

   installation
   auto_examples/index
   tutorials/index
   api/api

