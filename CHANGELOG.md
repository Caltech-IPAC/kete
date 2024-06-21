# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [Unreleased]

### Added

- Added `Time.year_float` which converts the `Time` object to the Year as a float.
- Added `SimultaneousState.obs_vecs` which creates vectors from the observer to the
  objects contained within the state.

### Changed

- Significant rewrite of the SPICE kernel file management, this rewrite is required
  so that future work can enable writing SPICE kernel files.
- Moved much of the python documentation into the rust wrappers, and removed the
  remaining empty python files.

### Fixed

- Orbital Elements now correctly computes `true_anomaly`, `mean_anomaly`,
  `eccentric_anomaly`, `semi_major_axis`, and `mean_motion` for parabolic and
  hyperbolic orbits.
- Typo in Observability example was leading to an incorrect Magnitude value.
- Astropy WCS warnings are now suppressed.
- `SpiceKernels.moon_illuminated_frac` now agrees with JPL Horizons, this was doing a
  geometric calculation which was correct, but with a different interpretation than
  desired.

### Removed

- Removed support for SPK Files of type 3, these should be a trivial change from type 2
  however I do not have access to a file of type 3 for testing. Because it cannot be
  validated at the moment, the code has been removed.


## [0.2.3] - 2024 - 6 - 12

### Added

- Added computation for Tisserand's Parameter.
- Added support for NEOS Visit FOVs, which are joint FOVs containing 4 rectangles.
- Added python interface to WISE Color Correction functions.
- Added support for querying static sky sources in FOVs.
- SpiceKernels can now load SPK files directly, without needing the files to be in the
  cache.
- SpiceKernels now supports displaying the contents of the headers of SPk/PCK files.
- SpiceKernels can now return the available loaded SPK segments of an object.

### Changed

- Renamed "HorizonsCovariance" to "Covariance" and generalized it to support both
  cometary and cartesian representations.
- Made api more consistent for conversion to State objects by renaming all instances of
  `.as_state` to `.state`.
- Improved performance of integration and propagation by about 10%.
- Restructured the Rust FOVs to be organized by observatory.
- Renamed `SpiceKernels.cache_kernel_reload` to `kernel_reload` and changed it's args.

### Fixed

- Astropy will no longer warning about deprecated WCS header values for NEOWISE images.


## [0.2.2] - 2024 - 5 - 20

### Added

- Added support for downloading ZTF full fields FOV information into the neospy cache.
- Added simple lookup for observatory ecliptic state using MPC observatory codes.
- Added 'Tutorials' to the documentation, these are larger worked examples which do not
  build at the same time as the standard docs. They are designed to be examples which
  take significant compute, and may run for many minutes.
- Added Tutorials for 'KONA' and 'WISE Precovery'.

### Changed

- Moved population definitions to the base level and renamed it `population`.
- Changed references to diameters in the broken power law sampler.
- FOV propagation tests now return a flatten list of states.
- Renamed `data` to `cache` to more accurately represent its function.

### Fixed

- Coordinate frame conversion was resulting in incorrect coordinate positions when
  passing Equatorial based frames to some FOV related functions. Specifically it was
  failing to convert the Equatorial frame to Ecliptic frame before performing orbit
  calculations.

### Removed

- Removed construction of populations entirely.
- Removed folder for `population`, the remaining contents were moved up to the base
  level of the package.
- Removed PDS related tools.


## [0.2.1]  - 2024 - 5 - 13

### Added

- Added support for querying SPK kernels of a specific epoch of fit from JPL Horizons.
- Added support for looking up the MPC preferred unpacked designation for asteroids.

### Changed

- Change documentation to clarify that designations are unpacked.
- Updated MPC Obs codes to include NEO Surveyor with code C58.
- Updated NAIF ID list to include new designations for about 20-30 comets.
- Renamed underlying `_rust` binary to `_core`.
- Moved NAIF ID list to a dedicated CSV file which is read during compilation.
- Rename `population.diameter` to `population.power_law`.

### Removed

- Removed main belt construction tools, out of scope and not accurate enough.
- Removed redundant MPC name resolver function.


## [0.2.0]  - 2024-3-16

Initial Release


- Initial version release of NEOSpy!
- NEOSpy's primary goal is to enable simulations of NEOs, however it also supports any 
  asteroid or comet. Included are n-body orbit propagation tools, thermal and optical
  modeling, tools for computing what minor planets can be seen by an observer.
  Along with many helpful interfaces to web tools such as JPL Horizons or IPAC's IRSA.


[Unreleased]: https://github.com/IPAC-SW/neospy/tree/main
[0.2.3]: https://github.com/IPAC-SW/neospy/releases/tag/v0.2.3
[0.2.2]: https://github.com/IPAC-SW/neospy/releases/tag/v0.2.2
[0.2.1]: https://github.com/IPAC-SW/neospy/releases/tag/v0.2.1
[0.2.0]: https://github.com/IPAC-SW/neospy/releases/tag/v0.2.0
