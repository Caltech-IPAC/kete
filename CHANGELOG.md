# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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


[0.2.2]: https://github.com/IPAC-SW/neospy/tree/main
[0.2.1]: https://github.com/IPAC-SW/neospy/releases/tag/v0.2.1
[0.2.0]: https://github.com/IPAC-SW/neospy/releases/tag/v0.2.0
