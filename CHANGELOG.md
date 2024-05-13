# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
[0.2.1]: https://github.com/IPAC-SW/neospy/releases/tag/v0.2.1
[0.2.0]: https://github.com/IPAC-SW/neospy/releases/tag/v0.2.0
