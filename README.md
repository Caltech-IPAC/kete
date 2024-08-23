# Apohele

![Github Actions](https://github.com/IPAC-SW/apohele/actions/workflows/test-lint.yml/badge.svg?branch=main)

The repository for Apohele, tools for Solar System object Surveyor Simulation.

Apohele comes from the Hawaiian word for 'Orbit', Apo - Circle, Hele - to go. The first
detected Atira asteroid was 1998 DK36, which was proposed to be called an Apohele
asteroid, before Atira was discovered. 1998 DK36 remains a lost asteroid since its
initial discovery, though is very likely to be observed by the NEO Surveyor, most
likely either during the first year of operation, or in the first half of 2032
(the uncertainty in the measured orbit make this difficult to predict).

## Introduction

The Apohele tools are intended to enable the simulation of all-sky surveys of minor
planets. This includes multi-body physics orbital dynamics, thermal and optical modeling
of the objects, as well as field of view and light delay corrections. These tools in
conjunction with the Minor Planet Centers (MPC) database of known asteroids can be used
to not only plan surveys but can also be used to predict what objects are visible for
existing or past surveys.

The primary goal for apohele is to enable a set of tools that can operate on the entire
MPC catalog at once, without having to do queries on specific objects. It has been
used to simulate over 10 years of survey time for the NEO Surveyor mission using 10
million main-belt and near-Earth asteroids.

## Installation

If apohele is built from source, the rust compiler must be installed. Installation
instructions may be found here: 

https://www.rust-lang.org/learn/get-started

Ensure that your Python is up to date, this code runs on Python 3.9+.
``` bash
python --version
```

Ensure that your pip is up to date, this should be at least version `22.0.0`.
``` bash
pip --version
```

This can be updated using:
``` bash
python -m pip install "pip>=22.0.0" --upgrade
pip install setuptools --upgrade
```

### Units and Reference Frame

apohele uses the ICRF Reference frame as the base coordinate frame, with units of AU,
with time in JD with Barycentric Dynamical Time (TDB) scaling. Internally this frame
converted to an Ecliptic coordinate system which is defined by the Obliquity Angle
definition used by JPL Horizons, which is the defined IAU76/80 model in the J2000 frame.

      - https://en.wikipedia.org/wiki/Axial_tilt#Short_term
      - https://ssd.jpl.nasa.gov/horizons/manual.html#defs

Both JPL Horizons and the Minor Planet Center (MPC) use this coordinate frame, which is
essentially equivalent to J2000 Ecliptic coordinates. Conversion tools are available in
apohele which enable conversion to the Equatorial frame and to various flavors of time.

### Cache directory

Many operations in apohele result in downloading various files. These files are cached
automatically, the directory where this cache is stored may be set by setting the
environment variable `apohele_CACHE_DIR`. The default directory is `~/.apohele/`.
``` bash
export apohele_CACHE_DIR="~/.apohele/"
```

### Development
If you plan on doing development, it is recommended to install with the following:
``` bash
pip install '.[dev]'
```
The `[dev]` in that line has pip install a number of optional dependencies which
are useful for development. Including pytest and documentation tools.

### Building Documentation

In order for documentation to be built, some additional Python libraries are needed.
These can be installed with:
``` bash
pip install sphinx sphinx_gallery autodoc
```
After this has been installed, the documentation can be built by running inside the
apohele directory.
``` bash
(cd docs && make html && open html/index.html&)
```
Once this has completed running, open the file `apohele/docs/html/index.html` for access
to the HTML documentation.

To clean the previous docs build:
``` bash
(cd docs && make clean)
```

Documentation tests may be run with:
``` bash
(cd docs && make doctest)
```

### Running tests

Running tests require that the `pytest` and `pytest-cov` packages be installed.

Open a terminal in the base of this folder and run the following command:
``` bash
pytest --cov-report term-missing --cov=apohele   
```

Another coverage report type is HTML, this will generate a folder called `htmlcov`
in the directory where the command was run, then you can open the `htmlcov/index.html`
file. This is a user-friendly website representation of the code coverage.
``` bash
pytest --cov-report html --cov=apohele   
```

### Running Tutorials

Tutorials are computationally expensive examples which are more indicative of typical
expected use. Since these examples are so expensive to run, they are not run unless
manually performed. A convenience python script has been provided to do just this.

``` bash
cd docs
python utils.py
```

### Running Benchmarks

There are a test suite of micro-benchmarks in the rust backend of apohele. These require
`gnuplot` to be installed, and may be run using the following command:

``` bash
cargo bench
open target/criterion/report/index.html
```

Additionally, Flamegraphs may be produced using the following:

``` bash
cargo bench --bench propagation -- --profile-time=5
cargo bench --bench spice -- --profile-time=5
cargo bench --bench thermal -- --profile-time=5
```

These flamegraphs will be put in `target/criterion/*/profile/flamegraph.svg`. Opening
these files in a web browser will show what functions are being used during the bench.
