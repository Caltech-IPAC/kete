# Kete

![Github Actions](https://github.com/IPAC-SW/kete/actions/workflows/test-lint.yml/badge.svg?branch=main)
[![DOI](https://zenodo.org/badge/787588564.svg)](https://zenodo.org/badge/latestdoi/787588564)

The repository for Kete, Solar System Survey Simulation Software.
'Kete' comes from ancient greek mythology, meaning sea monsters, and is the root word
for Cetaceans (Whales).

## Note:

Original maintainer (Dar Dahlen) is beginning a PhD and will no longer be performing
maintenance, as he is no longer affiliated with CalTech.

## Introduction

The kete tools are intended to enable the simulation of all-sky surveys of minor
planets. This includes multi-body physics orbital dynamics, thermal and optical modeling
of the objects, as well as field of view and light delay corrections. These tools in
conjunction with the Minor Planet Centers (MPC) database of known asteroids can be used
to not only plan surveys but can also be used to predict what objects are visible for
existing or past surveys.

The primary goal for kete is to enable a set of tools that can operate on the entire
MPC catalog at once, without having to do queries on specific objects. It has been
used to simulate over 10 years of survey time for the NEO Surveyor mission using 10
million main-belt and near-Earth asteroids.

[Documentation](https://caltech-ipac.github.io/kete/)  
   - [Examples](https://caltech-ipac.github.io/kete/auto_examples/index.html)
   - [Tutorials](https://caltech-ipac.github.io/kete/tutorials/index.html)


https://github.com/user-attachments/assets/a48491d8-9c15-4659-9022-1767a3aa1e94

Here is a simulation of what the ZTF survey would observe during the entirety of 2023.
This is every position of every numbered asteroid, along with a calculation of the
expected V-band magnitudes. If the expected magnitude is less than ZTF's reported
magnitude limit for the specific frame, then the object will flash light grey.

This took about 50 minutes on a desktop computer to compute, and about 40 minutes
to generate the movie.




## Installation

Kete may be installed using pip:

``` bash
pip install kete
```

## Installation - From Source

If kete is built from source, the rust compiler must be installed. Installation
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

Kete uses the ICRF Reference frame as the base coordinate frame, with units of AU,
with time in JD with Barycentric Dynamical Time (TDB) scaling. Internally this frame
converted to an Ecliptic coordinate system which is defined by the Obliquity Angle
definition used by JPL Horizons, which is the defined IAU76/80 model in the J2000 frame.

      - https://en.wikipedia.org/wiki/Axial_tilt#Short_term
      - https://ssd.jpl.nasa.gov/horizons/manual.html#defs

Both JPL Horizons and the Minor Planet Center (MPC) use this coordinate frame, which is
essentially equivalent to J2000 Ecliptic coordinates. Conversion tools are available in
kete which enable conversion to the Equatorial frame and to various flavors of time.

### Cache directory

Many operations in kete result in downloading various files. These files are cached
automatically, the directory where this cache is stored may be set by setting the
environment variable `KETE_CACHE_DIR`. The default directory is `~/.kete/`.
``` bash
export KETE_CACHE_DIR="~/.kete/"
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
kete directory.
``` bash
(cd docs && make html && open html/index.html&)
```
Once this has completed running, open the file `kete/docs/html/index.html` for access
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
pytest --cov-report term-missing --cov=kete   
```

Another coverage report type is HTML, this will generate a folder called `htmlcov`
in the directory where the command was run, then you can open the `htmlcov/index.html`
file. This is a user-friendly website representation of the code coverage.
``` bash
pytest --cov-report html --cov=kete   
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

There are a test suite of micro-benchmarks in the rust backend of kete. These require
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
