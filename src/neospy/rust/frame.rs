use neospy_core::{frames::*, state::State};
use pyo3::prelude::*;

/// Defined inertial frames supported by the python side of NEOSpy
#[pyclass(frozen, name = "Frames", module = "neospy")]
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PyFrames {
    Ecliptic,
    Equatorial,
    Galactic,
    FK4,
}

#[derive(Clone, Debug)]
pub enum FrameState {
    Ecliptic(State<Ecliptic>),
    Equatorial(State<Equatorial>),
    Galactic(State<Galactic>),
    FK4(State<FK4>),
}

/// Compute a ECEF position from WCS84 Geodetic latitude/longitude/height.
///
/// This returns the X/Y/Z coordinates in km from the geocenter of the Earth.
///
/// Parameters
/// ----------
/// lat :
///     Latitude in degrees.
/// lon :
///     Longitude in degrees.
/// h :
///     Height above the surface of the Earth in km.
#[pyfunction]
#[pyo3(name = "wgs_lat_lon_to_ecef")]
pub fn wgs_lat_lon_to_ecef(lat: f64, lon: f64, h: f64) -> (f64, f64, f64) {
    geodetic_lat_lon_to_ecef(lat.to_radians(), lon.to_radians(), h)
}
