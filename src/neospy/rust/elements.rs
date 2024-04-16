use neospy_core::elements;
use neospy_core::prelude;
use pyo3::{pyclass, pymethods, PyResult};

use crate::state::PyState;

/// Cometary Elements class made accessible to python.
///
/// Angles must be in degrees, distances in AU.
///
/// Parameters
/// ----------
/// desig:
///     The designations of the object.
/// epoch:
///     The epoch time for the orbital elements.
/// eccentricity:
///     The eccentricity of the orbit.
/// inclination:
///     The inclination, must be in degrees.
/// peri_dist:
///     The perihelion distance in AU.
/// peri_arg:
///     The argument of perihelion, must be in degrees.
/// peri_time:
///     The JD time of perihelion.
/// lon_of_ascending:
///     The longitude of ascending node, in degrees.
#[pyclass(module = "neospy", frozen, name = "CometElements")]
#[derive(Clone, Debug)]
pub struct PyCometElements(pub elements::CometElements);

#[pymethods]
impl PyCometElements {
    #[new]
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        desig: String,
        epoch: f64,
        eccentricity: f64,
        inclination: f64,
        peri_dist: f64,
        peri_arg: f64,
        peri_time: f64,
        lon_of_ascending: f64,
    ) -> Self {
        Self(elements::CometElements {
            desig: prelude::Desig::Name(desig),
            frame: prelude::Frame::Ecliptic,
            epoch,
            eccentricity,
            inclination: inclination.to_radians(),
            lon_of_ascending: lon_of_ascending.to_radians(),
            peri_time,
            peri_arg: peri_arg.to_radians(),
            peri_dist,
        })
    }

    #[staticmethod]
    pub fn from_state(state: &PyState) -> Self {
        Self(elements::CometElements::from_state(&state.0))
    }

    #[getter]
    pub fn epoch(&self) -> f64 {
        self.0.epoch
    }

    #[getter]
    pub fn desig(&self) -> String {
        match &self.0.desig {
            prelude::Desig::Name(s) => s.clone(),
            prelude::Desig::Naif(s) => {
                neospy_core::spice::try_name_from_id(*s).unwrap_or(s.to_string())
            }
            prelude::Desig::Perm(s) => format!("{:?}", s),
            prelude::Desig::Prov(s) => s.clone(),
            prelude::Desig::Empty => "None".into(),
        }
    }

    #[getter]
    pub fn eccentricity(&self) -> f64 {
        self.0.eccentricity
    }

    #[getter]
    pub fn inclination(&self) -> f64 {
        self.0.inclination.to_degrees()
    }

    #[getter]
    pub fn lon_of_ascending(&self) -> f64 {
        self.0.lon_of_ascending.to_degrees()
    }

    #[getter]
    pub fn peri_time(&self) -> f64 {
        self.0.peri_time
    }

    #[getter]
    pub fn peri_arg(&self) -> f64 {
        self.0.peri_arg.to_degrees()
    }

    #[getter]
    pub fn peri_dist(&self) -> f64 {
        self.0.peri_dist
    }

    #[getter]
    pub fn semi_major(&self) -> f64 {
        self.0.semi_major()
    }

    #[getter]
    pub fn mean_motion(&self) -> f64 {
        self.0.mean_motion().to_degrees()
    }

    #[getter]
    pub fn orbital_period(&self) -> f64 {
        self.0.orbital_period()
    }

    #[getter]
    pub fn as_state(&self) -> PyResult<PyState> {
        Ok(self.0.try_to_state()?.into())
    }

    #[getter]
    pub fn eccentric_anomaly(&self) -> PyResult<f64> {
        Ok(self.0.eccentric_anomaly().map(|x| x.to_degrees())?)
    }

    #[getter]
    pub fn mean_anomaly(&self) -> f64 {
        self.0.mean_anomaly().to_degrees()
    }

    #[getter]
    pub fn true_anomaly(&self) -> PyResult<f64> {
        Ok(self.0.true_anomaly().map(|x| x.to_degrees())?)
    }

    fn __repr__(&self) -> String {
        format!("CometElements(desig={:?}, epoch={}, eccentricity={}, inclination={}, lon_of_ascending={}, peri_time={}, peri_arg={}, peri_dist={})", self.desig(), self.epoch(), self.eccentricity(), self.inclination(), self.lon_of_ascending(), self.peri_time(), self.peri_arg(), self.peri_dist())
    }
}
