use std::fmt::Debug;

use crate::elements::PyCometElements;
use crate::state::PyState;
use neospy_core::{io::FileIO, prelude};
use pyo3::prelude::*;
use serde::{Deserialize, Serialize};

/// Horizons object properties
/// Physical, orbital, and observational properties of a solar system object as recorded in JPL Horizons.
#[pyclass(frozen, get_all, module = "neospy")]
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct HorizonsCovariance {
    /// Epoch of the covariance matrix fit.
    epoch: f64,

    /// Name and best estimate of the parameters in the fit.
    params: Vec<(String, f64)>,

    /// The covariance matrix, where the order of the array corresponds to the parameters.
    cov_matrix: Vec<Vec<f64>>,
}

impl neospy_core::io::FileIO for HorizonsCovariance {}

#[pymethods]
impl HorizonsCovariance {
    #[new]
    #[allow(clippy::too_many_arguments)]
    pub fn new(epoch: f64, params: Vec<(String, f64)>, cov_matrix: Vec<Vec<f64>>) -> Self {
        Self {
            epoch,
            params,
            cov_matrix,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "HorizonsCovariance(epoch={:?}, params={:?}, cov_matrix={:?})",
            self.epoch, self.params, self.cov_matrix,
        )
    }

    /// Save the covariance matrix to a file.
    #[pyo3(name = "save")]
    pub fn py_save(&self, filename: String) -> PyResult<usize> {
        Ok(self.save(filename)?)
    }

    /// Load a covariance matrix from a file.
    #[staticmethod]
    #[pyo3(name = "load")]
    pub fn py_load(filename: String) -> PyResult<Self> {
        Ok(Self::load(filename)?)
    }
}

/// Horizons object properties
/// Physical, orbital, and observational properties of a solar system object as recorded in JPL Horizons.
#[pyclass(frozen, get_all, module = "neospy")]
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct HorizonsProperties {
    /// The MPC packed designation of the object.
    desig: String,

    /// An optional group name to associate the object with a group.
    group: Option<String>,

    /// The epoch during which the orbital elements listed are accurate, in JD, TDB.
    epoch: Option<f64>,

    /// The eccentricity of the orbit.
    eccentricity: Option<f64>,

    /// The inclination of the orbit in degrees.
    inclination: Option<f64>,

    /// The longitudinal node of the orbit in degrees.
    lon_of_ascending: Option<f64>,

    /// The argument of perihelion in degrees.
    peri_arg: Option<f64>,

    /// The perihelion distance in AU.
    peri_dist: Option<f64>,

    /// The time of perihelion in JD, TDB scaled time.
    peri_time: Option<f64>,

    /// The H magnitude of the object.
    h_mag: Option<f64>,

    /// The visible albedo of the object, between 0 and 1.
    vis_albedo: Option<f64>,

    /// The diameter of the object in km.
    diameter: Option<f64>,

    /// The minimum orbital intersection distance between the object and Earth in AU.
    moid: Option<f64>,

    /// The g parameter of the object.
    g_phase: Option<f64>,

    /// If the object was previously known, this lists the length of time of the
    /// observations of the object in days.
    arc_len: Option<f64>,

    /// Covariance values in the orbit fit.
    covariance: Option<HorizonsCovariance>,
}

impl neospy_core::io::FileIO for HorizonsProperties {}

#[pymethods]
impl HorizonsProperties {
    #[new]
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        desig: String,
        group: Option<String>,
        epoch: Option<f64>,
        eccentricity: Option<f64>,
        inclination: Option<f64>,
        lon_of_ascending: Option<f64>,
        peri_arg: Option<f64>,
        peri_dist: Option<f64>,
        peri_time: Option<f64>,
        h_mag: Option<f64>,
        vis_albedo: Option<f64>,
        diameter: Option<f64>,
        moid: Option<f64>,
        g_phase: Option<f64>,
        arc_len: Option<f64>,
        covariance: Option<HorizonsCovariance>,
    ) -> Self {
        Self {
            desig,
            group,
            vis_albedo,
            diameter,
            moid,
            peri_dist,
            eccentricity,
            inclination,
            lon_of_ascending,
            peri_arg,
            peri_time,
            h_mag,
            g_phase,
            epoch,
            arc_len,
            covariance,
        }
    }

    /// Cometary orbital elements.
    #[getter]
    pub fn elements(&self) -> PyResult<PyCometElements> {
        Ok(PyCometElements(prelude::CometElements {
            desig: prelude::Desig::Name(self.desig.clone()),
            epoch: self
                .epoch
                .ok_or_else(|| prelude::NEOSpyError::ValueError("No Epoch defined".into()))?,
            eccentricity: self.eccentricity.ok_or_else(|| {
                prelude::NEOSpyError::ValueError("No Eccentricity defined".into())
            })?,
            inclination: self
                .inclination
                .ok_or_else(|| prelude::NEOSpyError::ValueError("No Inclination defined".into()))?
                .to_radians(),
            peri_arg: self
                .peri_arg
                .ok_or_else(|| prelude::NEOSpyError::ValueError("No peri_arg defined".into()))?
                .to_radians(),
            peri_dist: self
                .peri_dist
                .ok_or_else(|| prelude::NEOSpyError::ValueError("No peri_dist defined".into()))?,
            peri_time: self
                .peri_time
                .ok_or_else(|| prelude::NEOSpyError::ValueError("No peri_time defined".into()))?,
            lon_of_ascending: self
                .lon_of_ascending
                .ok_or_else(|| {
                    prelude::NEOSpyError::ValueError(
                        "No longitude of ascending node defined".into(),
                    )
                })?
                .to_radians(),
            frame: prelude::Frame::Ecliptic,
        }))
    }

    /// Convert the orbital elements of the object to a State.
    #[getter]
    pub fn state(&self) -> PyResult<PyState> {
        self.elements()?.as_state()
    }

    fn __repr__(&self) -> String {
        fn cleanup<T: Debug>(opt: Option<T>) -> String {
            match opt {
                None => "None".into(),
                Some(val) => format!("{:?}", val),
            }
        }

        let cov = match self.covariance {
            Some(_) => "<present>",
            None => "None",
        };

        format!(
            "HorizonsObject(desig={:?}, group={:}, epoch={:}, eccentricity={:}, inclination={:}, \
            lon_of_ascending={:}, peri_arg={:}, peri_dist={:}, peri_time={:}, h_mag={:}, \
            vis_albedo={:}, diameter={:}, moid={:}, g_phase={:}, arc_len={:}, \
            covariance={:})",
            self.desig,
            cleanup(self.group.clone()),
            cleanup(self.epoch),
            cleanup(self.eccentricity),
            cleanup(self.inclination),
            cleanup(self.lon_of_ascending),
            cleanup(self.peri_arg),
            cleanup(self.peri_dist),
            cleanup(self.peri_time),
            cleanup(self.h_mag),
            cleanup(self.vis_albedo),
            cleanup(self.diameter),
            cleanup(self.moid),
            cleanup(self.g_phase),
            cleanup(self.arc_len),
            cov,
        )
    }

    /// Save the horizons query to a file.
    #[pyo3(name = "save")]
    pub fn py_save(&self, filename: String) -> PyResult<usize> {
        Ok(self.save(filename)?)
    }

    /// Load the horizons query from a file.
    #[staticmethod]
    #[pyo3(name = "load")]
    pub fn py_load(filename: String) -> PyResult<Self> {
        Ok(Self::load(filename)?)
    }
}
