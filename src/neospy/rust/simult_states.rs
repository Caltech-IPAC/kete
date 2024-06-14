use neospy_core::io::FileIO;
use neospy_core::simult_states;
use neospy_core::state::State;
use pyo3::exceptions;
use pyo3::prelude::*;
use pyo3::{pyclass, pymethods, PyResult};

use crate::{fovs::AllowedFOV, frame::PyFrames, state::PyState};

/// Polymorphic support
#[derive(FromPyObject)]
pub enum SimulStateLike {
    Vec(Vec<PyState>),
    Simul(Py<PySimultaneousStates>),
}

impl SimulStateLike {
    /// Convert state-like object into a simultaneous state.
    pub fn into_simul_state(self, py: Python<'_>) -> PyResult<simult_states::SimultaneousStates> {
        Ok(match self {
            SimulStateLike::Vec(state_vec) => simult_states::SimultaneousStates::new_exact(
                state_vec.into_iter().map(|x| x.0).collect(),
                None,
            )?,

            SimulStateLike::Simul(p) => *p.extract::<PySimultaneousStates>(py)?.0,
        })
    }

    /// Convert a collection of state-like objects into a collection of exact states.
    pub fn into_states(self, py: Python<'_>) -> PyResult<Vec<State>> {
        Ok(match self {
            SimulStateLike::Vec(state_vec) => state_vec.into_iter().map(|x| x.0).collect(),

            SimulStateLike::Simul(p) => p.extract::<PySimultaneousStates>(py)?.0.states,
        })
    }
}

/// Representation of a collection of [`State`] at a single point in time.
///
/// The main value in this is that also includes an optional Field of View.
/// If the FOV is provided, it is implied that the states which are present
/// in this file were objects seen by the FOV.
///
#[pyclass(module = "neospy", frozen, sequence, name = "SimultaneousStates")]
#[derive(Clone, Debug)]
pub struct PySimultaneousStates(pub Box<simult_states::SimultaneousStates>);

impl From<simult_states::SimultaneousStates> for PySimultaneousStates {
    fn from(value: simult_states::SimultaneousStates) -> Self {
        Self(Box::new(value))
    }
}

#[pymethods]
impl PySimultaneousStates {
    /// Create a new collection of States at a specific time.
    ///
    ///
    /// Parameters
    /// ----------
    /// states :
    ///     List of States to include.
    /// fov :
    ///     An optional FOV, if this is provided it is expected that the states provided
    ///     are what have been seen by this FOV. This is not checked.
    #[new]
    pub fn new(states: Vec<PyState>, fov: Option<AllowedFOV>) -> PyResult<Self> {
        let states: Vec<_> = states.into_iter().map(|x| x.0).collect();
        let fov = fov.map(|x| x.unwrap());
        Ok(simult_states::SimultaneousStates::new_exact(states, fov)
            .map(|x| PySimultaneousStates(Box::new(x)))?)
    }

    /// The FOV if it exists.
    #[getter]
    pub fn fov(&self) -> Option<AllowedFOV> {
        self.0.fov.clone().map(|x| x.into())
    }

    /// States contained within.
    #[getter]
    pub fn states(&self) -> Vec<PyState> {
        self.0.states.iter().map(|x| x.clone().into()).collect()
    }

    /// The time of the simultaneous states.
    #[getter]
    pub fn jd(&self) -> f64 {
        self.0.jd
    }

    /// The reference center NAIF ID for this state.
    #[getter]
    pub fn center_id(&self) -> i32 {
        self.0.center_id
    }

    /// Coordinate Frame.
    #[getter]
    pub fn frame(&self) -> PyFrames {
        self.0.frame.into()
    }

    /// Load a single SimultaneousStates from a file.
    #[staticmethod]
    pub fn load(filename: String) -> PyResult<Self> {
        Ok(PySimultaneousStates(Box::new(
            neospy_core::simult_states::SimultaneousStates::load(filename)?,
        )))
    }

    /// Save a single SimultaneousStates to a file.
    pub fn save(&self, filename: String) -> PyResult<()> {
        self.0.save(filename)?;
        Ok(())
    }

    pub fn __len__(&self) -> usize {
        self.0.states.len()
    }

    pub fn __getitem__(&self, idx: usize) -> PyResult<PyState> {
        if idx >= self.__len__() {
            return Err(PyErr::new::<exceptions::PyIndexError, _>(""));
        }
        Ok(self.0.states[idx].clone().into())
    }

    fn __repr__(&self) -> String {
        let n_states = self.0.states.len();
        let fov_str = match self.fov() {
            None => "None".into(),
            Some(f) => f.__repr__(),
        };
        format!(
            "SimultaneousStates(states=<{} States>, fov={})",
            n_states, fov_str
        )
    }

    /// Save a list to a binary file.
    ///
    /// Note that this saves a list of SimultaneousStates, meaning it is a list of a list of States.
    #[staticmethod]
    #[pyo3(name = "save_list")]
    pub fn py_save_list(vec: Vec<Self>, filename: String) -> PyResult<()> {
        let vec: Vec<_> = vec.into_iter().map(|x| *x.0).collect();
        Ok(neospy_core::simult_states::SimultaneousStates::save_vec(
            &vec, filename,
        )?)
    }

    /// Load a list from a binary file.
    ///
    /// Note that this loads a list of SimultaneousStates, meaning it is a list of a list of States.
    #[staticmethod]
    #[pyo3(name = "load_list")]
    pub fn py_load_list(filename: String) -> PyResult<Vec<Self>> {
        let res = neospy_core::simult_states::SimultaneousStates::load_vec(filename)?;
        Ok(res.into_iter().map(|x| Self(Box::new(x))).collect())
    }
}
