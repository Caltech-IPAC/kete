use itertools::Itertools;
use nalgebra::Vector3;
use neospy_core::state::State;
use neospy_core::{constants, propagation};
use pyo3::{exceptions, PyErr};
use pyo3::{pyfunction, PyResult};
use rayon::prelude::*;

use crate::state::PyState;
use crate::vector::Vector;

/// Solve kepler's equation for the Eccentric Anomaly.
///
/// Parameters
/// ----------
/// ecc :
///     Eccentricity, must be non-negative.
/// mean_anom :
///     Mean Anomaly between 0 and 2*pi.
/// peri_dist :
///     Perihelion distance in AU.
#[pyfunction]
#[pyo3(name = "compute_eccentric_anomaly")]
pub fn compute_eccentric_anomaly_py(
    ecc: Vec<f64>,
    mean_anom: Vec<f64>,
    peri_dist: Vec<f64>,
) -> PyResult<Vec<f64>> {
    if ecc.len() != mean_anom.len() || ecc.len() != peri_dist.len() {
        return Err(PyErr::new::<exceptions::PyValueError, _>(
            "Input lengths must all match.",
        ));
    }
    Ok(ecc
        .iter()
        .zip(mean_anom)
        .zip(peri_dist)
        .collect_vec()
        .par_iter()
        .map(|((e, anom), peri)| {
            propagation::compute_eccentric_anomaly(**e, *anom, *peri).unwrap_or(f64::NAN)
        })
        .collect())
}

/// Given a collection of states, and an amount of time (days). Move the states
/// in time using two body mechanics by the desired `dt`.
///
/// This is a multi-core operation.
///
/// Parameters
/// ----------
/// state :
///     List of states, which are in units of AU from the Sun and velocity is in AU/Day.
/// jd :
///     Time to integrate to in JD days with TDB scaling.
/// sun2obj :
///     Position of the observer in AU.
#[pyfunction]
#[pyo3(name = "propagate_two_body")]
pub fn propagation_kepler_py(
    states: Vec<PyState>,
    jd: f64,
    sun2obs: Option<Vector>,
) -> Vec<PyState> {
    states
        .par_iter()
        .map(|state| {
            let center = state.center_id();

            let Some(state) = state.change_center(10).ok() else {
                return State::new_nan(state.0.desig.clone(), jd, state.0.frame, center).into_frame();
            };

            let Some(mut new_state) = propagation::propagate_two_body(&state.0, jd).ok() else {
                return State::new_nan(state.0.desig.clone(), jd, state.0.frame, center).into_frame();
            };

            if let Some(sun2obs) = &sun2obs {
                let sun2obs = Vector3::<f64>::from(sun2obs.raw);
                let delay =
                    -(Vector3::from(new_state.pos) - sun2obs).norm() / constants::C_AU_PER_DAY;

                new_state = match propagation::propagate_two_body(&new_state, new_state.jd + delay)
                {
                    Ok(state) => state,
                    Err(_) => {
                        return State::new_nan(state.0.desig.clone(), jd, state.0.frame, center)
                            .into_frame()
                    }
                };
            }
            PyState(new_state).change_center(center).unwrap_or(
                State::new_nan(state.0.desig.clone(), jd, state.0.frame, state.0.center_id).into_frame(),
            )
        })
        .collect()
}
