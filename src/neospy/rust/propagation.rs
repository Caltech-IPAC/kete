use itertools::Itertools;
use neospy_core::{
    errors::NEOSpyError,
    propagation,
    spice::{self, get_spk_singleton},
    state::State,
};
use pyo3::{pyfunction, PyErr, PyResult, Python};
use rayon::prelude::*;

use crate::simult_states::SimulStateLike;
use crate::state::PyState;

/// python wrapper over the propagation_n_body_spk function.
#[pyfunction]
#[pyo3(name = "propagate_n_body_spk")]
pub fn propagation_n_body_spk_py(
    py: Python<'_>,
    states: SimulStateLike,
    jd_final: f64,
    include_asteroids: bool,
    a_terms: Vec<Option<(f64, f64, f64, bool)>>,
    suppress_errors: Option<bool>,
) -> PyResult<Vec<PyState>> {
    let suppress_errors = suppress_errors.unwrap_or(true);
    let states = states.into_states(py)?;

    let mut res: Vec<PyState> = Vec::new();

    // propagation is broken into chunks of 1000 states, every time a chunk is completed
    // python is checked for signals. This allows keyboard interrupts to be caught
    // and the process interrupted.

    for chunk in states.into_iter().zip(a_terms).collect_vec().chunks(1000) {
        py.check_signals()?;

        let mut proc_chunk = chunk
            .to_vec()
            .into_par_iter()
            .map(|(mut state, a_term)| {
                let a_term = a_term.map(|(a1, a2, a3, comet)| {
                    if comet {
                        propagation::NonGravModel::new_jpl_comet_default(a1, a2, a3)
                    } else {
                        propagation::NonGravModel::new_dust(a1, a2)
                    }
                });
                let spk = get_spk_singleton().try_read().unwrap();
                let center = state.center_id;
                if let Err(e) = spk.try_change_center(&mut state, 0) {
                    if !suppress_errors {
                        Err(e)?;
                    };
                    return Ok::<PyState, PyErr>(
                        State::new_nan(state.desig.clone(), jd_final, state.frame, center).into(),
                    );
                };

                // if the input has a NAN in it, skip the propagation entirely and return
                // the nans.
                if state.pos.iter().any(|x| !x.is_finite())
                    || state.vel.iter().any(|x| !x.is_finite())
                {
                    if !suppress_errors {
                        Err(NEOSpyError::ValueError("Input state contains NaNs.".into()))?;
                    };
                    return Ok(
                        State::new_nan(state.desig.clone(), jd_final, state.frame, center).into(),
                    );
                }
                let desig = state.desig.clone();
                let frame = state.frame;
                match propagation::propagate_n_body_spk(state, jd_final, include_asteroids, a_term)
                {
                    Ok(mut state) => {
                        if let Err(er) = spk.try_change_center(&mut state, center) {
                            if !suppress_errors {
                                Err(er)?;
                            }
                            return Ok(State::new_nan(desig, jd_final, frame, center).into());
                        };
                        Ok(state.into())
                    }
                    Err(er) => {
                        if !suppress_errors {
                            Err(er)?
                        } else {
                            if let neospy_core::errors::NEOSpyError::Impact(id, time) = er {
                                eprintln!(
                                    "Impact detected between {:?} <-> {} at time {}",
                                    desig,
                                    spice::try_name_from_id(id).unwrap_or(id.to_string()),
                                    time
                                );
                            };
                            Ok(State::new_nan(desig, jd_final, frame, center).into())
                        }
                    }
                }
            })
            .collect::<PyResult<Vec<_>>>()?;
        res.append(&mut proc_chunk);
    }

    Ok(res)
}
