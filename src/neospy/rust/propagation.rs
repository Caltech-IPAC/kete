use itertools::Itertools;
use neospy_core::{
    errors::Error,
    propagation::{self, NonGravModel},
    spice::{self, get_spk_singleton},
    state::State,
    time::{scales::TDB, Time},
};
use pyo3::{pyfunction, PyErr, PyResult, Python};
use rayon::prelude::*;

use crate::state::PyState;
use crate::{nongrav::PyNonGravModel, time::PyTime};

/// Propagate the provided :class:`~neospy.State` using N body mechanics to the
/// specified times, no approximations are made, this can be very CPU intensive.
///
/// This does not compute light delay, however it does include corrections for general
/// relativity due to the Sun.
///
/// Parameters
/// ----------
/// states:
///     The initial states, this is a list of multiple State objects.
/// jd:
///     A JD to propagate the initial states to.
/// include_asteroids:
///     If this is true, the computation will include the largest 5 asteroids.
///     The asteroids are: Ceres, Pallas, Interamnia, Hygiea, and Vesta.
/// non_gravs:
///     A list of non-gravitational terms for each object. If provided, then every
///     object must have an associated :class:`~NonGravModel` or `None`.
/// suppress_errors:
///     If True, errors during propagation will return NaN for the relevant state
///     vectors, but propagation will continue.
///
/// Returns
/// -------
/// Iterable
///     A :class:`~neospy.State` at the new time.
#[pyfunction]
#[pyo3(name = "propagate_n_body", signature = (states, jd, include_asteroids=false,
    non_gravs=None, suppress_errors=true))]
pub fn propagation_n_body_spk_py(
    py: Python<'_>,
    states: Vec<PyState>,
    jd: PyTime,
    include_asteroids: bool,
    non_gravs: Option<Vec<Option<PyNonGravModel>>>,
    suppress_errors: bool,
) -> PyResult<Vec<PyState>> {
    let states: Vec<State> = states.into_iter().map(|x| x.0).collect();
    let non_gravs = non_gravs.unwrap_or(vec![None; states.len()]);

    if states.len() != non_gravs.len() {
        Err(Error::ValueError(
            "non_gravs must be the same length as states.".into(),
        ))?;
    }

    let mut res: Vec<PyState> = Vec::new();
    let jd = jd.jd();

    // propagation is broken into chunks of 1000 states, every time a chunk is completed
    // python is checked for signals. This allows keyboard interrupts to be caught
    // and the process interrupted.

    for chunk in states
        .into_iter()
        .zip(non_gravs.into_iter())
        .collect_vec()
        .chunks(1000)
    {
        py.check_signals()?;

        let mut proc_chunk = chunk
            .to_owned()
            .into_par_iter()
            .map(|(mut state, model)| {
                let model = model.map(|x| x.0);
                let spk = get_spk_singleton().try_read().unwrap();
                let center = state.center_id;
                if let Err(e) = spk.try_change_center(&mut state, 0) {
                    if !suppress_errors {
                        Err(e)?;
                    };
                    return Ok::<PyState, PyErr>(
                        State::new_nan(state.desig, jd, state.frame, center).into(),
                    );
                };

                // if the input has a NAN in it, skip the propagation entirely and return
                // the nans.
                if state.pos.iter().any(|x| !x.is_finite())
                    || state.vel.iter().any(|x| !x.is_finite())
                {
                    if !suppress_errors {
                        Err(Error::ValueError("Input state contains NaNs.".into()))?;
                    };
                    return Ok(State::new_nan(state.desig, jd, state.frame, center).into());
                }
                let desig = state.desig.clone();
                let frame = state.frame;
                match propagation::propagate_n_body_spk(state, jd, include_asteroids, model) {
                    Ok(mut state) => {
                        if let Err(er) = spk.try_change_center(&mut state, center) {
                            if !suppress_errors {
                                Err(er)?;
                            }
                            return Ok(State::new_nan(desig, jd, frame, center).into());
                        };
                        Ok(state.into())
                    }
                    Err(er) => {
                        if !suppress_errors {
                            Err(er)?
                        } else {
                            if let neospy_core::errors::Error::Impact(id, time) = er {
                                let time_full: Time<TDB> = Time::new(time);
                                eprintln!(
                                    "Impact detected between {:?} <-> {} at time {} ({})",
                                    desig,
                                    spice::try_name_from_id(id).unwrap_or(id.to_string()),
                                    time,
                                    time_full.utc().to_iso().unwrap()
                                );
                            };
                            Ok(State::new_nan(desig, jd, frame, center).into())
                        }
                    }
                }
            })
            .collect::<PyResult<Vec<_>>>()?;
        res.append(&mut proc_chunk);
    }

    Ok(res)
}

/// Propagate the provided :class:`~neospy.State` using N body mechanics to the
/// specified times, no approximations are made, this can be very CPU intensive.
///
/// This does not compute light delay, however it does include corrections for general
/// relativity due to the Sun.
///
/// Parameters
/// ----------
/// states:
///     The initial states, this is a list of multiple State objects.
/// jd:
///     A JD to propagate the initial states to.
/// include_asteroids:
///     If this is true, the computation will include the largest 5 asteroids.
///     The asteroids are: Ceres, Pallas, Interamnia, Hygiea, and Vesta.
/// non_gravs:
///     A list of non-gravitational terms for each object. If provided, then every
///     object must have an associated :class:`~NonGravModel` or `None`.
/// suppress_errors:
///     If True, errors during propagation will return NaN for the relevant state
///     vectors, but propagation will continue.
///
/// Returns
/// -------
/// Iterable
///     A :class:`~neospy.State` at the new time.
#[pyfunction]
#[pyo3(name = "propagate_n_body_vec", signature = (states, jd_final, non_gravs=None))]
pub fn propagation_n_body_py(
    states: Vec<PyState>,
    jd_final: PyTime,
    non_gravs: Option<Vec<PyNonGravModel>>,
) -> PyResult<Vec<PyState>> {
    let states: Vec<State> = states.into_iter().map(|x| x.0).collect();
    let non_gravs: Option<Vec<NonGravModel>> =
        non_gravs.map(|x| x.into_iter().map(|x| x.0).collect());

    let jd = jd_final.jd();
    let res = propagation::propagate_n_body_vec(states, jd, non_gravs)
        .map(|x| x.into_iter().map(PyState::from).collect::<Vec<_>>())?;
    Ok(res)
}
