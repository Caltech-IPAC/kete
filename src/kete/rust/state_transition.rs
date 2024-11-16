//! State Transition matrix computation
use kete_core::frames::Equatorial;
use kete_core::propagation::compute_state_transition;
use kete_core::state::State;
use pyo3::pyfunction;

use crate::state::PyState;
use crate::time::PyTime;

/// Compute an approximate STM
#[pyfunction]
#[pyo3(name = "compute_stm")]
pub fn compute_stm_py(
    state: PyState,
    jd_end: PyTime,
    central_mass: f64,
) -> ([[f64; 3]; 2], [[f64; 6]; 6]) {
    let mut state: State<Equatorial> = state.0;

    let (final_state, stm) = compute_state_transition(&mut state, jd_end.jd(), central_mass);

    (final_state, stm.into())
}
