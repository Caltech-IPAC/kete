use super::*;
use neospy_core::propagation::propagate_n_body_spk;
use pyo3::prelude::*;
use rayon::prelude::*;

use crate::{simult_states::{PySimultaneousStates, SimulStateLike}, vector::{Vector, VectorLike}};

#[pyfunction]
#[pyo3(name = "fov_checks")]
pub fn fov_checks_py(
    py: Python<'_>,
    obj_state: SimulStateLike,
    fovs: FOVListLike,
    dt_limit: f64,
) -> PyResult<Vec<PySimultaneousStates>> {
    let fovs = fovs.into_sorted_vec_fov();

    // This is only here for a check to verify the states are valid
    let pop = obj_state.into_simul_state(py)?;

    let mut jd = pop.jd;
    let mut big_jd = jd;
    let mut states = pop.states;
    let mut big_step_states = states.clone();

    let mut visible = Vec::new();
    for fov in fovs.into_iter() {
        // Take large steps which are 10x the smaller steps, this helps long term numerical stability
        if (fov.observer().jd - big_jd).abs() >= dt_limit * 50.0 {
            big_jd = fov.observer().jd;
            big_step_states = big_step_states
                .into_par_iter()
                .filter_map(|state| propagate_n_body_spk(state, jd, true, None).ok())
                .collect();
        };
        // Take small steps based off of the large steps.
        if (fov.observer().jd - jd).abs() >= dt_limit {
            if (jd - big_jd) >= dt_limit * 25.0 {
                states = big_step_states.clone();
            }
            jd = fov.observer().jd;
            states = states
                .into_par_iter()
                .filter_map(|state| propagate_n_body_spk(state, jd, true, None).ok())
                .collect();
        };

        let vis: Vec<PySimultaneousStates> = fov
            .check_visible(&states, dt_limit)
            .into_iter()
            .filter_map(|pop| pop.map(|p| PySimultaneousStates(Box::new(p))))
            .collect();
        if !vis.is_empty() {
            visible.push(vis);
        }

        py.check_signals()?;
    }
    Ok(visible.into_iter().flatten().collect())
}

#[pyfunction]
#[pyo3(name = "fov_spk_checks")]
pub fn fov_spk_checks_py(obj_ids: Vec<isize>, fovs: FOVListLike) -> Vec<PySimultaneousStates> {
    let fovs = fovs.into_sorted_vec_fov();

    fovs.into_par_iter()
        .filter_map(|fov| {
            let vis: Vec<_> = fov
                .check_spks(&obj_ids)
                .into_iter()
                .filter_map(|pop| pop.map(|p| PySimultaneousStates(Box::new(p))))
                .collect();
            match vis.is_empty() {
                true => None,
                false => Some(vis),
            }
        })
        .flatten()
        .collect()
}


/// Check if a list of static sky positions are present in the given Field of View list.
/// 
/// This returns a list of tuples, where the first entry in the tuple is the vector of
/// all of the points in the provided FOV, and the second entry is the original FOV.
/// 
/// Parameters
/// ----------
/// pos :
///     Collection of Vectors defining sky positions from the point of view of the observer.
///     These vectors are automatically converted to the Ecliptic frame, results will be
///     returned in that frame as well.
/// fovs :
///     Collection of Field of Views to check.
#[pyfunction]
#[pyo3(name = "fov_static_checks")]
pub fn fov_static_checks_py(pos: Vec<VectorLike>, fovs: FOVListLike) -> Vec<(Vec<Vector>, AllowedFOV)> {
    let fovs = fovs.into_sorted_vec_fov();
    let pos:Vec<_> = pos.into_iter().map(|p| p.into_vec(crate::frame::PyFrames::Ecliptic)).collect();

    fovs.into_par_iter()
        .filter_map(|fov| {
            let vis: Vec<_> = fov
                .check_statics(&pos)
                .into_iter()
                .filter_map(|pop| {
                    pop.map(|(p_vec, fov)| {
                        let p_vec = p_vec.into_iter().map(|p| Vector::new(p.into(), crate::frame::PyFrames::Ecliptic)).collect();
                        (p_vec, fov.into())
                    })
            })
                .collect();
            match vis.is_empty() {
                true => None,
                false => Some(vis),
            }
        })
        .flatten()
        .collect()
}
