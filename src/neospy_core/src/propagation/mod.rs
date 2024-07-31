//! # Propagation
//! The motion of objects (represented by a [`State`]) as a function of time.
//! There are multiple levels of precision available, each with different pros/cons
//! (usually performance related).

use crate::constants::{MASSES, PLANETS, SIMPLE_PLANETS};
use crate::errors::Error;
use crate::frames::Frame;
use crate::prelude::{Desig, NeosResult};
use crate::spice::get_spk_singleton;
use crate::state::State;
use nalgebra::{DVector, Vector3};

mod acceleration;
mod kepler;
mod nongrav;
mod radau;
mod runge_kutta;
mod state_transition;

// expose the public methods in spk to the outside world.
pub use acceleration::*;
pub use kepler::*;
pub use nongrav::*;
pub use radau::*;
pub use runge_kutta::*;
pub use state_transition::*;

/// Using the Radau 15th order integrator, integrate the position and velocity of an
/// object assuming two body mechanics with a central object located at 0, 0 with the
/// mass of the sun.
///
/// This is primarily intended for testing the numerical integrator, it is strongly
/// recommended to use the Kepler analytic equation to do actual two body calculations.
pub fn propagate_two_body_radau(dt: f64, pos: &[f64; 3], vel: &[f64; 3]) -> ([f64; 3], [f64; 3]) {
    let res = RadauIntegrator::integrate(
        &central_accel,
        Vector3::from(*pos),
        Vector3::from(*vel),
        0.0,
        dt,
        CentralAccelMeta::default(),
    )
    .unwrap();
    (res.0.into(), res.1.into())
}

/// Propagate the state of an object, only considering linear motion.
///
/// This is a very poor approximation over more than a few minutes/hours, however it
/// is very fast.
pub fn propagate_linear(state: &State, jd_final: f64) -> NeosResult<State> {
    let dt = jd_final - state.jd;
    let mut pos: Vector3<f64> = state.pos.into();
    pos.iter_mut()
        .zip(&state.vel)
        .for_each(|(p, v)| *p += v * dt);

    Ok(State::new(
        state.desig.to_owned(),
        jd_final,
        pos,
        state.vel.into(),
        state.frame,
        state.center_id,
    ))
}

/// Propagate an object using full N-Body physics with the Radau 15th order integrator.
pub fn propagate_n_body_spk(
    mut state: State,
    jd_final: f64,
    include_extended: bool,
    non_grav_model: Option<NonGravModel>,
) -> NeosResult<State> {
    let center = state.center_id;
    let frame = state.frame;
    let spk = get_spk_singleton().try_read().unwrap();
    spk.try_change_center(&mut state, 0)?;
    state.try_change_frame_mut(Frame::Equatorial)?;

    let mass_list = {
        if include_extended {
            MASSES
        } else {
            PLANETS
        }
    };

    let metadata = AccelSPKMeta {
        close_approach: None,
        non_grav_model,
        massive_obj: mass_list,
    };

    let (pos, vel, _meta) = {
        RadauIntegrator::integrate(
            &spk_accel,
            state.pos.into(),
            state.vel.into(),
            state.jd,
            jd_final,
            metadata,
        )?
    };

    let mut new_state = State::new(
        state.desig.to_owned(),
        jd_final,
        pos,
        vel,
        Frame::Equatorial,
        0,
    );
    spk.try_change_center(&mut new_state, center)?;
    new_state.try_change_frame_mut(frame)?;
    Ok(new_state)
}

/// Propagate an object using two body mechanics.
/// This is a brute force way to solve the kepler equations of motion as it uses Radau
/// as an integrator.
///
/// It is *strongly recommended* to use the `kepler.rs` code for this, as
/// it will be much more computationally efficient.
pub fn propagation_central(state: &State, jd_final: f64) -> NeosResult<[[f64; 3]; 2]> {
    let pos: Vector3<f64> = state.pos.into();
    let vel: Vector3<f64> = state.vel.into();
    let (pos, vel, _meta) = RadauIntegrator::integrate(
        &central_accel,
        pos,
        vel,
        state.jd,
        jd_final,
        CentralAccelMeta::default(),
    )?;
    Ok([pos.into(), vel.into()])
}

/// Propagate using n-body mechanics but skipping SPK queries.
/// This will propagate all planets and the Moon, so it may vary from SPK states slightly.
pub fn propagate_n_body_vec(
    states: Vec<State>,
    jd_final: f64,
    planet_states: Option<Vec<State>>,
    non_gravs: Vec<Option<NonGravModel>>,
) -> NeosResult<(Vec<State>, Vec<State>)> {
    if states.is_empty() {
        Err(Error::ValueError(
            "State vector is empty, propagation cannot continue".into(),
        ))?;
    }

    if non_gravs.len() != states.len() {
        Err(Error::ValueError(
            "Number of non-grav models doesnt match the number of provided objects.".into(),
        ))?;
    }

    let jd_init = states.first().unwrap().jd;

    let mut pos: Vec<f64> = Vec::new();
    let mut vel: Vec<f64> = Vec::new();
    let mut desigs: Vec<Desig> = Vec::new();

    let planet_states = planet_states.unwrap_or_else(|| {
        let spk = get_spk_singleton().try_read().unwrap();
        let mut planet_states = Vec::with_capacity(SIMPLE_PLANETS.len());
        for obj in SIMPLE_PLANETS.iter() {
            let planet = spk
                .try_get_state(obj.naif_id, jd_init, 10, Frame::Equatorial)
                .expect("Failed to find state for the provided initial jd");
            planet_states.push(planet);
        }
        planet_states
    });

    if planet_states.len() != SIMPLE_PLANETS.len() {
        Err(Error::ValueError(
            "Input planet states must contain the correct number of states.".into(),
        ))?;
    }
    if planet_states.first().unwrap().jd != jd_init {
        Err(Error::ValueError(
            "Planet states JD must match JD of input state.".into(),
        ))?;
    }
    for planet_state in planet_states.into_iter() {
        pos.append(&mut planet_state.pos.into());
        vel.append(&mut planet_state.vel.into());
        desigs.push(planet_state.desig);
    }

    for mut state in states.into_iter() {
        state.try_change_frame_mut(Frame::Equatorial)?;
        if jd_init != state.jd {
            Err(Error::ValueError(
                "All input states must have the same JD".into(),
            ))?;
        }
        if state.center_id != 10 {
            Err(Error::ValueError(
                "Center of all states must be 10 (the Sun).".into(),
            ))?
        }
        pos.append(&mut state.pos.into());
        vel.append(&mut state.vel.into());
        desigs.push(state.desig);
    }

    let meta = AccelVecMeta {
        non_gravs,
        massive_obj: SIMPLE_PLANETS,
    };

    let (pos, vel, _) = {
        RadauIntegrator::integrate(
            &vec_accel,
            DVector::from(pos),
            DVector::from(vel),
            jd_init,
            jd_final,
            meta,
        )?
    };
    let sun_pos = pos.fixed_rows::<3>(0);
    let sun_vel = vel.fixed_rows::<3>(0);
    let mut all_states: Vec<State> = Vec::new();
    for (idx, desig) in desigs.into_iter().enumerate() {
        let pos = pos.fixed_rows::<3>(idx * 3) - sun_pos;
        let vel = vel.fixed_rows::<3>(idx * 3) - sun_vel;
        let state = State::new(desig, jd_final, pos, vel, Frame::Equatorial, 10);
        all_states.push(state);
    }
    let final_states = all_states.split_off(SIMPLE_PLANETS.len());
    Ok((final_states, all_states))
}
