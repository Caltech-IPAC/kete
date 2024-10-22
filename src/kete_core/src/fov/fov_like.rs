//! # Field of View like trait
//! This trait defines field of view checks for portions of the sky.

use super::*;
use nalgebra::Vector3;
use rayon::prelude::*;

use crate::constants::C_AU_PER_DAY_INV;
use crate::prelude::*;

/// Field of View like objects.
/// These may contain multiple unique sky patches, so as a result the expected
/// behavior is to return the index as well as the [`Contains`] for the closest
/// sky patch.
pub trait FovLike: Sync + Sized {
    /// Return the FOV of the patch at the specified index.
    /// This will panic if the index is out of allowed bounds.
    fn get_fov(&self, index: usize) -> FOV;

    /// Position of the observer.
    fn observer(&self) -> &State;

    /// Is the specified vector contained within this FOVLike object.
    /// A ['Contains'] is returned for each sky patch.
    fn contains(&self, obs_to_obj: &Vector3<f64>) -> (usize, Contains);

    /// Number of sky patches contained within this FOV.
    fn n_patches(&self) -> usize;

    /// Change the target frame to the new frame.
    fn try_frame_change_mut(&mut self, new_frame: Frame) -> KeteResult<()>;

    /// Check if a static source is visible. This assumes the vector passed in is at an
    /// infinite distance from the observer.
    #[inline]
    fn check_static(&self, pos: &Vector3<f64>) -> (usize, Contains) {
        self.contains(pos)
    }

    /// Assuming the object undergoes linear motion, check to see if it is within the
    /// field of view.
    #[inline]
    fn check_linear(&self, state: &State) -> (usize, Contains, State) {
        let pos: Vector3<_> = state.pos.into();
        let vel = state.vel.into();
        let obs = self.observer();

        let obs_pos: Vector3<_> = obs.pos.into();

        let rel_pos = pos - obs_pos;

        // This also accounts for first order light delay.
        let dt = obs.jd - state.jd - rel_pos.norm() * C_AU_PER_DAY_INV;
        let new_pos = pos + vel * dt;
        let new_rel_pos = new_pos - obs_pos;
        let (idx, contains) = self.contains(&new_rel_pos);
        let new_state = State::new(
            state.desig.clone(),
            obs.jd + dt,
            new_pos,
            vel,
            obs.frame,
            obs.center_id,
        );
        (idx, contains, new_state)
    }

    /// Assuming the object undergoes two-body motion, check to see if it is within the
    /// field of view.
    #[inline]
    fn check_two_body(&self, state: &State) -> KeteResult<(usize, Contains, State)> {
        let obs = self.observer();
        let obs_pos: Vector3<_> = obs.pos.into();

        // bring state up to observer time.
        let final_state = propagate_two_body(state, obs.jd)?;

        // correct for light delay
        let dt = -(Vector3::from(final_state.pos) - obs_pos).norm() * C_AU_PER_DAY_INV;
        let final_state = propagate_two_body(&final_state, obs.jd + dt)?;
        let rel_pos = Vector3::from(final_state.pos) - obs_pos;

        let (idx, contains) = self.contains(&rel_pos);
        Ok((idx, contains, final_state))
    }

    /// Assuming the object undergoes n-body motion, check to see if it is within the
    /// field of view.
    #[inline]
    fn check_n_body(
        &self,
        state: &State,
        include_asteroids: bool,
    ) -> KeteResult<(usize, Contains, State)> {
        let obs = self.observer();
        let obs_pos = Vector3::from(obs.pos);

        let exact_state = propagate_n_body_spk(state.clone(), obs.jd, include_asteroids, None)?;

        // correct for light delay
        let dt = -(Vector3::from(exact_state.pos) - obs_pos).norm() * C_AU_PER_DAY_INV;
        let final_state = propagate_two_body(&exact_state, obs.jd + dt)?;
        let rel_pos = Vector3::from(final_state.pos) - obs_pos;

        let (idx, contains) = self.contains(&rel_pos);

        Ok((idx, contains, final_state))
    }

    /// Given a list of states, check to see if the objects are visible at the desired time.
    ///
    /// Only the final observed states are returned, if the object was not seen it will not
    /// be returned.
    ///
    /// This does progressively more exact checks.
    /// - If the propagation time is under the specified `dt_limit`, it is assumed that the
    ///   objects obey two body propagation within a small error. In this case they are
    ///   first checked using the assumption that everything is linear, and are checked
    ///   using the rough field of view, this field of view should be larger than the
    ///   exact fov. The linear check assumes the object and observer are moving in straight
    ///   lines.
    /// - If the propagation time between the input state and the fov time is greater than
    ///   the `dt_limit` then two body propagation is used for the rough field of view
    ///   check. If this rough field of view check is passed, then use n-body propagation to
    ///   get exact position.
    ///
    /// # Arguments
    ///
    /// * `state` - A vector of States which define the objects, the center ID should be set
    ///             to 10 (the Sun).
    /// * `dt_limit` - Length of time in days where two body motion is considered valid.
    /// * `include_asteroids` - Include the 5 largest asteroids during the computation.
    ///
    fn check_visible(
        &self,
        states: &[State],
        dt_limit: f64,
        include_asteroids: bool,
    ) -> Vec<Option<SimultaneousStates>> {
        let obs_state = self.observer();

        let final_states: Vec<(usize, State)> = states
            .iter()
            .filter_map(|state: &State| {
                // assuming linear motion, how far can the object have moved relative
                // to the observer? Then add a factor of 2 for safety
                let max_dist = (Vector3::from(state.vel) - Vector3::from(obs_state.vel)).norm()
                    * dt_limit
                    * 2.0;

                if (state.jd - obs_state.jd).abs() < dt_limit {
                    let (_, contains, _) = self.check_linear(state);
                    if let Contains::Outside(dist) = contains {
                        if dist > max_dist {
                            return None;
                        }
                    }
                    let (idx, contains, state) = self.check_two_body(state).ok()?;
                    match contains {
                        Contains::Inside => Some((idx, state)),
                        _ => None,
                    }
                } else {
                    let (_, contains, _) = self.check_two_body(state).ok()?;
                    if let Contains::Outside(dist) = contains {
                        if dist > max_dist {
                            return None;
                        }
                    }
                    let (idx, contains, state) =
                        self.check_n_body(state, include_asteroids).ok()?;
                    match contains {
                        Contains::Inside => Some((idx, state)),
                        _ => None,
                    }
                }
            })
            .collect();

        let mut detector_states = vec![Vec::<State>::new(); self.n_patches()];
        for (idx, state) in final_states.into_iter() {
            detector_states[idx].push(state)
        }

        detector_states
            .into_iter()
            .enumerate()
            .map(|(idx, states)| {
                SimultaneousStates::new_exact(states, Some(self.get_fov(idx))).ok()
            })
            .collect()
    }

    /// Given an object ID, attempt to load the object from the SPKs.
    /// This will fail silently if the object is not found.
    fn check_spks(&self, obj_ids: &[i64]) -> Vec<Option<SimultaneousStates>> {
        let obs = self.observer();
        let spk = get_spk_singleton().try_read().unwrap();

        let mut visible: Vec<Vec<State>> = vec![Vec::new(); self.n_patches()];

        let states: Vec<_> = obj_ids
            .into_par_iter()
            .filter_map(|&obj_id| {
                match spk.try_get_state(obj_id, obs.jd, obs.center_id, obs.frame) {
                    Ok(state) => match self.check_two_body(&state) {
                        Ok((idx, Contains::Inside, state)) => Some((idx, state)),
                        _ => None,
                    },
                    _ => None,
                }
            })
            .collect();

        states
            .into_iter()
            .for_each(|(patch_idx, state)| visible[patch_idx].push(state));

        visible
            .into_iter()
            .enumerate()
            .map(|(idx, states_patch)| {
                SimultaneousStates::new_exact(states_patch, Some(self.get_fov(idx))).ok()
            })
            .collect()
    }

    /// Given a collection of static positions, return the index of the input vector
    /// which was visible.
    fn check_statics(&self, pos: &[Vector3<f64>]) -> Vec<Option<(Vec<usize>, FOV)>> {
        let mut visible: Vec<Vec<usize>> = vec![Vec::new(); self.n_patches()];

        pos.iter().enumerate().for_each(|(vec_idx, p)| {
            if let (patch_idx, Contains::Inside) = self.check_static(p) {
                visible[patch_idx].push(vec_idx)
            }
        });

        visible
            .into_iter()
            .enumerate()
            .map(|(idx, vis_patch)| {
                if vis_patch.is_empty() {
                    None
                } else {
                    Some((vis_patch, self.get_fov(idx)))
                }
            })
            .collect()
    }
}
