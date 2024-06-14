//! Compute the acceleration of test particles with regard to massive objects.
//! This is used by the propagation code along with the numerical integrator to
//! calculated orbital dynamics.
//!
//! There are several functions defined here, which enable various levels of accuracy.
//!
//! These functions have a strict function signature, which is defined inside of the
//! radau integrator class. This function signature contains 4 terms:
//!
//! (time, x, x_der, &mut MetaData, exact_eval) -> Result<x_der_der, Error>
//!
//! Where `x` and its derivative `x_der` are vectors. This also accepts a mutable
//! reference to a metadata collection. Metadata may include things like object
//! specific orbit parameters such as the non-grav A terms, or keep track of close
//!
//!
//! `exact_eval` is a bool which is passed when the integrator is passing values where
//! the `x` and `x_der` are being evaluated at true locations. IE: where the integrator
//! thinks that the object should actually be. These times are when close encounter
//! information should be recorded.
//!
use crate::spice::get_spk_singleton;
use crate::{constants::*, errors::NEOSpyError, frames::Frame};
use itertools::Itertools;
use nalgebra::allocator::Allocator;
use nalgebra::{DefaultAllocator, Dim, Matrix, Matrix3, OMatrix, OVector, Vector3, U1, U2};
use std::ops::AddAssign;

/// Metadata object used by the [`central_accel`] function below.
#[derive(Debug)]
pub struct CentralAccelMeta {
    /// A vector of times where the central accel function was evaluated at.
    pub times: Vec<f64>,

    /// The position where the central accel function was evaluated.
    pub pos: Vec<Vector3<f64>>,

    /// The velocity where the central accel function was evaluated.
    pub vel: Vec<Vector3<f64>>,

    /// Scaling factor for central mass.
    pub mass_scaling: f64,
}

impl Default for CentralAccelMeta {
    fn default() -> Self {
        CentralAccelMeta {
            times: Vec::new(),
            pos: Vec::new(),
            vel: Vec::new(),
            mass_scaling: 1.0,
        }
    }
}

/// Compute the accel on an object which experiences acceleration due to the Sun only.
/// Integrating this with Radau should result in the same values as two-body analytic
/// integration.
///
/// # Arguments
///
/// * `time` - Time of the evaluation. This is saved in the metadata but otherwise
///            unused.
/// * `pos` - A vector which defines the position with respect to the Sun in AU.
/// * `vel` - A vector which defines the velocity with respect to the Sun in AU/Day.
/// * `meta` - Metadata object which records values at integration steps.
pub fn central_accel(
    time: f64,
    pos: &Vector3<f64>,
    vel: &Vector3<f64>,
    meta: &mut CentralAccelMeta,
    exact_eval: bool,
) -> Result<Vector3<f64>, NEOSpyError> {
    if exact_eval {
        meta.times.push(time);
        meta.pos.push(*pos);
        meta.vel.push(*vel);
    }

    Ok(-pos * pos.norm().powi(-3) * GMS)
}

/// Non-gravitational forces are modeled as defined on page 139 of the Comets II
/// textbook.
///
/// This model adds 3 "A" terms to the acceleration which the object feels. These
/// A terms represent additional radial, tangential, and normal forces on the object.
///
/// accel_additional = A_1 * g(r) * r_vec + A_2 * g(r) * t_vec + A_3 * g(r) * n_vec
/// Where r_vec, t_vec, n_vec are the radial, tangential, and normal unit vectors for
/// the object.
///
/// The g(r) function is defined by the equation:
/// g(r) = alpha (r / r0) ^ -m * (1 + (r / r0) ^ n) ^ -k
///
/// When alpha=1.0, n=0.0, k=0.0, r0=1.0, and m=2.0, this is equivalent to a 1/r^2
/// correction.
#[derive(Debug)]
pub struct NonGravModel {
    /// Constant for the radial non-gravitational force.
    pub a1: f64,

    /// Constant for the tangential non-gravitational force.
    pub a2: f64,

    /// Constant for the normal non-gravitational force.
    pub a3: f64,

    /// Coefficients for the g(r) function defined above.
    pub alpha: f64,
    /// Coefficients for the g(r) function defined above.
    pub r_0: f64,
    /// Coefficients for the g(r) function defined above.
    pub m: f64,
    /// Coefficients for the g(r) function defined above.
    pub n: f64,
    /// Coefficients for the g(r) function defined above.
    pub k: f64,
}

impl NonGravModel {
    /// Construct a new non-grav model, manually specifying all parameters.
    /// Consider using the other constructors if this is a simple object.
    #[allow(clippy::too_many_arguments)]
    pub fn new(a1: f64, a2: f64, a3: f64, alpha: f64, r_0: f64, m: f64, n: f64, k: f64) -> Self {
        Self {
            a1,
            a2,
            a3,
            alpha,
            r_0,
            m,
            n,
            k,
        }
    }

    /// Construct a new non-grav model which follows the 1/r^2 drop-off.
    pub fn new_r2(a1: f64, a2: f64, a3: f64) -> Self {
        Self {
            a1,
            a2,
            a3,
            alpha: 1.0,
            r_0: 1.0,
            m: 2.0,
            n: 0.0,
            k: 0.0,
        }
    }

    /// Construct a new non-grav model which follows the default comet drop-off.
    pub fn new_comet_default(a1: f64, a2: f64, a3: f64) -> Self {
        Self {
            a1,
            a2,
            a3,
            alpha: 0.111262,
            r_0: 2.808,
            m: 2.15,
            n: 5.093,
            k: 4.6142,
        }
    }

    /// Compute the non-gravitational acceleration vector when provided the position
    /// and velocity vector with respect to the sun.
    pub fn accel_vec(&self, pos: &Vector3<f64>, vel: &Vector3<f64>) -> Vector3<f64> {
        let pos_norm = pos.normalize();
        let rr0 = pos.norm() / self.r_0;
        let scale = self.alpha * rr0.powf(-self.m) * (1.0 + rr0.powf(self.n)).powf(-self.k);
        let t_vec = (vel - pos_norm * vel.dot(&pos_norm)).normalize();
        let n_vec = t_vec.cross(&pos_norm).normalize();
        let mut accel = pos_norm * (scale * self.a1);
        accel += t_vec * (scale * self.a2);
        accel += n_vec * (scale * self.a3);
        accel
    }
}

/// Metadata for the [`spk_accel`] function defined below.
#[derive(Debug)]
pub struct AccelSPKMeta<'a> {
    /// Closest approach to a massive object.
    /// This records the ID of the object, time, and distance in AU.
    pub close_approach: Option<(i32, f64, f64)>,

    /// `A` terms of the non-gravitational forces.
    /// If this is not provided, only standard gravitational model is applied.
    /// If these values are provided, then the effects of the A terms are added in
    /// addition to standard forces.
    pub non_grav_a: Option<NonGravModel>,

    /// The list of massive objects to apply during SPK computation.
    /// This list contains the ID of the object in the SPK along with the mass and
    /// radius of the object. Mass is given in fractions of solar mass and radius is
    /// in AU.
    pub massive_obj: &'a [(i32, f64, f64)],
}

/// Compute the accel on an object which experiences acceleration due to all massive
/// objects contained within the Spice Kernel SPKs. This uses the planets and the Moon,
/// and applies a General Relativity correction for the Sun and Jupiter.
///
/// Whatever objects are present in the metadata struct will be used as sources of mass.
/// These objects' states are queried from the loaded SPICE kernels, so this will fail
/// if the SPICE kernels are not present/loaded.
///
/// Typically this relies on DE440s, which contains the years ~1800-2200.
///
/// Metadata:
///     Metadata here records the closest approach to any planet. This value is updated
///     every time the function is called.
///
/// # Arguments
///
/// * `time` - Time of the evaluation in JD in TDB scaled time multiplied by SUN_GMS_SQRT.
/// * `pos` - A vector which defines the position with respect to the Sun in AU.
/// * `vel` - A vector which defines the velocity with respect to the Sun in AU/Day multiplied by SUN_GMS_SQRT.
/// * `meta` - Metadata object [`AccelSPKMeta`] which records values at each integration step.
pub fn spk_accel(
    time: f64,
    pos: &Vector3<f64>,
    vel: &Vector3<f64>,
    meta: &mut AccelSPKMeta,
    exact_eval: bool,
) -> Result<Vector3<f64>, NEOSpyError> {
    let mut accel = Vector3::<f64>::zeros();

    if exact_eval {
        if let Some(close_approach) = meta.close_approach.as_mut() {
            if close_approach.2 == 0.0 {
                *close_approach = (-1, 1000000.0, 0.0)
            }
        }
    }

    let spk = get_spk_singleton().try_read().unwrap();

    for (id, mass, radius) in meta.massive_obj.iter() {
        let state = spk.try_get_state(*id, time, 0, Frame::Equatorial)?;
        let rel_pos: Vector3<f64> = pos - Vector3::from(state.pos);
        let rel_pos_norm = rel_pos.normalize();
        let r = rel_pos.norm();

        if exact_eval {
            if let Some(close_approach) = meta.close_approach.as_mut() {
                if close_approach.2 == r {
                    *close_approach = (*id, r, time);
                }
            }
        }
        if r <= *radius {
            return Err(NEOSpyError::Impact(*id, time));
        }

        let r2_inv = r.powi(-2);
        accel -= rel_pos_norm * r2_inv * *mass * GMS;

        // if it is the sun or jupiter, apply a correction for GR
        if *id == 10 || *id == 5 {
            let r3_inv = r.powi(-3);
            let rel_vel: Vector3<f64> = vel - Vector3::from(state.vel);

            let r_v = 4.0 * rel_pos.dot(&rel_vel);

            let rel_v2: f64 = rel_vel.norm_squared();
            let gr_const: f64 = mass * C_AU_PER_DAY_INV_SQUARED * r3_inv * GMS;
            let c: f64 = 4. * mass * GMS / r - rel_v2;

            accel += gr_const * (c * rel_pos + r_v * rel_vel);

            if *id == 10 {
                // J2 for the Sun
                let coef = SUN_J2 * GMS * *mass * r.powi(-5) * 1.5 * radius.powi(2);
                let z2 = 5.0 * rel_pos_norm.z.powi(2);
                accel[0] -= rel_pos.x * coef * (z2 - 1.0);
                accel[1] -= rel_pos.y * coef * (z2 - 1.0);
                accel[2] -= rel_pos.z * coef * (z2 - 3.0);

                // non-gravitational forces
                if let Some(model) = &meta.non_grav_a {
                    accel += model.accel_vec(&rel_pos, &rel_vel)
                }
            }
        }
    }
    Ok(accel)
}

/// Metadata for the [`vec_accel`] function defined below.
#[derive(Debug)]
pub struct AccelVecMeta<'a, D: Dim>
where
    DefaultAllocator: Allocator<f64, D, U2>,
{
    /// `A` terms of the non-gravitational forces.
    /// If this is not provided, only standard gravitational model is applied.
    /// If these values are provided, then the effects of the A terms are added in
    /// addition to standard forces.
    pub non_grav_a: Option<OMatrix<f64, D, U2>>,

    /// The list of massive objects to apply during SPK computation.
    /// This list contains the ID of the object in the SPK along with the mass and
    /// radius of the object. Mass is given in fractions of solar mass and radius is
    /// in AU.
    pub massive_obj: &'a [(i32, f64, f64)],
}

/// Compute the accel on an object which experiences acceleration due to all massive
/// objects. This assumes that the first N objects match the objects in the metadata
/// list in order. IE: if MASSIVE_OBJECTS from the constants file is used in the meta
/// data, then those objects are assumed to be in the same order in the pos/vel vectors
/// provided.
///
/// # Arguments
///
/// * `time` - Time is not used in this.
/// * `pos` - A vector which defines the position with respect to the Sun in AU.
/// * `vel` - A vector which defines the velocity with respect to the Sun in AU/Day.
/// * `meta` - Metadata.
pub fn vec_accel<D: Dim>(
    time: f64,
    pos: &OVector<f64, D>,
    vel: &OVector<f64, D>,
    meta: &mut AccelVecMeta<D>,
    _exact_eval: bool,
) -> Result<OVector<f64, D>, NEOSpyError>
where
    DefaultAllocator: Allocator<f64, D> + Allocator<f64, D, U2>,
{
    // objects in the pos/vel vectors are setup like so
    // (x, y, z, x, y, z, ...

    let n_objects = pos.len() / 3;

    let (dim, _) = pos.shape_generic();
    let mut accel = Matrix::zeros_generic(dim, U1);

    for idx in 0..n_objects {
        let pos_idx = pos.rows(idx * 3, 3);

        for (idy, (id, mass, radius)) in meta.massive_obj.iter().enumerate() {
            if idx == idy {
                continue;
            }
            let pos_idy = pos.rows(idy * 3, 3);
            let rel_pos = pos_idx - pos_idy;
            let rel_pos_norm = rel_pos.normalize();
            let r = rel_pos.norm();

            if r <= *radius {
                return Err(NEOSpyError::Impact(*id, time));
            }

            let r2_inv = r.powi(-2);
            accel
                .rows_mut(idx * 3, 3)
                .add_assign(-&rel_pos_norm * r2_inv * *mass * GMS);

            // if it is the sun or jupiter, apply a correction for GR
            if *id == 10 || *id == 5 {
                let vel_idx = vel.rows(idx * 3, 3);
                let vel_idy = vel.rows(idy * 3, 3);
                let r3_inv = r.powi(-3);
                let rel_vel = vel_idx - vel_idy;

                let r_v = 4.0 * rel_pos.dot(&rel_vel);

                let rel_v2: f64 = rel_vel.norm_squared();
                let gr_const: f64 = mass * C_AU_PER_DAY_INV_SQUARED * r3_inv * GMS;
                let c: f64 = 4. * mass * GMS / r - rel_v2;

                accel
                    .rows_mut(idx * 3, 3)
                    .add_assign(gr_const * (c * rel_pos + r_v * &rel_vel));

                // Add non-grav forces if defined.
                if *id == 10 {
                    if let Some(a_matrix) = meta.non_grav_a.as_ref() {
                        let a_row = a_matrix.row(idx);
                        let (a1, a2) = a_row.iter().collect_tuple().unwrap();
                        accel
                            .rows_mut(idx * 3, 3)
                            .add_assign(&rel_pos_norm * r2_inv * *mass * GMS * *a1);
                        accel
                            .rows_mut(idx * 3, 3)
                            .add_assign(rel_vel.normalize() * *mass * GMS * *a2);
                    }
                }
            }
        }
    }
    Ok(accel)
}

/// Calculate the Jacobian for the central_accel function.
///
/// This enables the computation of the STM.
pub fn central_accel_grad(
    _time: f64,
    pos: &Vector3<f64>,
    _vel: &Vector3<f64>,
    meta: &mut CentralAccelMeta,
) -> Matrix3<f64> {
    let zeros = Vector3::<f64>::zeros();
    accel_grad(pos, _vel, &zeros, &zeros, meta.mass_scaling)
}

/// Calculate the Jacobian for the central_accel function.
///
/// This enables the computation of the STM.
pub fn accel_grad(
    obj_pos: &Vector3<f64>,
    _obj_vel: &Vector3<f64>,
    mass_pos: &Vector3<f64>,
    _mass_vel: &Vector3<f64>,
    mass: f64,
) -> Matrix3<f64> {
    let pos = obj_pos - mass_pos;
    let r = pos.norm();
    let r_2 = r.powi(2);
    let r_5_inv = r.powi(5);
    Matrix3::<f64>::new(
        r_2 - 3.0 * pos.x.powi(2),
        -3.0 * pos.x * pos.y,
        -3.0 * pos.x * pos.z,
        -3.0 * pos.x * pos.y,
        r_2 - 3.0 * pos.y.powi(2),
        -3.0 * pos.y * pos.z,
        -3.0 * pos.x * pos.z,
        -3.0 * pos.y * pos.z,
        r_2 - 3.0 * pos.z.powi(2),
    ) / (2.0 * r_5_inv)
        * mass
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use super::*;

    #[test]
    fn check_accelerations_equal() {
        let spk = get_spk_singleton().try_read().unwrap();
        let jd = 2451545.0;
        let mut pos: Vec<f64> = Vec::new();
        let mut vel: Vec<f64> = Vec::new();

        for (id, _mass, _radius) in MASSIVE_OBJECTS.iter() {
            let planet = spk.try_get_state(*id, jd, 0, Frame::Equatorial).unwrap();
            pos.append(&mut planet.pos.into());
            vel.append(&mut planet.vel.into());
        }

        pos.append(&mut [0.0, 0.0, 0.5].into());
        vel.append(&mut [0.0, 0.0, 1.0].into());

        let accel = vec_accel(
            jd,
            &pos.into(),
            &vel.into(),
            &mut AccelVecMeta {
                non_grav_a: None,
                massive_obj: MASSIVE_OBJECTS,
            },
            false,
        )
        .unwrap()
        .iter()
        .copied()
        .skip(MASSIVE_OBJECTS.len() * 3)
        .collect_vec();

        let accel2 = spk_accel(
            jd,
            &[0.0, 0.0, 0.5].into(),
            &[0.0, 0.0, 1.0].into(),
            &mut AccelSPKMeta {
                close_approach: None,
                non_grav_a: None,
                massive_obj: MASSIVE_OBJECTS,
            },
            false,
        )
        .unwrap();
        assert!((accel[0] - accel2[0]).abs() < 1e-10);
        assert!((accel[1] - accel2[1]).abs() < 1e-10);
        assert!((accel[2] - accel2[2]).abs() < 1e-10);
    }
}
