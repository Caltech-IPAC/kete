//! Functions related to two-body motion.
//! These are very fast to compute, however are not very accurate in multi-body systems
//! such as the solar system.
//!

use crate::constants::*;
use crate::errors::NEOSpyError;
use crate::fitting::newton_raphson;
use crate::state::State;
use nalgebra::{ComplexField, Vector3};
use std::f64::consts::TAU;

/// Compute the eccentric anomaly for all orbital classes.
///
/// # Arguments
///
/// * `ecc` - The eccentricity, must be non-negative.
/// * `mean_anomaly` - Mean anomaly, between 0 and 2 pi.
/// * `peri_dist` - Perihelion Distance in AU, only used for parabolic orbits.
///
pub fn compute_eccentric_anomaly(
    ecc: f64,
    mean_anom: f64,
    peri_dist: f64,
) -> Result<f64, NEOSpyError> {
    match ecc {
        ecc if ecc < 0.0 => Err(NEOSpyError::ValueError(
            "Eccentricity must be greater than 0".into(),
        )),
        ecc if ecc < 1e-6 => Ok(mean_anom),
        ecc if ecc < 0.9999 => {
            // Elliptical
            let f = |ecc_anom: f64| -ecc * ecc_anom.sin() + ecc_anom - mean_anom;
            let d = |ecc_anom: f64| 1.0 - 1.0 * ecc * ecc_anom.cos();
            Ok(newton_raphson(f, d, mean_anom % TAU, 1e-11)? % TAU)
        }
        ecc if ecc < 1.0001 => {
            // Parabolic
            let f = |ecc_anom: f64| -mean_anom + peri_dist * ecc_anom + ecc_anom.powi(3) / 6.0;
            let d = |ecc_anom: f64| peri_dist + ecc_anom.powi(2) / 2.0;
            Ok(newton_raphson(f, d, mean_anom.rem_euclid(TAU), 1e-11)? % TAU)
        }
        ecc => {
            // Hyperbolic
            let f = |ecc_anom: f64| ecc * ecc_anom.sinh() - ecc_anom - mean_anom;
            let d = |ecc_anom: f64| -1.0 + ecc * ecc_anom.cosh();
            Ok(newton_raphson(f, d, mean_anom.rem_euclid(TAU), 1e-11)?)
        }
    }
}

/// This function is used by the kepler universal solver below.
/// They are defined by the referenced paper.
///
///  Wisdom, Jack, and David M. Hernandez.
///  "A fast and accurate universal Kepler solver without Stumpff series."
///  Monthly Notices of the Royal Astronomical Society 453.3 (2015): 3015-3023.
///  https://arxiv.org/abs/1508.02699
fn g_1(s: f64, beta: f64) -> f64 {
    // the limit of this equation as beta approaches 0 is s
    if beta.abs() < 1e-14 {
        return s;
    }
    let beta_sqrt = beta.abs().sqrt();
    if beta >= 0.0 {
        (beta_sqrt * s).sin() / beta_sqrt
    } else {
        (beta_sqrt * s).sinh() / beta_sqrt
    }
}

/// This function is used by the kepler universal solver below.
/// They are defined by the referenced paper.
///
///  Wisdom, Jack, and David M. Hernandez.
///  "A fast and accurate universal Kepler solver without Stumpff series."
///  Monthly Notices of the Royal Astronomical Society 453.3 (2015): 3015-3023.
///  https://arxiv.org/abs/1508.02699
fn g_2(s: f64, beta: f64) -> f64 {
    // the limit of this equation as beta approaches 0 is s^2 / 2
    if beta.abs() < 1e-14 {
        return s.powi(2) / 2.0;
    }
    let beta_sqrt = beta.abs().sqrt() * s;
    if beta >= 0.0 {
        2.0 * (beta_sqrt / 2.0).sin().powf(2.0) / beta
    } else {
        (1.0 - (beta_sqrt).cosh()) / beta
    }
}

/// Solve the kepler equation for a universal formulation.
///
/// This roughly follows the outline of the algorithm described here:
///
///  Wisdom, Jack, and David M. Hernandez.
///  "A fast and accurate universal Kepler solver without Stumpff series."
///  Monthly Notices of the Royal Astronomical Society 453.3 (2015): 3015-3023.
///  https://arxiv.org/abs/1508.02699
///
/// # Arguments
///
/// * `dt` - The step size in days.
/// * `r0` - Distance from the central body (AU).
/// * `v0` - Velocity with respect to the central body (AU/day)
/// * `rv0` - R vector dotted with the V vector, not normalized.
///
fn solve_kepler_universal(
    mut dt: f64,
    r0: f64,
    v0: f64,
    rv0: f64,
) -> Result<(f64, f64), NEOSpyError> {
    // beta is GMS / semi_major
    let beta = 2.0 * GMS / r0 - v0.powi(2);
    let b_sqrt = beta.abs().sqrt();

    let res = {
        if beta.abs() < 1e-10 {
            // This is for parabolic orbits.
            // solve a cubic of the form:
            // x^3 + a2 * x^2 + a1 * x + a0 = 0
            let a0 = -6.0 * dt / GMS;
            let a1 = 6.0 * r0 / GMS;
            let a2 = 3.0 * rv0 / GMS;
            let p = (3.0 * a1 - a2.powi(2)) / 3.0;
            let q = (9.0 * a1 * a2 - 27.0 * a0 - 2.0 * a2.powi(3)) / 27.0;
            let w = (q / 2.0 + (q.powi(2) / 4.0 + p.powi(3) / 27.0).sqrt()).powf(1.0 / 3.0);
            let ans = w - p / (3.0 * w);
            Ok(ans)
        } else if beta > 0.0 {
            let period = GMS * b_sqrt.powi(-3);
            dt %= period * TAU;
            // elliptical orbits
            let f = |x: f64| {
                let g1 = (b_sqrt * x).sin() / b_sqrt;
                let g2 = 2.0 * (b_sqrt * x / 2.0).sin().powf(2.0) / beta;
                let g3 = (x - g1) / beta;
                r0 * g1 + rv0 * g2 + GMS * g3 - dt
            };
            let d = |x: f64| {
                let (g2d, g1d) = (b_sqrt * x).sin_cos();
                let g3d = (1.0 - g1d) / beta;
                r0 * g1d + rv0 * g2d + GMS * g3d
            };
            // Significantly better initial guess.
            let guess = b_sqrt.recip() * dt / period;
            newton_raphson(f, d, guess, 1e-11)
        } else {
            // hyperbolic orbits
            let f = |x: f64| {
                let g1 = (b_sqrt * x).sinh() / b_sqrt;
                let g2 = (1.0 - (b_sqrt * x).cosh()) / beta;
                let g3 = (x - g1) / beta;
                r0 * g1 + rv0 * g2 + GMS * g3 - dt
            };
            let d = |x: f64| {
                let (g2d, g1d) = (b_sqrt * x).sinh_cosh();
                let g3d = (1.0 - g1d) / beta;
                r0 * g1d + rv0 * g2d + GMS * g3d
            };
            newton_raphson(f, d, GMS_SQRT * beta.abs() * dt, 1e-11)
        }
    }?;
    Ok((res, beta))
}

/// Propagate an object forward in time by the specified amount assuming only 2 body
/// mechanics.
///
/// This will recursively call itself if it fails to converge, attempting to compute for
/// half intervals of time.
///
/// This roughly follows the outline of the algorithm described here:
///
///  Wisdom, Jack, and David M. Hernandez.
///  "A fast and accurate universal Kepler solver without Stumpff series."
///  Monthly Notices of the Royal Astronomical Society 453.3 (2015): 3015-3023.
///  <https://arxiv.org/abs/1508.02699>
///
/// # Arguments
///
/// * `time` - Time to propagate in days.
/// * `pos` - Starting position, from the center of the sun, in AU.
/// * `vel` - Starting velocity, in AU/Day.
/// * `depth` - Current Recursion depth, this should be called with `None` initially.
pub fn analytic_2_body(
    time: f64,
    pos: &Vector3<f64>,
    vel: &Vector3<f64>,
    depth: Option<usize>,
) -> Result<(Vector3<f64>, Vector3<f64>), NEOSpyError> {
    let mut depth = depth.to_owned().unwrap_or(0);
    if depth >= 10 {
        return Err(NEOSpyError::Convergence(
            "Two body recursion depth reached.".into(),
        ));
    } else {
        depth += 1;
    }
    if time.abs() < 1e-10 {
        return Ok((*pos, *vel));
    }
    let r0 = pos.norm();
    let v0 = vel.norm();
    let rv0 = pos.dot(vel);

    if rv0.is_nan() || rv0.is_infinite() {
        return Err(NEOSpyError::Convergence(
            "Input included infinity or NAN.".into(),
        ));
    }

    match solve_kepler_universal(time, r0, v0, rv0) {
        Ok((universal_s, beta)) => {
            let g1 = g_1(universal_s, beta);
            let g2 = g_2(universal_s, beta);
            let f = 1.0 - GMS / r0 * g2;
            let g = r0 * g1 + rv0 * g2;
            let new_pos = pos * f + vel * g;
            let new_r0 = new_pos.norm();
            let f_dot = -GMS / (new_r0 * r0) * g1;
            let g_dot = 1.0 - (GMS / new_r0) * g2;
            let new_vel = pos * f_dot + vel * g_dot;
            Ok((new_pos, new_vel))
        }
        Err(_) => {
            let (inter_pos, inter_vel) = analytic_2_body(0.5 * time, pos, vel, Some(depth))?;
            analytic_2_body(time * 0.5, &inter_pos, &inter_vel, Some(depth))
        }
    }
}

/// Propagate a [`State`] object using the analytic two body equations of motion.
///
/// Orbits are calculated where the central mass is located at (0, 0, 0) and has the
/// mass of the Sun. It is important that the appropriate center ID is chosen, it is
/// recommended to use a center of the Sun (10) for this computation.
///
/// This is the fastest method for getting a relatively good estimate of the orbits.
pub fn propagate_two_body(state: &State, jd_final: f64) -> Result<State, NEOSpyError> {
    let (pos, vel) = analytic_2_body(
        jd_final - state.jd,
        &state.pos.into(),
        &state.vel.into(),
        None,
    )?;

    Ok(State::new(
        state.desig.to_owned(),
        jd_final,
        pos,
        vel,
        state.frame,
        state.center_id,
    ))
}

#[cfg(test)]
mod tests {
    use std::f64::consts::TAU;

    use super::*;
    use nalgebra::Vector3;

    use super::compute_eccentric_anomaly;

    #[test]
    fn test_kepler_circular() {
        for r in [0.2, 0.5, 1.0, 2.0f64].iter() {
            let pos = Vector3::new(0.0, *r, 0.0);
            let vel = Vector3::new(-GMS_SQRT / r.sqrt(), 0.0, 0.0);
            let year = TAU / GMS_SQRT * r.powf(3.0 / 2.0);
            let res = analytic_2_body(year, &pos, &vel, None).unwrap();
            assert!((res.0 - pos).norm() < 1e-8);
            assert!((res.1 - vel).norm() < 1e-8);
        }
    }

    #[test]
    fn test_compute_eccentric_anom_hyperbolic() {
        let _ = compute_eccentric_anomaly(2.0, 63.21151553950512, 0.1).unwrap();
    }

    #[test]
    fn test_kepler_parabolic() {
        let pos = Vector3::new(0.0, 2.0, 0.0);
        let vel = Vector3::new(-GMS_SQRT, 0.0, 0.0);
        let year = TAU / GMS_SQRT;
        let res = analytic_2_body(year, &pos, &vel, None).unwrap();
        let pos_exp = Vector3::new(-4.448805955479905, -0.4739843046525608, 0.0);
        let vel_exp = Vector3::new(-0.00768983428326951, -0.008552645144187791, 0.0);
        assert!((res.0 - pos_exp).norm() < 1e-8);
        assert!((res.1 - vel_exp).norm() < 1e-8);
    }

    #[test]
    fn test_kepler_hyperbolic() {
        let pos = Vector3::new(0.0, 3.0, 0.0);
        let vel = Vector3::new(-GMS_SQRT, 0.0, 0.0);
        let year = TAU / GMS_SQRT;
        let res = analytic_2_body(year, &pos, &vel, None).unwrap();
        let pos_exp = Vector3::new(-5.556785268950049, 1.6076633958058089, 0.0);
        let vel_exp = Vector3::new(-0.013061655543084886, -0.005508140023183166, 0.0);
        assert!((res.0 - pos_exp).norm() < 1e-8);
        assert!((res.1 - vel_exp).norm() < 1e-8);
    }
}