use nalgebra::Vector3;

use crate::constants::{C_AU_PER_DAY_INV_SQUARED, GMS};

/// Non-Gravitational models.
/// These are used during integration to model non-gravitational forces on particles in
/// the solar system.
#[derive(Debug, Clone)]
pub enum NonGravModel {
    /// JPL's non-gravitational forces are modeled as defined on page 139 of the Comets II
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
    JplComet {
        /// Constant for the radial non-gravitational force.
        a1: f64,
        /// Constant for the tangential non-gravitational force.
        a2: f64,
        /// Constant for the normal non-gravitational force.
        a3: f64,
        /// Coefficients for the g(r) function defined above.
        alpha: f64,
        /// Coefficients for the g(r) function defined above.
        r_0: f64,
        /// Coefficients for the g(r) function defined above.
        m: f64,
        /// Coefficients for the g(r) function defined above.
        n: f64,
        /// Coefficients for the g(r) function defined above.
        k: f64,
    },

    /// Dust model, including Solar Radiation Pressure (SRP) and the Poynting-Robertson
    /// effect.
    ///
    /// SRP acts as an effective reduction in the gravitational force of the
    /// Sun, reducing the central acceleration force in the radial direction.
    ///
    /// Poynting-Robertson acts as a drag force, in the opposite direction of motion.
    Dust {
        /// Beta Parameter
        beta: f64,
    },
}

impl NonGravModel {
    /// Construct a new non-grav model, manually specifying all parameters.
    /// Consider using the other constructors if this is a simple object.
    #[allow(clippy::too_many_arguments)]
    pub fn new_jpl(
        a1: f64,
        a2: f64,
        a3: f64,
        alpha: f64,
        r_0: f64,
        m: f64,
        n: f64,
        k: f64,
    ) -> Self {
        Self::JplComet {
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

    /// Construct a new non-grav dust model.
    pub fn new_dust_raw(beta: f64) -> Self {
        Self::Dust { beta }
    }

    /// Construct a new non-grav model which follows the default comet drop-off.
    pub fn new_jpl_comet_default(a1: f64, a2: f64, a3: f64) -> Self {
        Self::JplComet {
            a1,
            a2,
            a3,
            alpha: 0.1112620426,
            r_0: 2.808,
            m: 2.15,
            n: 5.093,
            k: 4.6142,
        }
    }

    /// Compute the non-gravitational acceleration vector when provided the position
    /// and velocity vector with respect to the sun.
    pub fn accel_vec(&self, pos: &Vector3<f64>, vel: &Vector3<f64>) -> Vector3<f64> {
        match self {
            Self::Dust { beta } => {
                let pos_norm = pos.normalize();
                let r_dot = &pos_norm.dot(vel);
                let norm2_inv = pos.norm_squared().recip();
                let scaling = GMS * beta * norm2_inv;
                scaling
                    * ((1.0 - r_dot * C_AU_PER_DAY_INV_SQUARED) * pos_norm
                        - vel * C_AU_PER_DAY_INV_SQUARED)
            }

            Self::JplComet {
                a1,
                a2,
                a3,
                alpha,
                r_0,
                m,
                n,
                k,
            } => {
                let pos_norm = pos.normalize();
                let rr0 = pos.norm() / r_0;
                let scale = alpha * rr0.powf(-m) * (1.0 + rr0.powf(*n)).powf(-k);
                let t_vec = (vel - pos_norm * vel.dot(&pos_norm)).normalize();
                let n_vec = t_vec.cross(&pos_norm).normalize();
                let mut accel = pos_norm * (scale * a1);
                accel += t_vec * (scale * a2);
                accel += n_vec * (scale * a3);
                accel
            }
        }
    }
}
