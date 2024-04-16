use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// Reflected light properties of a comet using the MK magnitude system.
///
/// <https://en.wikipedia.org/wiki/Absolute_magnitude#Cometary_magnitudes>
///
/// This model additionally includes a correction for phase effects.
#[derive(Debug, Deserialize, Serialize)]
pub struct CometMKParams {
    /// Designation (name) of the object.
    pub desig: String,

    /// M1 and K1 if defined.
    pub mk_1: Option<[f64; 2]>,

    /// M2 and K2 if defined.
    pub mk_2: Option<[f64; 2]>,

    /// Phase correction coefficient in units of Mag/Deg
    pub phase_corr_coef: f64,
}

impl CometMKParams {
    /// Create a new CometMKParams object.
    pub fn new(
        desig: String,
        mk_1: Option<[f64; 2]>,
        mk_2: Option<[f64; 2]>,
        phase_corr: f64,
    ) -> Self {
        Self {
            desig,
            mk_1,
            mk_2,
            phase_corr_coef: phase_corr,
        }
    }

    /// Compute the apparent total flux including both coma and nucleus of the comet.
    /// This includes an additional 0.035 Mag/Deg phase correction.
    pub fn apparent_total_flux(
        &self,
        sun2obs: &Vector3<f64>,
        sun2obj: &Vector3<f64>,
    ) -> Option<f64> {
        let [m1, k1] = self.mk_1?;
        let obj2obs = -sun2obj + sun2obs;
        let obs_dist = obj2obs.norm();
        let helio_dist = sun2obj.norm();
        let phase_corr = self.phase_corr_coef * obj2obs.angle(&-sun2obj).to_degrees();
        Some(m1 + k1 * helio_dist.log10() + 5.0 * obs_dist.log10() + phase_corr)
    }

    /// Compute the apparent nuclear flux of the comet, not including the coma.
    /// This includes an additional 0.035 Mag/Deg phase correction.
    pub fn apparent_nuclear_flux(
        &self,
        sun2obs: &Vector3<f64>,
        sun2obj: &Vector3<f64>,
    ) -> Option<f64> {
        let [m2, k2] = self.mk_2?;
        let obj2obs = -sun2obj + sun2obs;
        let obs_dist = obj2obs.norm();
        let helio_dist = sun2obj.norm();
        let phase_corr = self.phase_corr_coef * obj2obs.angle(&-sun2obj).to_degrees();
        Some(m2 + k2 * helio_dist.log10() + 5.0 * obs_dist.log10() + phase_corr)
    }
}
