//! # Field of View
//! On-Sky field of view checks.
pub mod fov_like;
pub mod generic;
pub mod neos;
pub mod patches;
pub mod wise;
pub mod ztf;

pub use fov_like::*;
pub use generic::*;
pub use neos::*;
pub use patches::*;
pub use wise::*;
pub use ztf::*;

use serde::{Deserialize, Serialize};

use crate::prelude::*;

/// Allowed FOV objects, either contiguous or joint.
/// Many of these exist solely to carry additional metadata.
#[derive(Debug, Clone, Deserialize, Serialize)]
pub enum FOV {
    /// WISE or NEOWISE FOV.
    Wise(WiseCmos),

    /// NEOS single cmos FOV.
    NeosCmos(NeosCmos),

    /// ZTF Single Quad of single CCD FOV.
    ZtfCcdQuad(ZtfCcdQuad),

    /// Generic cone FOV without any additional metadata.
    GenericCone(GenericCone),

    /// Generic rectangle FOV without any additional metadata.
    GenericRectangle(GenericRectangle),

    /// Full ZTF field of up to 64 individual files.
    ZtfField(ZtfField),

    /// NEOS Visit.
    NeosVisit(NeosVisit),
}

impl FOV {
    /// Check if a collection of states are visible to this FOV using orbital propagation
    pub fn check_visible(self, states: &[State], dt_limit: f64) -> Vec<Option<SimultaneousStates>> {
        match self {
            FOV::Wise(fov) => fov.check_visible(states, dt_limit),
            FOV::NeosCmos(fov) => fov.check_visible(states, dt_limit),
            FOV::ZtfCcdQuad(fov) => fov.check_visible(states, dt_limit),
            FOV::GenericCone(fov) => fov.check_visible(states, dt_limit),
            FOV::GenericRectangle(fov) => fov.check_visible(states, dt_limit),
            FOV::ZtfField(fov) => fov.check_visible(states, dt_limit),
            FOV::NeosVisit(fov) => fov.check_visible(states, dt_limit),
        }
    }

    /// Observer position in this FOV
    pub fn observer(&self) -> &State {
        match self {
            FOV::Wise(fov) => fov.observer(),
            FOV::NeosCmos(fov) => fov.observer(),
            FOV::ZtfCcdQuad(fov) => fov.observer(),
            FOV::GenericCone(fov) => fov.observer(),
            FOV::GenericRectangle(fov) => fov.observer(),
            FOV::ZtfField(fov) => fov.observer(),
            FOV::NeosVisit(fov) => fov.observer(),
        }
    }

    /// Check if any loaded SPK objects are visible to this FOV
    pub fn check_spks(&self, obj_ids: &[isize]) -> Vec<Option<SimultaneousStates>> {
        match self {
            FOV::Wise(fov) => fov.check_spks(obj_ids),
            FOV::NeosCmos(fov) => fov.check_spks(obj_ids),
            FOV::ZtfCcdQuad(fov) => fov.check_spks(obj_ids),
            FOV::GenericCone(fov) => fov.check_spks(obj_ids),
            FOV::GenericRectangle(fov) => fov.check_spks(obj_ids),
            FOV::ZtfField(fov) => fov.check_spks(obj_ids),
            FOV::NeosVisit(fov) => fov.check_spks(obj_ids),
        }
    }
}
