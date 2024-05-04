//! # NEOS field of views
use super::{closest_inside, Contains, FovLike, OnSkyRectangle, SkyPatch, FOV};
use crate::constants::{NEOS_HEIGHT, NEOS_WIDTH};
use crate::prelude::*;
use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// NEOS frame data, a single detector on a single band
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct NeosCmos {
    /// State of the observer
    observer: State,

    /// Patch of sky
    pub patch: OnSkyRectangle,

    /// Rotation of the FOV.
    pub rotation: f64,

    /// Side ID
    pub side_id: u16,

    /// Stack ID
    pub stack_id: u8,

    /// Quad ID
    pub quad_id: u8,

    /// Loop ID
    pub loop_id: u8,

    /// Subloop ID
    pub subloop_id: u8,

    /// Exposure ID
    pub exposure_id: u8,

    /// Wavelength band, either 1 or 2 for NC1 or NC2
    pub band: u8,

    /// CMOS ID
    /// ID number of the CMOS chip, 0, 1, 2, or 3
    pub cmos_id: u8,
}

impl NeosCmos {
    /// Create a NEOS FOV
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        pointing: Vector3<f64>,
        rotation: f64,
        observer: State,
        side_id: u16,
        stack_id: u8,
        quad_id: u8,
        loop_id: u8,
        subloop_id: u8,
        exposure_id: u8,
        cmos_id: u8,
        band: u8,
    ) -> Self {
        let patch =
            OnSkyRectangle::new(pointing, rotation, NEOS_WIDTH, NEOS_HEIGHT, observer.frame);
        Self {
            observer,
            patch,
            side_id,
            stack_id,
            quad_id,
            loop_id,
            subloop_id,
            exposure_id,
            cmos_id,
            band,
            rotation,
        }
    }
}

impl FovLike for NeosCmos {
    fn get_fov(&self, index: usize) -> FOV {
        if index != 0 {
            panic!("FOV only has a single patch")
        }
        FOV::NeosCmos(self.clone())
    }

    #[inline]
    fn observer(&self) -> &State {
        &self.observer
    }

    #[inline]
    fn contains(&self, obs_to_obj: &Vector3<f64>) -> (usize, Contains) {
        (0, self.patch.contains(obs_to_obj))
    }

    #[inline]
    fn n_patches(&self) -> usize {
        1
    }
}

/// NEOS frame data, all 8 chips of a visit.
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct NeosVisit {
    /// Individual CMOS fields
    chips: Box<[NeosCmos; 8]>,

    /// Observer position
    observer: State,

    /// Rotation of the FOV.
    pub rotation: f64,

    /// Side ID
    pub side_id: u16,

    /// Stack ID
    pub stack_id: u8,

    /// Quad ID
    pub quad_id: u8,

    /// Loop ID
    pub loop_id: u8,

    /// Subloop ID
    pub subloop_id: u8,

    /// Exposure ID
    pub exposure_id: u8,
}

impl NeosVisit {
    /// Construct a new NeosVisit from a list of cmos fovs.
    /// These cmos fovs must be from the same metadata when appropriate.
    pub fn new(chips: Vec<NeosCmos>) -> Result<Self, NEOSpyError> {
        if chips.len() != 8 {
            return Err(NEOSpyError::ValueError(
                "Visit must contains 8 NeosCmos fovs".into(),
            ));
        }
        let chips: Box<[NeosCmos; 8]> = Box::new(chips.try_into().unwrap());

        let first = chips.first().unwrap();

        let observer = first.observer().clone();
        let side_id = first.side_id;
        let stack_id = first.stack_id;
        let quad_id = first.quad_id;
        let loop_id = first.loop_id;
        let subloop_id = first.subloop_id;
        let exposure_id = first.exposure_id;
        let rotation = first.rotation;

        for ccd in chips.iter() {
            if ccd.side_id != side_id
                || ccd.stack_id != stack_id
                || ccd.quad_id != quad_id
                || ccd.loop_id != loop_id
                || ccd.subloop_id != subloop_id
                || ccd.exposure_id != exposure_id
                || ccd.rotation != rotation
                || ccd.observer().jd != observer.jd
            {
                return Err(NEOSpyError::ValueError(
                    "All NeosCmos must have matching values.".into(),
                ));
            }
        }
        Ok(Self {
            chips,
            observer,
            rotation,
            side_id,
            stack_id,
            quad_id,
            loop_id,
            subloop_id,
            exposure_id,
        })
    }
}

impl FovLike for NeosVisit {
    fn get_fov(&self, index: usize) -> FOV {
        FOV::NeosCmos(self.chips[index].clone())
    }

    fn observer(&self) -> &State {
        &self.observer
    }

    fn contains(&self, obs_to_obj: &Vector3<f64>) -> (usize, Contains) {
        closest_inside(
            &self
                .chips
                .iter()
                .map(|x| x.contains(obs_to_obj).1)
                .collect::<Vec<_>>(),
        )
    }

    fn n_patches(&self) -> usize {
        8
    }
}
