//! # Definitions of contiguous field of views
//! These field of views are made up of single contiguous patches of sky, typically single image sensors.

use std::fmt::Debug;

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use super::{Contains, FovLike, OnSkyRectangle, SkyPatch, SphericalCone, FOV};
use crate::{
    constants::{NEOS_HEIGHT, NEOS_WIDTH, WISE_WIDTH},
    state::State,
};

/// WISE or NEOWISE frame data, all bands
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct WiseCmos {
    /// State of the observer
    observer: State,

    /// Patch of sky
    pub patch: OnSkyRectangle,

    /// Rotation of the FOV.
    pub rotation: f64,

    /// Frame number of the fov
    pub frame_num: usize,

    /// Scan ID of the fov
    pub scan_id: Box<str>,
}

impl WiseCmos {
    /// Create a Wise fov
    pub fn new(
        pointing: Vector3<f64>,
        rotation: f64,
        observer: State,
        frame_num: usize,
        scan_id: Box<str>,
    ) -> Self {
        let patch = OnSkyRectangle::new(pointing, rotation, WISE_WIDTH, WISE_WIDTH, observer.frame);
        Self {
            patch,
            observer,
            frame_num,
            rotation,
            scan_id,
        }
    }
}

impl FovLike for WiseCmos {
    #[inline]
    fn get_fov(&self, index: usize) -> FOV {
        if index != 0 {
            panic!("Wise FOV only has a single patch")
        }
        FOV::Wise(self.clone())
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

    fn try_frame_change_mut(
        &mut self,
        target_frame: super::Frame,
    ) -> Result<(), super::NEOSpyError> {
        self.observer.try_change_frame_mut(target_frame)?;
        self.patch = self.patch.try_frame_change(target_frame)?;
        Ok(())
    }
}

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

    fn try_frame_change_mut(
        &mut self,
        target_frame: super::Frame,
    ) -> Result<(), super::NEOSpyError> {
        self.observer.try_change_frame_mut(target_frame)?;
        self.patch = self.patch.try_frame_change(target_frame)?;
        Ok(())
    }
}

/// ZTF frame data, single quad of a single chip
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ZtfCcdQuad {
    /// State of the observer
    observer: State,

    /// Patch of sky
    pub patch: OnSkyRectangle,

    /// Field ID
    pub field: u32,

    /// File Frac Day
    /// String representation of the filename for this frame.
    pub filefracday: u64,

    /// Magnitude limit of this frame
    pub maglimit: f64,

    /// Filter ID
    pub fid: usize,

    /// Filter code used for the frame
    pub filtercode: Box<str>,

    /// Image Type Code
    pub imgtypecode: Box<str>,

    /// Which CCID was the frame taken with
    pub ccdid: u8,

    /// Quadrant ID
    pub qid: u8,
}

impl ZtfCcdQuad {
    /// Create a ZTF field of view
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        corners: [Vector3<f64>; 4],
        observer: State,
        field: u32,
        filefracday: u64,
        ccdid: u8,
        filtercode: Box<str>,
        imgtypecode: Box<str>,
        qid: u8,
        maglimit: f64,
        fid: usize,
    ) -> Self {
        let patch = OnSkyRectangle::from_corners(corners, observer.frame);
        Self {
            patch,
            observer,
            field,
            filefracday,
            ccdid,
            filtercode,
            imgtypecode,
            qid,
            maglimit,
            fid,
        }
    }
}

impl FovLike for ZtfCcdQuad {
    fn get_fov(&self, index: usize) -> FOV {
        if index != 0 {
            panic!("FOV only has a single patch")
        }
        FOV::ZtfCcdQuad(self.clone())
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

    fn try_frame_change_mut(
        &mut self,
        target_frame: super::Frame,
    ) -> Result<(), super::NEOSpyError> {
        self.observer.try_change_frame_mut(target_frame)?;
        self.patch = self.patch.try_frame_change(target_frame)?;
        Ok(())
    }
}

/// Generic rectangular FOV
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct GenericRectangle {
    observer: State,

    /// Patch of sky
    pub patch: OnSkyRectangle,
}

impl GenericRectangle {
    /// Create a new Generic Rectangular FOV
    pub fn new(
        pointing: Vector3<f64>,
        rotation: f64,
        lon_width: f64,
        lat_width: f64,
        observer: State,
    ) -> Self {
        let patch = OnSkyRectangle::new(pointing, rotation, lon_width, lat_width, observer.frame);
        Self { observer, patch }
    }

    /// Create a Field of view from a collection of corners.
    #[allow(clippy::too_many_arguments)]
    pub fn from_corners(corners: [Vector3<f64>; 4], observer: State) -> Self {
        let patch = OnSkyRectangle::from_corners(corners, observer.frame);
        Self { patch, observer }
    }

    /// Latitudinal width of the FOV.
    #[inline]
    pub fn lat_width(&self) -> f64 {
        self.patch.lat_width()
    }

    /// Longitudinal width of the FOV.
    #[inline]
    pub fn lon_width(&self) -> f64 {
        self.patch.lon_width()
    }
}

impl FovLike for GenericRectangle {
    #[inline]
    fn get_fov(&self, index: usize) -> FOV {
        if index != 0 {
            panic!("FOV only has a single patch")
        }
        FOV::GenericRectangle(self.clone())
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

    fn try_frame_change_mut(
        &mut self,
        target_frame: super::Frame,
    ) -> Result<(), super::NEOSpyError> {
        self.observer.try_change_frame_mut(target_frame)?;
        self.patch = self.patch.try_frame_change(target_frame)?;
        Ok(())
    }
}

/// Generic rectangular FOV
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct GenericCone {
    observer: State,

    /// Patch of sky
    pub patch: SphericalCone,
}
impl GenericCone {
    /// Create a new Generic Conic FOV
    pub fn new(pointing: Vector3<f64>, angle: f64, observer: State) -> Self {
        let patch = SphericalCone::new(&pointing, angle, observer.frame);
        Self { observer, patch }
    }
}

impl FovLike for GenericCone {
    #[inline]
    fn get_fov(&self, index: usize) -> FOV {
        if index != 0 {
            panic!("FOV only has a single patch")
        }
        FOV::GenericCone(self.clone())
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

    fn try_frame_change_mut(
        &mut self,
        target_frame: super::Frame,
    ) -> Result<(), super::NEOSpyError> {
        self.observer.try_change_frame_mut(target_frame)?;
        self.patch = self.patch.try_frame_change(target_frame)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::GMS_SQRT;
    use crate::prelude::*;
    use crate::state::Desig;

    #[test]
    fn test_check_visible() {
        let circular = State::new(
            Desig::Empty,
            2451545.0,
            [0.0, 1., 0.0].into(),
            [-GMS_SQRT, 0.0, 0.0].into(),
            Frame::Ecliptic,
            0,
        );
        let circular_back = State::new(
            Desig::Empty,
            2451545.0,
            [1.0, 0.0, 0.0].into(),
            [0.0, GMS_SQRT, 0.0].into(),
            Frame::Ecliptic,
            0,
        );

        for offset in [-10.0, -5.0, 0.0, 5.0, 10.0] {
            let off_state = propagate_n_body_spk(
                circular_back.clone(),
                circular_back.jd - offset,
                false,
                None,
            )
            .unwrap();

            let vec = Vector3::from(circular_back.pos) - Vector3::from(circular.pos);

            let fov = GenericRectangle::new(vec, 0.0001, 0.01, 0.01, circular.clone());
            assert!(fov.check_two_body(&off_state).is_ok());
            assert!(fov.check_n_body(&off_state).is_ok());

            assert!(fov
                .check_visible(&[off_state], 6.0)
                .first()
                .unwrap()
                .is_some());
        }
    }
}
