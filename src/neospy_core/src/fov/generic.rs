//! # Definitions of contiguous field of views
//! These field of views are made up of single contiguous patches of sky, typically single image sensors.

use std::fmt::Debug;

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use super::{Contains, FovLike, OnSkyRectangle, SkyPatch, SphericalCone, FOV, Frame, NEOSpyError};
use crate::state::State;

/// Generic rectangular FOV
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct GenericRectangle {
    observer: State,

    /// Patch of sky
    pub patch: OnSkyRectangle,

    /// Rotation of the FOV.
    pub rotation: f64,
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
        Self {
            observer,
            patch,
            rotation,
        }
    }

    /// Create a Field of view from a collection of corners.
    pub fn from_corners(corners: [Vector3<f64>; 4], observer: State) -> Self {
        let patch = OnSkyRectangle::from_corners(corners, observer.frame);
        Self { patch, observer , rotation:f64::NAN}
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
        target_frame: Frame,
    ) -> Result<(), NEOSpyError> {
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
        target_frame: Frame,
    ) -> Result<(), NEOSpyError> {
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
