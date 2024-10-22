//! # Definitions of contiguous field of views
//! These field of views are made up of single contiguous patches of sky, typically single image sensors.

use std::fmt::Debug;

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use super::{Contains, FovLike, Frame, KeteResult, OnSkyRectangle, SkyPatch, SphericalCone, FOV};
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
        Self {
            patch,
            observer,
            rotation: f64::NAN,
        }
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

    fn try_frame_change_mut(&mut self, target_frame: Frame) -> KeteResult<()> {
        self.observer.try_change_frame_mut(target_frame)?;
        self.patch = self.patch.try_frame_change(target_frame)?;
        Ok(())
    }
}

/// Generic rectangular FOV
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct OmniDirectional {
    observer: State,
}

impl OmniDirectional {
    /// Create a new Omni-Directional FOV
    pub fn new(observer: State) -> Self {
        Self { observer }
    }
}

impl FovLike for OmniDirectional {
    #[inline]
    fn get_fov(&self, index: usize) -> FOV {
        if index != 0 {
            panic!("FOV only has a single patch")
        }
        FOV::OmniDirectional(self.clone())
    }

    #[inline]
    fn observer(&self) -> &State {
        &self.observer
    }

    #[inline]
    fn contains(&self, _obs_to_obj: &Vector3<f64>) -> (usize, Contains) {
        (0, Contains::Inside)
    }

    #[inline]
    fn n_patches(&self) -> usize {
        1
    }

    fn try_frame_change_mut(&mut self, target_frame: Frame) -> KeteResult<()> {
        self.observer.try_change_frame_mut(target_frame)?;
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

    /// Angle of the cone from the central pointing vector.
    #[inline]
    pub fn angle(&self) -> &f64 {
        &self.patch.angle
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

    fn try_frame_change_mut(&mut self, target_frame: Frame) -> KeteResult<()> {
        self.observer.try_change_frame_mut(target_frame)?;
        self.patch = self.patch.try_frame_change(target_frame)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{self, GMS_SQRT};
    use crate::prelude::*;
    use crate::state::Desig;

    #[test]
    fn test_check_rectangle_visible() {
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
            assert!(fov.check_n_body(&off_state, false).is_ok());

            assert!(fov
                .check_visible(&[off_state], 6.0, false)
                .first()
                .unwrap()
                .is_some());
        }
    }

    /// Test the light delay computations for the different checks
    #[test]
    fn test_check_omni_visible() {
        // Build an observer, and check the observability of ceres with different offsets from the observer time.
        // this will exercise the position, velocity, and time offsets due to light delay.
        let spk = get_spk_singleton().read().unwrap();
        let observer = State::new(
            Desig::Empty,
            2451545.0,
            [0.0, 1., 0.0].into(),
            [-GMS_SQRT, 0.0, 0.0].into(),
            Frame::Ecliptic,
            10,
        );

        for offset in [-10.0, -5.0, 0.0, 5.0, 10.0] {
            let ceres = spk
                .try_get_state(20000001, observer.jd + offset, 10, Frame::Ecliptic)
                .unwrap();

            let fov = OmniDirectional::new(observer.clone());

            // Check two body approximation calculation
            let two_body = fov.check_two_body(&ceres);
            assert!(two_body.is_ok());
            let (_, _, two_body) = two_body.unwrap();
            let dist = (Vector3::from(two_body.pos) - Vector3::from(observer.pos)).norm();
            assert!((observer.jd - two_body.jd - dist * constants::C_AU_PER_DAY_INV).abs() < 1e-6);
            let ceres_exact = spk
                .try_get_state(20000001, two_body.jd, 10, Frame::Ecliptic)
                .unwrap();
            // check that we are within about 150km - not bad for 2 body
            assert!((Vector3::from(two_body.pos) - Vector3::from(ceres_exact.pos)).norm() < 1e-6);

            // Check n body approximation calculation
            let n_body = fov.check_n_body(&ceres, false);
            assert!(n_body.is_ok());
            let (_, _, n_body) = n_body.unwrap();
            assert!((observer.jd - n_body.jd - dist * constants::C_AU_PER_DAY_INV).abs() < 1e-6);
            let ceres_exact = spk
                .try_get_state(20000001, n_body.jd, 10, Frame::Ecliptic)
                .unwrap();
            // check that we are within about 150m
            assert!((Vector3::from(n_body.pos) - Vector3::from(ceres_exact.pos)).norm() < 1e-9);

            // Check spk queries
            let spk_check = &fov.check_spks(&[20000001])[0];
            assert!(spk_check.is_some());
            let spk_check = &spk_check.as_ref().unwrap().states[0];
            assert!((observer.jd - spk_check.jd - dist * constants::C_AU_PER_DAY_INV).abs() < 1e-6);
            let ceres_exact = spk
                .try_get_state(20000001, spk_check.jd, 10, Frame::Ecliptic)
                .unwrap();
            // check that we are within about 150 micron
            assert!((Vector3::from(spk_check.pos) - Vector3::from(ceres_exact.pos)).norm() < 1e-12);

            assert!(fov
                .check_visible(&[ceres], 6.0, false)
                .first()
                .unwrap()
                .is_some());
        }
    }
}
