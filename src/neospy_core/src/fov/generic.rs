//! # Definitions of contiguous field of views
//! These field of views are made up of single contiguous patches of sky, typically single image sensors.

use std::fmt::Debug;

use serde::{Deserialize, Serialize};

use super::{Contains, FovLike, OnSkyRectangle, SkyPatch, SphericalCone, FOV};
use crate::{
    frames::{Equatorial, InertialFrame, Vector},
    state::State,
};

/// Generic rectangular FOV
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct GenericRectangle<F: InertialFrame> {
    observer: State<F>,

    /// Patch of sky
    pub patch: OnSkyRectangle<F>,

    /// Rotation of the FOV.
    pub rotation: f64,
}

impl<F: InertialFrame> GenericRectangle<F> {
    /// Create a new Generic Rectangular FOV
    pub fn new(
        pointing: Vector<F>,
        rotation: f64,
        lon_width: f64,
        lat_width: f64,
        observer: State<F>,
    ) -> Self {
        let patch = OnSkyRectangle::new(pointing, rotation, lon_width, lat_width);
        Self {
            observer,
            patch,
            rotation,
        }
    }

    /// Create a Field of view from a collection of corners.
    pub fn from_corners(corners: [Vector<F>; 4], observer: State<F>) -> Self {
        let patch = OnSkyRectangle::from_corners(corners);
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

impl FovLike<Equatorial> for GenericRectangle<Equatorial> {
    #[inline]
    fn get_fov(&self, index: usize) -> FOV {
        if index != 0 {
            panic!("FOV only has a single patch")
        }
        FOV::GenericRectangle(self.clone())
    }

    #[inline]
    fn observer(&self) -> &State<Equatorial> {
        &self.observer
    }

    #[inline]
    fn contains(&self, obs_to_obj: &Vector<Equatorial>) -> (usize, Contains) {
        (0, self.patch.contains(obs_to_obj))
    }

    #[inline]
    fn n_patches(&self) -> usize {
        1
    }
}

/// Generic rectangular FOV
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct GenericCone<F: InertialFrame> {
    observer: State<F>,

    /// Patch of sky
    pub patch: SphericalCone<F>,
}
impl<F: InertialFrame> GenericCone<F> {
    /// Create a new Generic Conic FOV
    pub fn new(pointing: Vector<F>, angle: f64, observer: State<F>) -> Self {
        let patch = SphericalCone::new(&pointing, angle);
        Self { observer, patch }
    }
}

impl FovLike<Equatorial> for GenericCone<Equatorial> {
    #[inline]
    fn get_fov(&self, index: usize) -> FOV {
        if index != 0 {
            panic!("FOV only has a single patch")
        }
        FOV::GenericCone(self.clone())
    }

    #[inline]
    fn observer(&self) -> &State<Equatorial> {
        &self.observer
    }

    #[inline]
    fn contains(&self, obs_to_obj: &Vector<Equatorial>) -> (usize, Contains) {
        (0, self.patch.contains(obs_to_obj))
    }

    #[inline]
    fn n_patches(&self) -> usize {
        1
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::GMS_SQRT;
    use crate::frames::Ecliptic;
    use crate::prelude::*;

    #[test]
    fn test_check_visible() {
        let circular = State::<Ecliptic>::new(
            None,
            2451545.0,
            [0.0, 1., 0.0].into(),
            [-GMS_SQRT, 0.0, 0.0].into(),
            0,
        );
        let circular_back = State::<Ecliptic>::new(
            None,
            2451545.0,
            [1.0, 0.0, 0.0].into(),
            [0.0, GMS_SQRT, 0.0].into(),
            0,
        );

        for offset in [-10.0, -5.0, 0.0, 5.0, 10.0] {
            let off_state = propagate_n_body_spk(
                circular_back.clone().into_frame(),
                circular_back.jd - offset,
                false,
                None,
            )
            .unwrap()
            .into_frame();

            let vec = circular_back.pos - &circular.pos;

            let fov = GenericRectangle::new(
                vec.into_frame(),
                0.0001,
                0.01,
                0.01,
                circular.clone().into_frame(),
            );
            assert!(fov.check_two_body(&off_state.clone().into_frame()).is_ok());
            assert!(fov.check_n_body(&off_state).is_ok());

            assert!(fov
                .check_visible(&[off_state], 6.0)
                .first()
                .unwrap()
                .is_some());
        }
    }
}
