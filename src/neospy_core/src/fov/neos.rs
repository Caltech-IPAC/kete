//! # NEOS field of views
use super::{closest_inside, Contains, FovLike, OnSkyRectangle, SkyPatch, FOV};
use crate::constants::{NEOS_HEIGHT, NEOS_WIDTH};
// use crate::frames::rotate_around;
use crate::prelude::*;
use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// NEOS bands
#[derive(Debug, Clone, Copy, Deserialize, Serialize, PartialEq)]
pub enum NeosBand{
    /// No Band defined.
    Undefined,

    /// NEOS NC1 Band.
    NC1,

    /// NEOS NC2 Band.
    NC2,
}

/// Convert a NEOS band from u8
/// 1 is NC1 2 is NC2, everything else is Undefined.
impl From<u8> for NeosBand{
    fn from(value: u8) -> Self {
        match value{
            1=> NeosBand::NC1,
            2=> NeosBand::NC2,
            _=> NeosBand::Undefined
        }
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

    /// Wavelength band
    pub band: NeosBand,

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
        band: NeosBand,
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
        target_frame: Frame,
    ) -> Result<(), NEOSpyError> {
        self.observer.try_change_frame_mut(target_frame)?;
        self.patch = self.patch.try_frame_change(target_frame)?;
        Ok(())
    }
}

/// NEOS frame data, 4 chips of a visit.
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct NeosVisit {
    /// Individual CMOS fields
    chips: Box<[NeosCmos; 4]>,

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

    /// Wavelength band
    pub band: NeosBand,
}

impl NeosVisit {
    /// Construct a new NeosVisit from a list of cmos fovs.
    /// These cmos fovs must be from the same metadata when appropriate.
    pub fn new(chips: Vec<NeosCmos>) -> Result<Self, NEOSpyError> {
        if chips.len() != 4 {
            return Err(NEOSpyError::ValueError(
                "Visit must contains 4 NeosCmos fovs".into(),
            ));
        }
        let chips: Box<[NeosCmos; 4]> = Box::new(chips.try_into().unwrap());

        let first = chips.first().unwrap();

        let observer = first.observer().clone();
        let side_id = first.side_id;
        let stack_id = first.stack_id;
        let quad_id = first.quad_id;
        let loop_id = first.loop_id;
        let subloop_id = first.subloop_id;
        let exposure_id = first.exposure_id;
        let rotation = first.rotation;
        let band = first.band;

        for ccd in chips.iter() {
            if ccd.side_id != side_id
                || ccd.stack_id != stack_id
                || ccd.quad_id != quad_id
                || ccd.loop_id != loop_id
                || ccd.subloop_id != subloop_id
                || ccd.exposure_id != exposure_id
                || ccd.rotation != rotation
                || ccd.observer().jd != observer.jd
                || ccd.band != band
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
            band,
        })
    }

    // /// x_width is the longer dimension in radians
    // pub fn from_pointing(
    //     x_width: f64,
    //     y_width: f64,
    //     gap_fraction: f64,
    //     pointing: Vector3<f64>,
    //     rotation: f64,
    //     observer: State,
    //     side_id: u16,
    //     stack_id: u8,
    //     quad_id: u8,
    //     loop_id: u8,
    //     subloop_id: u8,
    //     exposure_id: u8,
    //     cmos_id: u8,
    //     band: NeosBand,
    // ) -> Self {
    //     let chip_x_width = (1.0 - 3.0 * gap_fraction) * x_width / 4.0;

    //     // Rotate the Z axis to match the defined rotation angle, this vector is not
    //     // orthogonal to the pointing vector, but is in the correct plane of the final
    //     // up vector.
    //     let up_vec = rotate_around(&Vector3::new(0.0, 0.0, 1.0), pointing, -rotation);

    //     // construct the vector orthogonal to the pointing and rotated z axis vectors.
    //     // left = cross(up, pointing)
    //     let left_vec = pointing.cross(&up_vec);

    //     // Given the new left vector, and the existing orthogonal pointing vector,
    //     // construct a new up vector which is in the same plane as it was before, but now
    //     // orthogonal to the two existing vectors.
    //     // up = cross(pointing, left)
    //     let up_vec = pointing.cross(&left_vec);

    //     // +------+-+------+-+------+-+------+   ^
    //     // |  1   |g|  2   |g|  3   |g|  4   |   |
    //     // |      |a|      |a|      |a|      |   y
    //     // |      |p|      |p|      |p|      |   |
    //     // +------+-+------+-+------+-+------+   _
    //     // <-cf->    x ->
    //     //
    //     // pointing vector is in the middle of the 'a' in the central gap.

    //     let outer: Vector3<f64> = rotate_around(&left_vec, up_vec, -lon_width / 2.0);
    //     let n2: Vector3<f64> = rotate_around(&(-left_vec), up_vec, lon_width / 2.0);

    //     let long_top: Vector3<f64> = rotate_around(&up_vec, left_vec, y_width / 2.0);
    //     let long_bottom: Vector3<f64> = rotate_around(&(-up_vec), left_vec, -y_width / 2.0);

    //     // construct the 4 normal vectors
    //     Self {
    //         edge_normals: [n1.into(), n2.into(), n3.into(), n4.into()],
    //         frame,
    //     }
    // }

    /// Return the central pointing vector.
    pub fn pointing(&self) -> Vector3<f64>{
        let mut pointing = Vector3::<f64>::zeros();
        self.chips.iter().for_each(|chip| pointing += *chip.patch.pointing());
        pointing.normalize()
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
        4
    }

    fn try_frame_change_mut(&mut self, target_frame: Frame) -> Result<(), NEOSpyError> {
        let _ = self
            .chips
            .iter_mut()
            .map(|ccd| ccd.try_frame_change_mut(target_frame))
            .collect::<Result<Vec<_>, _>>()?;
        self.observer.try_change_frame_mut(target_frame)?;
        Ok(())
    }
}
