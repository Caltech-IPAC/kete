//! # NEOS field of views
use super::{closest_inside, Contains, FovLike, OnSkyRectangle, SkyPatch, FOV};
use crate::constants::{NEOS_HEIGHT, NEOS_WIDTH};
use crate::prelude::*;
use serde::{Deserialize, Serialize};

/// NEOS frame data, a single detector on a single band
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct NeosCmos {
    /// State of the observer
    observer: State<Equatorial>,

    /// Patch of sky
    pub patch: OnSkyRectangle<Equatorial>,

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

    /// CMOS ID
    /// ID number of the CMOS chip, 0, 1, 2, or 3
    pub cmos_id: u8,
}

impl NeosCmos {
    /// Create a NEOS FOV
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        pointing: Vector<Equatorial>,
        rotation: f64,
        observer: State<Equatorial>,
        side_id: u16,
        stack_id: u8,
        quad_id: u8,
        loop_id: u8,
        subloop_id: u8,
        exposure_id: u8,
        cmos_id: u8,
    ) -> Self {
        let patch = OnSkyRectangle::new(pointing, rotation, NEOS_WIDTH, NEOS_HEIGHT);
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

/// NEOS frame data, 4 chips of a visit.
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct NeosVisit {
    /// Individual CMOS fields
    chips: Box<[NeosCmos; 4]>,

    /// Observer position
    observer: State<Equatorial>,

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
    pub fn new(chips: Vec<NeosCmos>) -> KeteResult<Self> {
        if chips.len() != 4 {
            Err(Error::ValueError(
                "Visit must contains 4 NeosCmos fovs".into(),
            ))?;
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
                Err(Error::ValueError(
                    "All NeosCmos must have matching values.".into(),
                ))?;
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

    /// x_width is the longer dimension in radians
    #[allow(clippy::too_many_arguments)]
    pub fn from_pointing(
        x_width: f64,
        y_width: f64,
        gap_angle: f64,
        pointing: Vector<Equatorial>,
        rotation: f64,
        observer: State<Equatorial>,
        side_id: u16,
        stack_id: u8,
        quad_id: u8,
        loop_id: u8,
        subloop_id: u8,
        exposure_id: u8,
    ) -> Self {
        // Rotate the Z axis to match the defined rotation angle, this vector is not
        // orthogonal to the pointing vector, but is in the correct plane of the final
        // up vector.

        let up_vec = Vector::<Equatorial>::new([0.0, 0.0, 1.0]).rotate_around(pointing, -rotation);

        // construct the vector orthogonal to the pointing and rotated z axis vectors.
        let left_vec = pointing.cross(&up_vec);

        // Given the new left vector, and the existing orthogonal pointing vector,
        // construct a new up vector which is in the same plane as it was before, but now
        // orthogonal to the two existing vectors.
        let up_vec = pointing.cross(&left_vec);

        // +-------+-+-------+-+-------+-+-------+   ^
        // |       |g|       |g|       |g|       |   |
        // |   1   |a|   2   |a|   3   |a|   4   |   y
        // |       |p|       |p|       |p|       |   |
        // +-------+-+-------+-+-------+-+-------+   -
        // |            <---- x ----->           |
        //
        // Pointing vector is in the middle of the 'a' in the central gap.

        // the Y direction is bounded by 2 planes, calculate them one time
        let y_top = up_vec.rotate_around(left_vec, y_width / 2.0);
        let y_bottom = (-up_vec).rotate_around(left_vec, -y_width / 2.0);

        let half_gap = gap_angle / 2.0;

        // chip width in the x direction:
        // 4 * chip_width + 3 * gap_angle = x_width
        // chip_width = (x_width - 3 * gap_angle) / 4
        let chip_width = (x_width - 3.0 * gap_angle) / 4.0;

        // for each chip calculate the x bounds
        let chip_1_a = left_vec.rotate_around(up_vec, -x_width / 2.0);
        let chip_1_b = -left_vec.rotate_around(up_vec, -x_width / 2.0 + chip_width);
        let chip_2_a = left_vec.rotate_around(up_vec, -chip_width - half_gap);
        let chip_2_b = -left_vec.rotate_around(up_vec, -half_gap);
        let chip_3_a = left_vec.rotate_around(up_vec, half_gap);
        let chip_3_b = -left_vec.rotate_around(up_vec, chip_width + half_gap);
        let chip_4_a = left_vec.rotate_around(up_vec, x_width / 2.0 - chip_width);
        let chip_4_b = -left_vec.rotate_around(up_vec, x_width / 2.0);

        // make the patches for each chip
        let chip_1_patch = OnSkyRectangle::from_normals([chip_1_a, y_top, chip_1_b, y_bottom]);
        let chip_2_patch = OnSkyRectangle::from_normals([chip_2_a, y_top, chip_2_b, y_bottom]);
        let chip_3_patch = OnSkyRectangle::from_normals([chip_3_a, y_top, chip_3_b, y_bottom]);
        let chip_4_patch = OnSkyRectangle::from_normals([chip_4_a, y_top, chip_4_b, y_bottom]);

        // make the chips
        let chip_1 = NeosCmos {
            observer: observer.clone(),
            patch: chip_1_patch,
            rotation,
            side_id,
            stack_id,
            quad_id,
            loop_id,
            subloop_id,
            exposure_id,
            cmos_id: 0,
        };
        let chip_2 = NeosCmos {
            observer: observer.clone(),
            patch: chip_2_patch,
            rotation,
            side_id,
            stack_id,
            quad_id,
            loop_id,
            subloop_id,
            exposure_id,
            cmos_id: 1,
        };
        let chip_3 = NeosCmos {
            observer: observer.clone(),
            patch: chip_3_patch,
            rotation,
            side_id,
            stack_id,
            quad_id,
            loop_id,
            subloop_id,
            exposure_id,
            cmos_id: 2,
        };
        let chip_4 = NeosCmos {
            observer: observer.clone(),
            patch: chip_4_patch,
            rotation,
            side_id,
            stack_id,
            quad_id,
            loop_id,
            subloop_id,
            exposure_id,
            cmos_id: 3,
        };

        // Put all the chips in a box for safe-keeping, try not to eat them all at once.
        let chips = Box::new([chip_1, chip_2, chip_3, chip_4]);
        Self {
            chips,
            observer,
            rotation,
            side_id,
            stack_id,
            quad_id,
            loop_id,
            subloop_id,
            exposure_id,
        }
    }

    /// Return the central pointing vector.
    pub fn pointing(&self) -> Vector<Equatorial> {
        let mut pointing = Vector::<Equatorial>::new([0.0; 3]);
        self.chips
            .iter()
            .for_each(|chip| pointing = pointing + chip.patch.pointing());
        pointing.normalize()
    }
}

impl FovLike for NeosVisit {
    fn get_fov(&self, index: usize) -> FOV {
        FOV::NeosCmos(self.chips[index].clone())
    }

    fn observer(&self) -> &State<Equatorial> {
        &self.observer
    }

    fn contains(&self, obs_to_obj: &Vector<Equatorial>) -> (usize, Contains) {
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
}
