//! # ZTF Fov definitions.

use super::{closest_inside, Contains, FovLike, OnSkyRectangle, SkyPatch, FOV};
use crate::{
    frames::{Equatorial, Vector},
    prelude::*,
};
use serde::{Deserialize, Serialize};

/// ZTF frame data, single quad of a single chip
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ZtfCcdQuad {
    /// State of the observer
    observer: State<Equatorial>,

    /// Patch of sky
    pub patch: OnSkyRectangle<Equatorial>,

    /// Field ID
    pub field: u32,

    /// File Frac Day
    /// String representation of the filename for this frame.
    pub filefracday: u64,

    /// Magnitude limit of this frame
    pub maglimit: f64,

    /// Filter ID
    pub fid: u64,

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
        corners: [Vector<Equatorial>; 4],
        observer: State<Equatorial>,
        field: u32,
        filefracday: u64,
        ccdid: u8,
        filtercode: Box<str>,
        imgtypecode: Box<str>,
        qid: u8,
        maglimit: f64,
        fid: u64,
    ) -> Self {
        let patch = OnSkyRectangle::from_corners(corners);
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

impl FovLike<Equatorial> for ZtfCcdQuad {
    fn get_fov(&self, index: usize) -> FOV {
        if index != 0 {
            panic!("FOV only has a single patch")
        }
        FOV::ZtfCcdQuad(self.clone())
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

/// ZTF frame data, single quad of a single chip
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ZtfField {
    /// Individual CCD quads
    ccd_quads: Vec<ZtfCcdQuad>,

    /// Observer position
    observer: State<Equatorial>,

    /// Field ID
    pub field: u32,

    /// Filter ID
    pub fid: u64,

    /// Filter code used for the frame
    pub filtercode: Box<str>,

    /// Image Type Code
    pub imgtypecode: Box<str>,
}

impl ZtfField {
    /// Construct a new ZtfField from a list of ccd quads.
    /// These ccd quads must be from the same field and having matching value as
    /// appropriate.
    pub fn new(ccd_quads: Vec<ZtfCcdQuad>) -> Result<Self, NEOSpyError> {
        if ccd_quads.is_empty() {
            return Err(NEOSpyError::ValueError(
                "Ztf Field must contains ZtfCcdQuads".into(),
            ));
        }

        let first = ccd_quads.first().unwrap();

        let observer = first.observer().clone();
        let field = first.field;
        let fid = first.fid;
        let filtercode = first.filtercode.clone();
        let imgtypecode = first.imgtypecode.clone();

        for ccd in ccd_quads.iter() {
            if ccd.field != field
                || ccd.fid != fid
                || ccd.filtercode != filtercode
                || ccd.imgtypecode != imgtypecode
                || ccd.observer().jd != observer.jd
            {
                return Err(NEOSpyError::ValueError(
                    "All ZtfCcdQuads must have matching values except CCD ID etc.".into(),
                ));
            }
        }
        Ok(Self {
            ccd_quads,
            observer,
            field,
            fid,
            filtercode,
            imgtypecode,
        })
    }
}

impl FovLike<Equatorial> for ZtfField {
    fn get_fov(&self, index: usize) -> FOV {
        FOV::ZtfCcdQuad(self.ccd_quads[index].clone())
    }

    fn observer(&self) -> &State<Equatorial> {
        &self.observer
    }

    fn contains(&self, obs_to_obj: &Vector<Equatorial>) -> (usize, Contains) {
        closest_inside(
            &self
                .ccd_quads
                .iter()
                .map(|x| x.contains(obs_to_obj).1)
                .collect::<Vec<_>>(),
        )
    }

    fn n_patches(&self) -> usize {
        self.ccd_quads.len()
    }
}
