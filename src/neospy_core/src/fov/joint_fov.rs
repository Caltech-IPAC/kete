//! # Definitios of joint field of views
//! These field of views are made up of multiple contiguous patches of sky, typically full image arrays.
use super::contiguous_fov::ZtfCcdQuad;
use super::{closest_inside, Contains, FovLike, FOV};
use crate::prelude::*;
use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// ZTF frame data, single quad of a single chip
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ZtfField {
    /// Individual CCD quads
    ccd_quads: Vec<ZtfCcdQuad>,

    /// Observer position
    observer: State,

    /// Field ID
    pub field: u32,

    /// Filter ID
    pub fid: usize,

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

impl FovLike for ZtfField {
    fn get_fov(&self, index: usize) -> FOV {
        FOV::ZtfCcdQuad(self.ccd_quads[index].clone())
    }

    fn observer(&self) -> &State {
        &self.observer
    }

    fn contains(&self, obs_to_obj: &Vector3<f64>) -> (usize, Contains) {
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
