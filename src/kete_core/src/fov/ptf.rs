//! # PTF Fov definitions.

use std::{fmt::Display, str::FromStr};

use super::{closest_inside, Contains, FovLike, OnSkyRectangle, SkyPatch, FOV};
use crate::prelude::*;
use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// PTF Filters used over the course of the survey.
#[derive(PartialEq, Clone, Copy, Debug, Serialize, Deserialize)]
pub enum PTFFilter {
    /// G Band Filter
    G,

    /// R Band Filter
    R,

    /// Hydrogen Alpha 656 nm Filter
    HA656,

    /// Hydrogen Alpha 663nm filter
    HA663,
}

impl Display for PTFFilter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PTFFilter::G => f.write_str("G"),
            PTFFilter::R => f.write_str("R"),
            PTFFilter::HA656 => f.write_str("HA656"),
            PTFFilter::HA663 => f.write_str("HA663"),
        }
    }
}

impl FromStr for PTFFilter {
    type Err = Error;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_uppercase().as_str() {
            "G" => Ok(PTFFilter::G),
            "R" => Ok(PTFFilter::R),
            "HA656" => Ok(PTFFilter::HA656),
            "HA663" => Ok(PTFFilter::HA663),
            _ => Err(Error::ValueError(
                "PTF Filter has to be one of ('G', 'R', 'HA656', 'HA663')".into(),
            )),
        }
    }
}

/// PTF frame data, single ccd
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct PtfCcd {
    /// State of the observer
    observer: State,

    /// Patch of sky
    pub patch: OnSkyRectangle,

    /// Field ID
    pub field: u32,

    /// Which CCID was the frame taken with
    pub ccdid: u8,

    /// Filter
    pub filter: PTFFilter,

    /// Filename of the processed image
    pub filename: Box<str>,

    /// Infobits flag
    pub info_bits: u32,

    /// FWHM seeing conditions
    pub seeing: f32,
}

impl PtfCcd {
    /// Create a Ptf field of view
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        corners: [Vector3<f64>; 4],
        observer: State,
        field: u32,
        ccdid: u8,
        filter: PTFFilter,
        filename: Box<str>,
        info_bits: u32,
        seeing: f32,
    ) -> Self {
        let patch = OnSkyRectangle::from_corners(corners, observer.frame);
        Self {
            patch,
            observer,
            field,
            ccdid,
            filter,
            filename,
            seeing,
            info_bits,
        }
    }
}

impl FovLike for PtfCcd {
    fn get_fov(&self, index: usize) -> FOV {
        if index != 0 {
            panic!("FOV only has a single patch")
        }
        FOV::PtfCcd(self.clone())
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

/// Ptf frame data, full collection of all CCDs
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct PtfField {
    /// Individual CCDs
    ccds: Vec<PtfCcd>,

    /// Observer position
    observer: State,

    /// Field ID
    pub field: u32,

    /// Filter
    pub filter: PTFFilter,
}

impl PtfField {
    /// Construct a new PtfField from a list of ccds.
    /// These ccds must be from the same field and having matching value as
    /// appropriate.
    pub fn new(ccds: Vec<PtfCcd>) -> KeteResult<Self> {
        if ccds.is_empty() {
            Err(Error::ValueError("Ptf Field must contains PtfCcd".into()))?;
        }

        let first = ccds.first().unwrap();

        let observer = first.observer().clone();
        let field = first.field;
        let filter = first.filter;

        for ccd in ccds.iter() {
            if ccd.field != field || ccd.filter != filter || ccd.observer().jd != observer.jd {
                Err(Error::ValueError(
                    "All PtfCcds must have matching values except CCD ID etc.".into(),
                ))?;
            }
        }
        Ok(Self {
            ccds,
            observer,
            field,
            filter,
        })
    }
}

impl FovLike for PtfField {
    fn get_fov(&self, index: usize) -> FOV {
        FOV::PtfCcd(self.ccds[index].clone())
    }

    fn observer(&self) -> &State {
        &self.observer
    }

    fn try_frame_change_mut(&mut self, target_frame: Frame) -> KeteResult<()> {
        let _ = self
            .ccds
            .iter_mut()
            .map(|ccd| ccd.try_frame_change_mut(target_frame))
            .collect::<Result<Vec<_>, _>>()?;
        self.observer.try_change_frame_mut(target_frame)?;
        Ok(())
    }

    fn contains(&self, obs_to_obj: &Vector3<f64>) -> (usize, Contains) {
        closest_inside(
            &self
                .ccds
                .iter()
                .map(|x| x.contains(obs_to_obj).1)
                .collect::<Vec<_>>(),
        )
    }

    fn n_patches(&self) -> usize {
        self.ccds.len()
    }
}
