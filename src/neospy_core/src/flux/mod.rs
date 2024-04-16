//! # Flux
//! Flux calculations including thermal and reflected light models.
//!
//! There are a few flux calculation models contained here:
//! [`hg_apparent_flux`] - Apparent flux of an object in Jy using the HG system.
//! [`hg_absolute_to_apparent_mag`] - Apparent Visible Mag of an object using the HG system.
//! [`neatm::neatm`] - The NEATM thermal model to compute black body flux in Jy.
//! [`frm::frm`] - The FRM thermal model to compute black body flux in Jy.
//!
//! Hybrid models:
//! [`neatm_refl_flux`] - Combination of the NEATM with the HG apparent flux.
//! [`frm_refl_flux`] - Combination of the FRM with the HG apparent flux.
//!
//! If working with NEOS in the near infrared, the hybrid models should be used.
//!
//! If working only in the visible, then the hg models should suffice.
mod comets;
mod common;
mod frm;
mod neatm;
mod neos;
mod reflected;
mod shapes;
mod sun;
mod wise;

pub use comets::*;
pub use common::*;
pub use frm::*;
pub use neatm::*;
pub use neos::*;
pub use reflected::*;
pub use shapes::*;
pub use sun::{solar_flux, solar_flux_black_body};
pub use wise::*;
