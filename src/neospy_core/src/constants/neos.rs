/// NEO Surveyor effective wavelength of NC1 and NC2 bands, for stellar sources
pub const NEOS_BANDS: [f64; 2] = [4700.0, 8000.0];

/// NEO Surveyor effective zero mag of NC1 and NC2 bands, for stellar sources
pub const NEOS_ZERO_MAG: [f64; 2] = [170.662, 64.13];

/// Width of wise FOV in radians, 7.10 degrees.
pub const NEOS_WIDTH: f64 = 0.12391837689159739;

/// Flux in the reflection model should be scaled by these values.
pub const NEOS_SUN_CORRECTION: [f64; 2] = [1.0, 1.0];

/// Height of wise FOV in radians, 1.68 degrees.
pub const NEOS_HEIGHT: f64 = 0.029321531433504737;
