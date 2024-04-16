use crate::{frame::PyFrames, vector::VectorLike};
use itertools::Itertools;
use nalgebra::UnitVector3;
use neospy_core::constants::C_V;
use neospy_core::flux::*;
use neospy_core::prelude::NEOSpyError;
use pyo3::{pyfunction, PyResult};
use rayon::prelude::*;

/// Calculate the visible flux at the observer assuming a convex faceted object made up
/// of a collection of lambertian surfaces.
///
/// Parameters
/// ----------
/// facet_flux :
///     Flux from each of the facets in Jy / Steradian per unit freq.
/// facet_normals :
///     List of surface normal vectors.
/// obs2obj :
///     Vector from the observer to the object in AU.
/// diameter :
///     Diameter of the object in km.
/// emissivity :
///     Emissivity of the lambertian surfaces (between 0 and 1, 0.9 is a good default).
#[pyfunction]
#[pyo3(name = "lambertian_flux")]
pub fn lambertian_flux_py(
    facet_flux: Vec<f64>,
    facet_normals: Vec<VectorLike>,
    obs2obj: VectorLike,
    diameter: f64,
    emissivity: f64,
) -> f64 {
    let obs2obj = obs2obj.into_vec(PyFrames::Ecliptic);
    let facet_normals: Vec<_> = facet_normals
        .into_iter()
        .map(|norm| norm.into_vec(PyFrames::Ecliptic))
        .collect();
    let obs2obj_r = obs2obj.norm();
    let obs2obj = UnitVector3::new_normalize(obs2obj);
    facet_normals
        .into_iter()
        .zip(facet_flux)
        .map(|(normal, flux)| {
            let normal = UnitVector3::new_normalize(normal);
            lambertian_flux(&normal, &obs2obj, &flux, &obs2obj_r, &diameter, &emissivity)
        })
        .sum()
}

/// Return the Solar flux in Jy / Steradian from the 2000 ASTM Standard Extraterrestrial
/// Spectrum Reference E-490-00:
/// <https://www.nrel.gov/grid/solar-resource/spectra-astm-e490.html>
///
/// Returned values are units Janskys / steradian per unit freq.
///
/// Parameters
/// ----------
///
/// dist:
///     Distance the object is from the Sun in au.
/// wavelength :
///     Wavelength in nm.
#[pyfunction]
#[pyo3(name = "solar_flux")]
pub fn solar_flux_py(dist: f64, wavelength: f64) -> PyResult<f64> {
    Ok(solar_flux(dist, wavelength).ok_or(NEOSpyError::ValueError(
        "Query is outside of the range of the dataset".into(),
    ))?)
}

/// Compute the subsolar surface temperature in kelvin.
///
/// Parameters
/// ----------
/// obj2sun :
///     Vector from the object to the sun in AU.
/// geom_albedo :
///     Geometric albedo.
/// g_param :
///     G Phase parameter in the HG system.
/// beaming :
///     Beaming parameter.
/// emissivity :
///     Emissivity.
#[pyfunction]
#[pyo3(name = "sub_solar_temperature")]
pub fn sub_solar_temperature_py(
    obj2sun: VectorLike,
    geom_albedo: f64,
    g_param: f64,
    beaming: f64,
    emissivity: f64,
) -> f64 {
    let obj2sun = obj2sun.into_vec(PyFrames::Ecliptic);
    sub_solar_temperature(&obj2sun, geom_albedo, g_param, beaming, emissivity)
}

/// Compute the black body flux at the specified temperature and wavelength.
///
/// Flux is in units Janskys / steradian per unit freq.
/// This can be converted to being per unit wavelength by multiplying the result of this
/// by `c / wavelength**2`.
///
/// This is limited by the size of float 32s, numbers exceeding 1e-300 will be rounded
/// to 0. This should not raise any warnings, and should be relatively numerically
/// stable and correct over most domains of interest.
///
/// Parameters
/// ----------
/// temps :
///     In units of Kelvin.
/// wavelength :
///     In units of nanometers.
#[pyfunction]
#[pyo3(name = "black_body_flux")]
pub fn black_body_flux_py(temps: Vec<f64>, wavelength: f64) -> Vec<f64> {
    temps
        .into_par_iter()
        .map(|t| black_body_flux(t, wavelength))
        .collect()
}

/// Calculate the temperatures of each facet assuming the NEATM model.
///
/// This is a multi-core operation.
///
/// Parameters
/// ----------
/// facet_normals :
///     List of surface normal vectors.
/// subsolar_temp :
///     Temperature at the sub-solar point in kelvin.
/// obj2sun :
///     Vector from the object to the sun in AU.
#[pyfunction]
#[pyo3(name = "neatm_facet_temps")]
pub fn neatm_facet_temperature_py(
    facet_normals: Vec<VectorLike>,
    subsolar_temp: f64,
    obj2sun: VectorLike,
) -> Vec<f64> {
    let obj2sun = UnitVector3::new_normalize(obj2sun.into_vec(PyFrames::Ecliptic));
    let facet_normals: Vec<_> = facet_normals
        .into_iter()
        .map(|norm| norm.into_vec(PyFrames::Ecliptic))
        .collect();
    facet_normals
        .into_iter()
        .map(|normal| {
            neatm_facet_temperature(
                &UnitVector3::new_normalize(normal),
                &obj2sun,
                &subsolar_temp,
            )
        })
        .collect_vec()
}

/// Calculate the temperatures of each facet assuming the FRM model.
///
/// This is a multi-core operation.
///
/// Parameters
/// ----------
/// facet_normals :
///     List of surface normal vectors.
/// subsolar_temp :
///     Temperature at the sub-solar point in kelvin.
/// obj2sun :
///     Vector from the object to the sun in AU.
#[pyfunction]
#[pyo3(name = "frm_facet_temps")]
pub fn frm_facet_temperature_py(
    facet_normals: Vec<VectorLike>,
    subsolar_temp: f64,
    obj2sun: VectorLike,
) -> Vec<f64> {
    let obj2sun = UnitVector3::new_normalize(obj2sun.into_vec(PyFrames::Ecliptic));
    let facet_normals: Vec<_> = facet_normals
        .into_iter()
        .map(|norm| norm.into_vec(PyFrames::Ecliptic))
        .collect();
    facet_normals
        .into_iter()
        .map(|normal| {
            frm_facet_temperature(&UnitVector3::new_normalize(normal), subsolar_temp, &obj2sun)
        })
        .collect_vec()
}

#[pyfunction]
#[pyo3(name = "neatm_thermal")]
#[allow(clippy::too_many_arguments)]
pub fn neatm_thermal_py(
    sun2obj: VectorLike,
    sun2obs: VectorLike,
    vis_albedo: f64,
    g_param: f64,
    beaming: f64,
    diameter: f64,
    wavelength: f64,
    emissivity: f64,
) -> f64 {
    let sun2obj = sun2obj.into_vec(PyFrames::Ecliptic);
    let sun2obs = sun2obs.into_vec(PyFrames::Ecliptic);

    let hg_params = HGParams::try_fill(
        "".into(),
        g_param,
        None,
        Some(C_V),
        Some(vis_albedo),
        Some(diameter),
    )
    .unwrap();
    let params = NeatmParams {
        obs_bands: ObserverBands::Generic {
            bands: vec![wavelength; 1],
            zero_mags: None,
            solar_correction: vec![1.0],
        },
        band_albedos: vec![0.0; 1],
        hg_params,
        emissivity,
        beaming,
    };
    params.apparent_thermal_flux(&sun2obj, &sun2obs).unwrap()[0]
}

#[pyfunction]
#[pyo3(name = "frm_thermal")]
#[allow(clippy::too_many_arguments)]
pub fn frm_thermal_py(
    sun2obj: VectorLike,
    sun2obs: VectorLike,
    vis_albedo: f64,
    g_param: f64,
    diameter: f64,
    wavelength: f64,
    emissivity: f64,
) -> f64 {
    let sun2obj = sun2obj.into_vec(PyFrames::Ecliptic);
    let sun2obs = sun2obs.into_vec(PyFrames::Ecliptic);
    let hg_params = HGParams::try_fill(
        "".into(),
        g_param,
        None,
        Some(C_V),
        Some(vis_albedo),
        Some(diameter),
    )
    .unwrap();

    let params = FrmParams {
        obs_bands: ObserverBands::Generic {
            bands: vec![wavelength; 1],
            zero_mags: None,
            solar_correction: vec![1.0],
        },
        band_albedos: vec![0.0; 1],
        hg_params,
        emissivity,
    };
    params.apparent_thermal_flux(&sun2obj, &sun2obs).unwrap()[0]
}

#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(name = "hg_apparent_flux")]
pub fn hg_apparent_flux_py(
    sun2obj: VectorLike,
    sun2obs: VectorLike,
    h_mag: f64,
    g_param: f64,
    c_hg: f64,
    diameter: f64,
    wavelength: f64,
    albedo: f64,
) -> f64 {
    let sun2obj = sun2obj.into_vec(PyFrames::Ecliptic);
    let sun2obs = sun2obs.into_vec(PyFrames::Ecliptic);
    let params = HGParams::try_fill(
        "".into(),
        g_param,
        Some(h_mag),
        Some(c_hg),
        None,
        Some(diameter),
    )
    .unwrap();
    params
        .apparent_flux(&sun2obj, &sun2obs, wavelength, albedo)
        .unwrap()
}

#[pyfunction]
#[pyo3(name = "hg_apparent_mag")]
pub fn hg_apparent_mag_py(
    sun2obj: VectorLike,
    sun2obs: VectorLike,
    h_mag: f64,
    g_param: f64,
) -> f64 {
    let sun2obj = sun2obj.into_vec(PyFrames::Ecliptic);
    let sun2obs = sun2obs.into_vec(PyFrames::Ecliptic);
    let params = HGParams::new("".into(), g_param, h_mag, None);
    params.apparent_mag(&sun2obj, &sun2obs)
}

#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(name = "comet_apparent_mags")]
pub fn comet_mags_py(
    sun2obj: VectorLike,
    sun2obs: VectorLike,
    mk_1: Option<[f64; 2]>,
    mk_2: Option<[f64; 2]>,
    phase_corr: Option<f64>,
) -> (Option<f64>, Option<f64>) {
    let sun2obj = sun2obj.into_vec(PyFrames::Ecliptic);
    let sun2obs = sun2obs.into_vec(PyFrames::Ecliptic);
    let corr = phase_corr.unwrap_or(0.035);
    let mk_params = CometMKParams::new("".into(), mk_1, mk_2, corr);
    (
        mk_params.apparent_total_flux(&sun2obs, &sun2obj),
        mk_params.apparent_nuclear_flux(&sun2obs, &sun2obj),
    )
}

#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(name = "neos_neatm_flux")]
pub fn neos_neatm_flux_py(
    sun2obj_vec: Vec<VectorLike>,
    sun2obs_vec: Vec<VectorLike>,
    vis_albedo: f64,
    ir_albedo: f64,
    g_param: f64,
    beaming: f64,
    diameter: f64,
    emissivity: f64,
) -> Vec<[Vec<f64>; 3]> {
    let hg_params = HGParams::try_fill(
        "".into(),
        g_param,
        None,
        Some(C_V),
        Some(vis_albedo),
        Some(diameter),
    )
    .unwrap();

    let params = NeatmParams::new_neos([ir_albedo; 2], beaming, hg_params, emissivity);

    sun2obs_vec
        .into_par_iter()
        .zip(sun2obj_vec)
        .map(|(sun2obs, sun2obj)| {
            let sun2obj = sun2obj.into_vec(PyFrames::Ecliptic);
            let sun2obs = sun2obs.into_vec(PyFrames::Ecliptic);
            let flux_res = params.apparent_total_flux(&sun2obj, &sun2obs).unwrap();
            let refl = flux_res.reflected_fraction();
            [flux_res.fluxes, refl, flux_res.magnitudes.unwrap()]
        })
        .collect()
}

#[pyfunction]
#[pyo3(name = "neos_frm_flux")]
pub fn neos_frm_flux_py(
    sun2obj_vec: Vec<VectorLike>,
    sun2obs_vec: Vec<VectorLike>,
    vis_albedo: f64,
    ir_albedo: f64,
    g_param: f64,
    diameter: f64,
    emissivity: f64,
) -> Vec<[Vec<f64>; 3]> {
    let hg_params = HGParams::try_fill(
        "".into(),
        g_param,
        None,
        Some(C_V),
        Some(vis_albedo),
        Some(diameter),
    )
    .unwrap();
    let params = FrmParams::new_neos([ir_albedo; 2], hg_params, emissivity);
    sun2obs_vec
        .into_par_iter()
        .zip(sun2obj_vec)
        .map(|(sun2obs, sun2obj)| {
            let sun2obj = sun2obj.into_vec(PyFrames::Ecliptic);
            let sun2obs = sun2obs.into_vec(PyFrames::Ecliptic);
            let flux_res = params.apparent_total_flux(&sun2obj, &sun2obs).unwrap();
            let refl = flux_res.reflected_fraction();
            [flux_res.fluxes, refl, flux_res.magnitudes.unwrap()]
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(name = "wise_neatm_flux")]
pub fn wise_neatm_flux_py(
    sun2obj_vec: Vec<VectorLike>,
    sun2obs_vec: Vec<VectorLike>,
    vis_albedo: f64,
    ir_albedo: f64,
    g_param: f64,
    beaming: f64,
    diameter: f64,
    emissivity: f64,
) -> Vec<[Vec<f64>; 3]> {
    let hg_params = HGParams::try_fill(
        "".into(),
        g_param,
        None,
        Some(C_V),
        Some(vis_albedo),
        Some(diameter),
    )
    .unwrap();
    let params = NeatmParams::new_wise([ir_albedo; 4], beaming, hg_params, emissivity);

    sun2obs_vec
        .into_par_iter()
        .zip(sun2obj_vec)
        .map(|(sun2obs, sun2obj)| {
            let sun2obj = sun2obj.into_vec(PyFrames::Ecliptic);
            let sun2obs = sun2obs.into_vec(PyFrames::Ecliptic);
            let flux_res = params.apparent_total_flux(&sun2obj, &sun2obs).unwrap();
            let refl = flux_res.reflected_fraction();
            [flux_res.fluxes, refl, flux_res.magnitudes.unwrap()]
        })
        .collect()
}

#[pyfunction]
#[pyo3(name = "wise_frm_flux")]
pub fn wise_frm_flux_py(
    sun2obj_vec: Vec<VectorLike>,
    sun2obs_vec: Vec<VectorLike>,
    vis_albedo: f64,
    ir_albedo: f64,
    g_param: f64,
    diameter: f64,
    emissivity: f64,
) -> Vec<[Vec<f64>; 3]> {
    let hg_params = HGParams::try_fill(
        "".into(),
        g_param,
        None,
        Some(C_V),
        Some(vis_albedo),
        Some(diameter),
    )
    .unwrap();

    let params = FrmParams::new_wise([ir_albedo; 4], hg_params, emissivity);
    sun2obs_vec
        .into_par_iter()
        .zip(sun2obj_vec)
        .map(|(sun2obs, sun2obj)| {
            let sun2obj = sun2obj.into_vec(PyFrames::Ecliptic);
            let sun2obs = sun2obs.into_vec(PyFrames::Ecliptic);
            let flux_res = params.apparent_total_flux(&sun2obj, &sun2obs).unwrap();
            let refl = flux_res.reflected_fraction();
            [flux_res.fluxes, refl, flux_res.magnitudes.unwrap()]
        })
        .collect()
}

#[pyfunction]
#[pyo3(name = "w1_color_correction")]
pub fn w1_color_correction_py(temps: f64) -> f64 {
    w1_color_correction(temps)
}

#[pyfunction]
#[pyo3(name = "w2_color_correction")]
pub fn w2_color_correction_py(temps: f64) -> f64 {
    w2_color_correction(temps)
}

#[pyfunction]
#[pyo3(name = "w3_color_correction")]
pub fn w3_color_correction_py(temps: f64) -> f64 {
    w3_color_correction(temps)
}

#[pyfunction]
#[pyo3(name = "w4_color_correction")]
pub fn w4_color_correction_py(temps: f64) -> f64 {
    w4_color_correction(temps)
}

#[pyfunction]
#[pyo3(name = "fib_lattice_vecs")]
pub fn fib_lattice_vecs_py(n_facets: usize) -> Vec<[f64; 3]> {
    let facets = ConvexShape::new_fibonacci_lattice(n_facets).facets;
    facets
        .iter()
        .map(|f| f.normal.into_inner().into())
        .collect()
}
