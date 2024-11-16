use super::{
    common::{black_body_flux, lambertian_vis_scale_factor, sub_solar_temperature, ModelResults},
    flux_to_mag, HGParams, ObserverBands, DEFAULT_SHAPE,
};
use crate::frames::{Equatorial, InertialFrame, UnitVector, Vector};
use crate::{constants::V_MAG_ZERO, io::FileIO};

use serde::{Deserialize, Serialize};

/// Using the NEATM thermal model, calculate the temperature of each facet given the
/// direction of the sun, the subsolar temperature and the facet normal vectors.
///
/// # Arguments
///
/// * `facet_normal` - The facet normal vector, these must be unit length.
/// * `subsolar_temp` - The temperature at the sub-solar point in kelvin.
/// * `obj2sun` - The vector from the object to the sun, unit length.
#[inline(always)]
pub fn neatm_facet_temperature<T: InertialFrame>(
    facet_normal: &UnitVector<T>,
    obj2sun: &UnitVector<T>,
    subsolar_temp: &f64,
) -> f64 {
    let tmp = facet_normal.into_inner().dot(&obj2sun.into_inner());
    if tmp > 0.0 {
        return tmp.sqrt().sqrt() * subsolar_temp;
    }
    0.0
}

/// NEATM input
#[derive(Deserialize, Serialize, Clone, Debug)]
pub struct NeatmParams {
    /// Wavelength band information of the observer.
    pub obs_bands: ObserverBands,

    /// Albedo of the object for each band.
    pub band_albedos: Vec<f64>,

    /// Beaming parameter.
    pub beaming: f64,

    /// HG parameters defining the HG reflected light model.
    pub hg_params: HGParams,

    /// Emissivity of the object.
    pub emissivity: f64,
}

impl FileIO for NeatmParams {}

impl NeatmParams {
    /// Create new NeatmParams with WISE band and zero mags
    pub fn new_wise(albedos: [f64; 4], beaming: f64, hg_params: HGParams, emissivity: f64) -> Self {
        NeatmParams {
            obs_bands: ObserverBands::Wise,
            band_albedos: albedos.to_vec(),
            beaming,
            hg_params,
            emissivity,
        }
    }

    /// Create new NeatmParams with NEOS band and zero mags
    pub fn new_neos(albedos: [f64; 2], beaming: f64, hg_params: HGParams, emissivity: f64) -> Self {
        NeatmParams {
            obs_bands: ObserverBands::Neos,
            band_albedos: albedos.to_vec(),
            beaming,
            hg_params,
            emissivity,
        }
    }

    /// Compute the Flux visible from an object using the NEATM model.
    ///
    /// # Arguments
    ///
    /// * `sun2obj` - Position of the object with respect to the Sun in AU.
    /// * `sun2obs` - Position of the Observer with respect to the Sun in AU.
    /// * `color_correction` - Optional function which defines the color correction. If this
    ///                        is provided, the function must accept a list of temperatures
    ///                        in kelvin and return a scaling factor for how much the flux
    ///                        gets scaled by for that specified temp.
    pub fn apparent_thermal_flux(
        &self,
        sun2obj: &Vector<Equatorial>,
        sun2obs: &Vector<Equatorial>,
    ) -> Option<Vec<f64>> {
        let obj2sun = -sun2obj;
        let obs2obj = sun2obj - sun2obs;
        let obs2obj_r = obs2obj.norm();
        let geom = &DEFAULT_SHAPE;
        let hg_params = &self.hg_params;
        let bands = self.obs_bands.band_wavelength();
        let color_correction = self.obs_bands.color_correction();

        let diameter = hg_params.diam()?;
        let ss_temp = sub_solar_temperature(
            &obj2sun,
            self.hg_params.vis_albedo()?,
            self.hg_params.g_param,
            self.beaming,
            self.emissivity,
        );

        let obj2sun = UnitVector::new_checked(obj2sun);
        let obs2obj = UnitVector::new_checked(obs2obj);

        let mut fluxes = vec![0.0; bands.len()];

        for facet in geom.facets.iter() {
            let temp = neatm_facet_temperature(&facet.normal, &obj2sun, &ss_temp);
            let obs_flux_scaling = lambertian_vis_scale_factor(
                &facet.normal,
                &obs2obj,
                &obs2obj_r,
                &diameter,
                &self.emissivity,
            );
            if temp == 0.0 || obs_flux_scaling == 0.0 {
                continue;
            }
            for (idx, (wavelength, flux)) in bands.iter().zip(&mut fluxes).enumerate() {
                let mut facet_flux = black_body_flux(temp, *wavelength);
                if let Some(funcs) = color_correction {
                    facet_flux *= funcs[idx](temp);
                };
                facet_flux *= facet.area;

                *flux += obs_flux_scaling * facet_flux;
            }
        }
        Some(fluxes)
    }

    /// Compute NEATM with an reflected reflection model added on.
    ///
    /// # Arguments
    ///
    /// * `sun2obj` - Position of the object with respect to the Sun in AU.
    /// * `sun2obs` - Position of the Observer with respect to the Sun in AU.
    pub fn apparent_total_flux(
        &self,
        sun2obj: &Vector<Equatorial>,
        sun2obs: &Vector<Equatorial>,
    ) -> Option<ModelResults> {
        let bands = self.obs_bands.band_wavelength();
        let mut fluxes = vec![0.0; bands.len()];
        let mut hg_fluxes = vec![0.0; bands.len()];

        let thermal_fluxes = self.apparent_thermal_flux(sun2obj, sun2obs)?;
        let sun_correction = self.obs_bands.solar_correction();

        for (idx, (wavelength, sun_corr)) in bands.iter().zip(sun_correction).enumerate() {
            let refl = self.hg_params.apparent_flux(
                sun2obj,
                sun2obs,
                *wavelength,
                self.band_albedos[idx],
            )? * sun_corr;
            hg_fluxes[idx] = refl;
            fluxes[idx] = thermal_fluxes[idx] + refl;
        }

        let v_band_magnitude = self.hg_params.apparent_mag(sun2obj, sun2obs);
        let v_band_flux = flux_to_mag(v_band_magnitude, V_MAG_ZERO);

        let magnitudes: Option<Vec<f64>> = self.obs_bands.zero_mags().map(|mags| {
            fluxes
                .iter()
                .zip(mags)
                .map(|(flux, z_mag)| flux_to_mag(*flux, *z_mag))
                .collect::<Vec<f64>>()
        });

        Some(ModelResults {
            fluxes,
            hg_fluxes,
            thermal_fluxes,
            v_band_magnitude,
            v_band_flux,
            magnitudes,
        })
    }
}

#[cfg(test)]
mod tests {

    use crate::{flux::*, frames::{Equatorial, UnitVector}};
    use std::f64::consts::PI;

    #[test]
    fn test_neatm_facet_temperature() {
        let obj2sun: UnitVector<Equatorial> = UnitVector::new_unchecked([1.0, 0.0, 0.0].into());
        let t = (PI / 4.0).cos().powf(0.25);

        let temp = neatm_facet_temperature(
            &UnitVector::new_unchecked([1.0, 0.0, 0.0].into()),
            &obj2sun,
            &1.0,
        );
        assert!((temp - 1.0).abs() < 1e-8);

        let temp = neatm_facet_temperature(
            &UnitVector::new_unchecked([0.0, 1.0, 0.0].into()),
            &obj2sun,
            &1.0,
        );
        assert!(temp.abs() < 1e-8);

        let temp = neatm_facet_temperature(
            &UnitVector::new_unchecked([-1.0, 0.0, 0.0].into()),
            &obj2sun,
            &1.0,
        );
        assert!(temp.abs() < 1e-8);

        let temp = neatm_facet_temperature(
            &UnitVector::new_checked([1.0, 1.0, 0.0].into()),
            &obj2sun,
            &1.0,
        );
        assert!((temp - t).abs() < 1e-8);

        let temp = neatm_facet_temperature(
            &UnitVector::new_checked([1.0, -1.0, 0.0].into()),
            &obj2sun,
            &1.0,
        );
        assert!((temp - t).abs() < 1e-8);
        let fib_n1024 = ConvexShape::new_fibonacci_lattice(1028);
        let fib_n2048 = ConvexShape::new_fibonacci_lattice(2048);

        // Test with different geometry, answer should converge
        let t1: f64 = fib_n2048
            .facets
            .iter()
            .map(|facet| neatm_facet_temperature(&facet.normal, &obj2sun, &1.0))
            .sum();
        let t2: f64 = fib_n1024
            .facets
            .iter()
            .map(|facet| neatm_facet_temperature(&facet.normal, &obj2sun, &1.0))
            .sum();
        let t1: f64 = t1 / fib_n2048.facets.len() as f64;
        let t2: f64 = t2 / fib_n1024.facets.len() as f64;
        assert!((t1 - t2).abs() < 1e-2);
    }
}
