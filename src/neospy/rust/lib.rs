use pyo3::prelude::{pymodule, wrap_pyfunction, Bound, PyModule, PyResult, Python};

mod elements;
mod fitting;
mod flux;
mod fovs;
mod frame;
mod horizons;
mod kepler;
mod propagation;
mod simult_states;
mod spice;
mod state;
mod state_transition;
mod vector;

// Due to the nature of this sort of interface, there is quite a bit of boiler-plate
// code which is difficult to avoid.
// There are many structs in this interface which are just thin wrappers over
// the underlying Rust objects. Examples include [`state::PyState`] which is a thin
// wrapper over [`neospy_core::state::State`]. These allow rust objects to exist in
// the python while still being reasonably fast. These wrappers also frequently expose
// the binary save/loading features written into the neospy_core code. Typically these
// wrappers are the same name as the underlying struct, with 'Py' appended on the
// front.

/// Python module which exposes all of the compiled rust functions.
#[pymodule]
fn _rust(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<frame::PyFrames>()?;
    m.add_class::<state::PyState>()?;
    m.add_class::<vector::Vector>()?;
    m.add_class::<elements::PyCometElements>()?;
    m.add_class::<simult_states::PySimultaneousStates>()?;

    m.add_class::<fovs::PyNeosCmos>()?;
    m.add_class::<fovs::PyWiseCmos>()?;
    m.add_class::<fovs::PyZtfCcdQuad>()?;
    m.add_class::<fovs::PyZtfField>()?;
    m.add_class::<fovs::PyGenericRectangle>()?;
    m.add_class::<fovs::FOVList>()?;

    m.add_class::<flux::PyNeatmParams>()?;
    m.add_class::<flux::PyFrmParams>()?;
    m.add_class::<flux::PyModelResults>()?;

    m.add_class::<horizons::HorizonsProperties>()?;
    m.add_class::<horizons::HorizonsCovariance>()?;

    m.add_function(wrap_pyfunction!(frame::frame_change_py, m)?)?;
    m.add_function(wrap_pyfunction!(frame::wgs_lat_lon_to_ecef, m)?)?;

    m.add_function(wrap_pyfunction!(kepler::compute_eccentric_anomaly_py, m)?)?;
    m.add_function(wrap_pyfunction!(kepler::propagation_kepler_py, m)?)?;

    m.add_function(wrap_pyfunction!(propagation::propagation_n_body_spk_py, m)?)?;

    m.add_function(wrap_pyfunction!(fovs::fov_checks_py, m)?)?;
    m.add_function(wrap_pyfunction!(fovs::fov_spk_checks_py, m)?)?;

    m.add_function(wrap_pyfunction!(flux::hg_apparent_flux_py, m)?)?;
    m.add_function(wrap_pyfunction!(flux::hg_apparent_mag_py, m)?)?;
    m.add_function(wrap_pyfunction!(flux::sub_solar_temperature_py, m)?)?;
    m.add_function(wrap_pyfunction!(flux::black_body_flux_py, m)?)?;
    m.add_function(wrap_pyfunction!(flux::neatm_thermal_py, m)?)?;
    m.add_function(wrap_pyfunction!(flux::frm_thermal_py, m)?)?;
    m.add_function(wrap_pyfunction!(flux::neatm_facet_temperature_py, m)?)?;
    m.add_function(wrap_pyfunction!(flux::frm_facet_temperature_py, m)?)?;
    m.add_function(wrap_pyfunction!(flux::lambertian_flux_py, m)?)?;
    m.add_function(wrap_pyfunction!(flux::w1_color_correction_py, m)?)?;
    m.add_function(wrap_pyfunction!(flux::w2_color_correction_py, m)?)?;
    m.add_function(wrap_pyfunction!(flux::w3_color_correction_py, m)?)?;
    m.add_function(wrap_pyfunction!(flux::w4_color_correction_py, m)?)?;
    m.add_function(wrap_pyfunction!(flux::wise_neatm_flux_py, m)?)?;
    m.add_function(wrap_pyfunction!(flux::neos_neatm_flux_py, m)?)?;
    m.add_function(wrap_pyfunction!(flux::wise_frm_flux_py, m)?)?;
    m.add_function(wrap_pyfunction!(flux::neos_frm_flux_py, m)?)?;
    m.add_function(wrap_pyfunction!(flux::comet_mags_py, m)?)?;
    m.add_function(wrap_pyfunction!(flux::fib_lattice_vecs_py, m)?)?;
    m.add_function(wrap_pyfunction!(flux::solar_flux_py, m)?)?;

    m.add_function(wrap_pyfunction!(spice::spk_load_py, m)?)?;
    m.add_function(wrap_pyfunction!(spice::spk_loaded_objects_py, m)?)?;
    m.add_function(wrap_pyfunction!(spice::spk_state_py, m)?)?;
    m.add_function(wrap_pyfunction!(spice::spk_raw_state_py, m)?)?;
    m.add_function(wrap_pyfunction!(spice::spk_reset_py, m)?)?;
    m.add_function(wrap_pyfunction!(spice::spk_get_name_from_id_py, m)?)?;
    m.add_function(wrap_pyfunction!(spice::spk_available_info_py, m)?)?;

    m.add_function(wrap_pyfunction!(spice::pck_reset_py, m)?)?;
    m.add_function(wrap_pyfunction!(spice::pck_load_py, m)?)?;
    m.add_function(wrap_pyfunction!(spice::pck_earth_frame_py, m)?)?;
    m.add_function(wrap_pyfunction!(spice::pck_state_to_earth, m)?)?;

    m.add_function(wrap_pyfunction!(spice::daf_header_info_py, m)?)?;

    m.add_function(wrap_pyfunction!(spice::obs_codes, m)?)?;

    m.add_function(wrap_pyfunction!(state_transition::compute_stm_py, m)?)?;

    m.add_function(wrap_pyfunction!(fitting::ks_test_py, m)?)?;
    m.add_function(wrap_pyfunction!(fitting::fit_chi2_py, m)?)?;

    Ok(())
}