use kete_core::{fitting, stats};
use pyo3::pyfunction;

#[pyfunction]
#[pyo3(name = "ks_test")]
pub fn ks_test_py(sample_a: Vec<f64>, sample_b: Vec<f64>) -> f64 {
    stats::two_sample_ks_statistic(&sample_a, &sample_b)
}

#[pyfunction]
#[pyo3(name = "fit_chi2")]
pub fn fit_chi2_py(data: Vec<f64>, sigmas: Vec<f64>) -> f64 {
    assert_eq!(data.len(), sigmas.len());
    fitting::fit_reduced_chi2(&data, &sigmas)
}
