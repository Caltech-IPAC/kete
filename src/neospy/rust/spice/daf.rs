use neospy_core::spice::DafHeader;
use pyo3::{pyfunction, PyResult};


/// Given a DAF file, return the comments contained within the header.
#[pyfunction]
#[pyo3(name = "daf_header_comments")]
pub fn daf_header_info_py(filename: &str) -> PyResult<String> {
    let mut file = std::fs::File::open(filename)?;
    let daf = DafHeader::try_load_header(&mut file)?;
    Ok(daf.comments)
}
