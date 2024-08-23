use pyo3::pyfunction;

#[pyfunction]
#[pyo3(name = "load_tle")]
pub fn load_tle(line1: &str, line2: &str) {
    let elm = sgp4::Elements::from_tle(None, line1.as_bytes(), line2.as_bytes()).unwrap();
    dbg!("{}", elm.epoch());
    let con = sgp4::Constants::from_elements(&elm).unwrap();
    dbg!(con);
}
