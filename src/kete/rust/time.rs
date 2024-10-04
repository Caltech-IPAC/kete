use kete_core::{
    errors::Error,
    time::{
        scales::{TAI, TDB, UTC},
        Time,
    },
};
use pyo3::prelude::*;

/// A representation of time, always in JD with TDB scaling.
///
/// Note that TDB is not the same as UTC, there is often about 60 seconds or more
/// offset between these time formats. This class enables fast conversion to and from
/// UTC however, via the :py:meth:`~Time.from_mjd`, and :py:meth:`~Time.from_iso`.
/// UTC can be recovered from this object through :py:meth:`~Time.utc_mjd`,
/// :py:meth:`~Time.utc_jd`, or :py:meth:`~Time.iso`.
///
/// UTC Leap seconds cannot be predicted, as a result of this, UTC becomes a bit fuzzy
/// when attempting to record future times. All conversion of future times ignores the
/// possibility of leap seconds.
///
/// This representation and conversion tools make some small tradeoff for performance
/// vs accuracy. Conversion between time scales is only accurate on the millisecond
/// scale, however internal representation accuracy is on the microsecond scale.
///
/// TDB is treated as equivalent to TT and TCB, because these times only differ by less
/// than milliseconds per century.
///
/// Parameters
/// ----------
/// jd:
///     Julian Date in days.
/// scaling:
///     Accepts 'tdb', 'tai', 'utc', 'tcb', and 'tt', but they are converted to TDB
///     immediately.
#[pyclass(frozen, module = "kete", name = "Time")]
#[derive(Debug)]
pub struct PyTime(Time<TDB>);

impl<'py> FromPyObject<'py> for PyTime {
    fn extract_bound(ob: &Bound<'py, PyAny>) -> PyResult<Self> {
        if let Ok(jd) = ob.extract::<f64>() {
            return Ok(PyTime(Time::new(jd)));
        }
        Ok(PyTime(ob.downcast_exact::<PyTime>()?.get().0))
    }
}

impl From<f64> for PyTime {
    fn from(value: f64) -> Self {
        PyTime(Time::new(value))
    }
}

#[pymethods]
impl PyTime {
    #[new]
    #[pyo3(signature = (jd, scaling="tdb"))]
    pub fn new(jd: f64, scaling: &str) -> PyResult<Self> {
        let scaling = scaling.to_lowercase();

        Ok(match scaling.as_str() {
            "tt" => PyTime(Time::<TDB>::new(jd)),
            "tdb" => PyTime(Time::<TDB>::new(jd)),
            "tcb" => PyTime(Time::<TDB>::new(jd)),
            "tai" => PyTime(Time::<TAI>::new(jd).tdb()),
            "utc" => PyTime(Time::<UTC>::new(jd).tdb()),
            s => Err(Error::ValueError(format!(
                "Scaling of type ({:?}) is not supported, must be one of: 'tt', 'tdb', 'tcb', 'tai', 'utc'",
                s
            )))?,
        })
    }

    /// Time from a modified julian date.
    ///
    /// Parameters
    /// ----------
    /// mjd:
    ///     Modified Julian Date in days.
    /// scaling:
    ///     Accepts 'tdb', 'tai', 'utc', and 'tt', but they are converted to TDB
    ///     immediately.
    #[staticmethod]
    #[pyo3(signature = (mjd, scaling="tdb"))]
    pub fn from_mjd(mjd: f64, scaling: &str) -> PyResult<Self> {
        let scaling = scaling.to_lowercase();

        Ok(match scaling.as_str() {
            "tt" => PyTime(Time::<TDB>::from_mjd(mjd)),
            "tdb" => PyTime(Time::<TDB>::from_mjd(mjd)),
            "tai" => PyTime(Time::<TAI>::from_mjd(mjd).tdb()),
            "utc" => PyTime(Time::<UTC>::from_mjd(mjd).tdb()),
            s => Err(Error::ValueError(format!(
                "Scaling of type ({:?}) is not supported, must be one of: 'tt', 'tdb', 'tai', 'utc'",
                s
            )))?,
        })
    }

    /// Time from an ISO formatted string.
    ///
    /// ISO formatted strings are assumed to be in UTC time scaling.
    ///
    /// Parameters
    /// ----------
    /// s:
    ///     ISO Formatted String.
    #[staticmethod]
    pub fn from_iso(s: &str) -> PyResult<Self> {
        Ok(PyTime(Time::<UTC>::from_iso(s)?.tdb()))
    }

    /// Create time object from the Year, Month, and Day.
    ///
    /// These times are assumed to be in UTC amd conversion is performed automatically.
    ///
    /// Parameters
    /// ----------
    /// year:
    ///     The Year, for example `2020`
    /// month:
    ///     The Month as an integer, 0 = January etc.
    /// day:
    ///     The day as an integer or float.
    #[staticmethod]
    pub fn from_ymd(year: i64, month: u32, day: f64) -> Self {
        let frac_day = day.rem_euclid(1.0);
        let day = day.div_euclid(1.0) as u32;
        PyTime(Time::<UTC>::from_year_month_day(year, month, day, frac_day).tdb())
    }

    /// Time in the current time.
    #[staticmethod]
    pub fn now() -> Self {
        PyTime(Time::<UTC>::now().unwrap().tdb())
    }

    /// Return (year, month, day), where day is a float.
    ///
    /// >>> kete.Time.from_ymd(2010, 1, 1).ymd
    /// (2010, 1, 1.0)
    #[getter]
    pub fn ymd(&self) -> (i64, u32, f64) {
        let (y, m, d, f) = self.0.utc().year_month_day();
        (y, m, d as f64 + f)
    }

    /// Julian Date in TDB scaled time.
    /// The difference between TT and TDB is never more than a few milliseconds
    /// per century, so these are treated as equivalent.
    #[getter]
    pub fn jd(&self) -> f64 {
        self.0.jd
    }

    /// Modified Julian Date in TDB scaled time.
    /// The difference between TT and TDB is never more than a few milliseconds
    /// per century, so these are treated as equivalent.
    #[getter]
    pub fn mjd(&self) -> f64 {
        self.0.mjd()
    }

    /// Julian Date in UTC scaled time.
    #[getter]
    pub fn utc_jd(&self) -> f64 {
        self.0.utc().jd
    }

    /// Modified Julian Date in UTC scaled time.
    #[getter]
    pub fn utc_mjd(&self) -> f64 {
        self.0.utc().mjd()
    }

    /// Time in the UTC ISO time format.
    #[getter]
    pub fn iso(&self) -> PyResult<String> {
        Ok(self.0.utc().to_iso()?)
    }

    #[staticmethod]
    pub fn j2000() -> Self {
        PyTime(Time::<TDB>::new(2451545.0))
    }

    fn __repr__(&self) -> String {
        format!("Time({})", self.0.jd)
    }
}
