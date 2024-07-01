use neospy_core::{errors::NEOSpyError, propagation::NonGravModel};
use pyo3::{pyclass, pymethods, PyResult};

/// Non-gravitational force models.
///
/// This is used optionally by the N-Body propagation methods to compute orbits
/// including non-gravitational forces, such as solar radiation pressure, or
/// poynting-robertson force.
///
/// There are two generic non-gravitational models available, one is specifically
/// intended for dust modeling, and includes the solar radiation pressure, the other
/// model is a mathematical match to the JPL Horizons comet model.
///
/// See :py:meth:`NonGravModel.new_dust` and :py:meth:`NonGravModel.new_comet` for more
/// details. Note that the Comet model can also represent asteroids which undergo the
/// Yarkovsky effect, see :py:meth:`NonGravModel.new_asteroid`, which is a convenience
/// function over the :py:meth:`NonGravModel.new_comet` method, but with 1/r^2 falloff.
///
#[pyclass(frozen, module = "neospy.propagation", name = "NonGravModel")]
#[derive(Debug, Clone)]
pub struct PyNonGravModel(pub NonGravModel);

#[pymethods]
impl PyNonGravModel {
    #[allow(clippy::new_without_default)]
    #[new]
    pub fn new() -> PyResult<Self> {
        Err(NEOSpyError::ValueError("Non-gravitational force models need to be constructed using either the new_dust, new_comet, or new_asteroid methods.".into()))?
    }

    /// New Dust model
    #[staticmethod]
    pub fn new_dust(radius: f64, density: f64) -> Self {
        Self(NonGravModel::Dust { radius, density })
    }

    /// JPL's non-gravitational forces are modeled as defined on page 139 of the Comets II
    /// textbook.
    ///
    /// This model adds 3 "A" terms to the acceleration which the object feels. These
    /// A terms represent additional radial, tangential, and normal forces on the object.
    ///
    /// The defaults of this method are the defaults that JPL Horizons uses for comets when
    /// they are not otherwise specified.
    ///
    /// accel_additional = A_1 * g(r) * r_vec + A_2 * g(r) * t_vec + A_3 * g(r) * n_vec
    /// Where r_vec, t_vec, n_vec are the radial, tangential, and normal unit vectors for
    /// the object.
    ///
    /// The g(r) function is defined by the equation:
    /// g(r) = alpha (r / r0) ^ -m * (1 + (r / r0) ^ n) ^ -k
    ///
    /// When alpha=1.0, n=0.0, k=0.0, r0=1.0, and m=2.0, this is equivalent to a 1/r^2
    /// correction.
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (a1, a2, a3, alpha=0.111262, r_0=2.808, m=2.15, n=5.093, k=4.6142))]
    #[staticmethod]
    pub fn new_comet(
        a1: f64,
        a2: f64,
        a3: f64,
        alpha: f64,
        r_0: f64,
        m: f64,
        n: f64,
        k: f64,
    ) -> Self {
        Self(NonGravModel::JplComet {
            a1,
            a2,
            a3,
            alpha,
            r_0,
            m,
            n,
            k,
        })
    }

    /// JPL's non-gravitational forces are modeled as defined on page 139 of the Comets II
    /// textbook.
    ///
    /// This model adds 3 "A" terms to the acceleration which the object feels. These
    /// A terms represent additional radial, tangential, and normal forces on the object.
    ///
    /// This is has the default values for a 1/r^2 decay of the A terms.
    ///
    /// accel_additional = A_1 * g(r) * r_vec + A_2 * g(r) * t_vec + A_3 * g(r) * n_vec
    /// Where r_vec, t_vec, n_vec are the radial, tangential, and normal unit vectors for
    /// the object.
    ///
    /// The g(r) function is defined by the equation:
    /// g(r) = alpha (r / r0) ^ -m * (1 + (r / r0) ^ n) ^ -k
    ///
    /// When alpha=1.0, n=0.0, k=0.0, r0=1.0, and m=2.0, this is equivalent to a 1/r^2
    /// correction.
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (a1, a2, a3, alpha=1.0, r_0=2.0, m= 2.0, n=1.0, k=0.0))]
    #[staticmethod]
    pub fn new_asteroid(
        a1: f64,
        a2: f64,
        a3: f64,
        alpha: f64,
        r_0: f64,
        m: f64,
        n: f64,
        k: f64,
    ) -> Self {
        Self(NonGravModel::JplComet {
            a1,
            a2,
            a3,
            alpha,
            r_0,
            m,
            n,
            k,
        })
    }

    pub fn __repr__(&self) -> String {
        match self.0 {
            NonGravModel::Dust { radius, density } => format!(
                "neospy.propagation.NonGravModel.new_dust(radius={:?}, density={:?})",
                radius, density
            ),
            NonGravModel::JplComet {
                a1,
                a2,
                a3,
                alpha,
                r_0,
                m,
                n,
                k,
            } => format!(
                "neospy.propagation.NonGravModel.new_comet(a1={:?}, a2={:?}, a3={:?}, alpha={:?}, r_0={:?}, m={:?}, n={:?}, k={:?})",
                a1, a2, a3, alpha, r_0, m, n, k,
            ),
        }
    }
}
