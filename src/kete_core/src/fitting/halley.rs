//! # Halley's method
//!
//! Third order root finding algorithm.
//! This is the next order method of newton-raphson.

use crate::{errors::Error, prelude::KeteResult};

/// Solve root using Halley's method.
///
/// This accepts a three functions, the first being a single input function for which
/// the root is desired. The second function being the derivative of the first with
/// respect to the input variable. The third is the third derivative.
///
/// ```
///     use kete_core::fitting::halley;
///     let f = |x| { 1.0 * x * x - 1.0 };
///     let d = |x| { 2.0 * x };
///     let dd = |_| { 2.0};
///     let root = halley(f, d, dd, 0.0, 1e-10).unwrap();
///     assert!((root - 1.0).abs() < 1e-12);
/// ```
///
#[inline(always)]
pub fn halley<Func, Der, SecDer>(
    func: Func,
    der: Der,
    sec_der: SecDer,
    start: f64,
    atol: f64,
) -> KeteResult<f64>
where
    Func: Fn(f64) -> f64,
    Der: Fn(f64) -> f64,
    SecDer: Fn(f64) -> f64,
{
    let mut x = start;

    // if the starting position has derivative of 0, nudge it a bit.
    if der(x).abs() < 1e-12 {
        x += 0.1;
    }

    let mut f_eval: f64;
    let mut d_eval: f64;
    let mut dd_eval: f64;
    let mut step: f64;
    for _ in 0..100 {
        f_eval = func(x);
        if f_eval.abs() < atol {
            return Ok(x);
        }
        d_eval = der(x);

        // Derivative is 0, cannot solve
        if d_eval.abs() < 1e-12 {
            Err(Error::Convergence(
                "Halley's root finding failed to converge due to zero derivative.".into(),
            ))?;
        }

        dd_eval = sec_der(x);

        if !dd_eval.is_finite() || !d_eval.is_finite() || !f_eval.is_finite() {
            Err(Error::Convergence(
                "Halley root finding failed to converge due to non-finite evaluations".into(),
            ))?;
        }
        step = f_eval / d_eval;
        step = step / (1.0 - step * dd_eval / (2.0 * d_eval));

        x -= step;
    }
    Err(Error::Convergence(
        "Halley's root finding hit iteration limit without converging.".into(),
    ))?
}

#[cfg(test)]
mod tests {
    use crate::fitting::halley;

    #[test]
    fn test_haley() {
        let f = |x| 1.0 * x * x - 1.0;
        let d = |x| 2.0 * x;
        let dd = |_| 2.0;
        let root = halley(f, d, dd, 0.0, 1e-10).unwrap();
        assert!((root - 1.0).abs() < 1e-12);
    }
}
