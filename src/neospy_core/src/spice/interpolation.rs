//! Interpolation methods used by Spice SPK Files.
//!
//! It is unlikely to be useful outside of reading these files.
//!
use crate::errors::NEOSpyError;
use nalgebra::DVector;

/// Given a list of chebyshev polynomial coefficients, compute the value of the function
/// and it's derivative.
///
/// This is useful for reading values out of JPL SPK format, be aware though that time
/// scaling is also important for that particular use case.
///
/// # Arguments
///
/// * `t`       - Time at which to evaluate the chebyshev polynomials.
/// * `coef`    - List of coefficients of the chebyshev polynomials.
///
#[inline(always)]
pub fn chebyshev_evaluate_both(x: f64, coef: &[f64]) -> Result<(f64, f64), NEOSpyError> {
    let n_coef = coef.len();

    if n_coef < 2 {
        return Err(NEOSpyError::IOError(
            "File not formatted correctly. Chebyshev polynomial must be greater than order 2."
                .into(),
        ));
    }
    let x2 = 2.0 * x;

    let mut val = 0.0;
    let mut second_t = 1.0;
    let mut last_t = x;
    let mut next_t;

    val += coef[0] * second_t;
    val += coef[1] * last_t;

    // The derivative of the first kind is defined by the recurrence relation:
    // d T_i / dx = i * U_{i-1}
    let mut der_vel;
    let mut second_u = 1.0;
    let mut last_u = x2;
    let mut next_u;

    der_vel = coef[1] * second_u;

    for (idx, &c) in coef.iter().enumerate().skip(2) {
        next_t = x2 * last_t - second_t;
        val += c * next_t;

        second_t = last_t;
        last_t = next_t;

        next_u = x2 * last_u - second_u;
        der_vel += c * last_u * (idx as f64);

        second_u = last_u;
        last_u = next_u;
    }

    Ok((val, der_vel))
}

/// Given a list of chebyshev polynomial coefficients, compute the value of the function
/// using chebyshev polynomials of the first type.
///
/// This is useful for reading values out of JPL SPK format, be aware though that time
/// scaling is also important for that particular use case.
///
/// # Arguments
///
/// * `t`       - Time at which to evaluate the chebyshev polynomials.
/// * `coef`    - List of coefficients of the chebyshev polynomials.
///
#[inline(always)]
pub fn chebyshev_evaluate_type1(t: f64, coef: &[f64]) -> Result<f64, NEOSpyError> {
    let mut val = 0.0;

    let mut second_t = 1.0;
    let mut last_t = t;

    let n_coef = coef.len();

    if n_coef == 0 {
        return Err(NEOSpyError::IOError(
            "File not formatted correctly. Chebyshev polynomial cannot be order 0.".into(),
        ));
    }

    val += coef[0] * second_t;

    if n_coef == 1 {
        return Ok(val);
    }

    val += coef[1] * last_t;

    let mut next_t;

    for &c in coef.iter().skip(2) {
        next_t = 2.0 * t * last_t - second_t;
        val += c * next_t;

        second_t = last_t;
        last_t = next_t;
    }

    Ok(val)
}

/// Interpolate using Hermite interpolation.
///
/// # Arguments
///
/// * `times` - Times where `x` and `d`x are evaluated at.
/// * `x` - The values of the function `f` evaluated at the specified times.
/// * `dx` - The values of the derivative of the function `f`.
/// * `eval_time` - Time at which to evaluate the interpolation function.
pub fn hermite_interpolation(times: &[&f64], x: &[f64], dx: &[f64], eval_time: f64) -> (f64, f64) {
    assert_eq!(times.len(), x.len());
    assert_eq!(times.len(), dx.len());

    let n = x.len();

    let mut work = DVector::<f64>::zeros(2 * x.len());
    let mut d_work = DVector::<f64>::zeros(2 * x.len());
    for (idx, (x0, dx0)) in x.iter().zip(dx).enumerate() {
        work[2 * idx] = *x0;
        work[2 * idx + 1] = *dx0;
    }

    for idx in 1..n {
        let c1 = times[idx] - eval_time;
        let c2 = eval_time - times[idx - 1];
        let denom = times[idx] - times[idx - 1];

        let prev = 2 * idx - 2;
        let cur = prev + 1;
        let next = cur + 1;

        d_work[prev] = work[cur];
        d_work[cur] = (work[next] - work[prev]) / denom;

        let tmp = work[cur] * (eval_time - times[idx - 1]) + work[prev];
        work[cur] = (c1 * work[prev] + c2 * work[next]) / denom;
        work[prev] = tmp;
    }

    d_work[2 * n - 2] = work[2 * n - 1];
    work[2 * n - 2] += work[2 * n - 1] * (eval_time - times[n - 1]);

    for idj in 2..(2 * n) {
        for idi in 1..(2 * n - idj + 1) {
            let xi = (idi + 1) / 2;
            let xij = (idi + idj + 1) / 2;
            let c1 = times[xij - 1] - eval_time;
            let c2 = eval_time - times[xi - 1];
            let denom = times[xij - 1] - times[xi - 1];

            d_work[idi - 1] =
                (c1 * d_work[idi - 1] + c2 * d_work[idi] + work[idi] - work[idi - 1]) / denom;
            work[idi - 1] = (c1 * work[idi - 1] + c2 * work[idi]) / denom;
        }
    }
    (work[0], d_work[0])
}
