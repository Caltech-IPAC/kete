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
/// This evaluates the coefficients at a single point of time, but for 3 sets of
/// coefficients at once. This is specifically done for performance reasons.
///
/// # Arguments
///
/// * `t`       - Time at which to evaluate the chebyshev polynomials.
/// * `coefx`    - Slice of coefficients of the chebyshev polynomials.
/// * `coefy`    - Slice of coefficients of the chebyshev polynomials.
/// * `coefz`    - Slice of coefficients of the chebyshev polynomials.
///
#[inline(always)]
pub fn chebyshev3_evaluate_both(
    x: f64,
    coefx: &[f64],
    coefy: &[f64],
    coefz: &[f64],
) -> Result<([f64; 3], [f64; 3]), NEOSpyError> {
    let n_coef = coefx.len();

    if n_coef < 2 {
        return Err(NEOSpyError::IOError(
            "File not formatted correctly. Chebyshev polynomial must be greater than order 2."
                .into(),
        ));
    }
    let x2 = 2.0 * x;

    let mut val = [
        coefx[0] + coefx[1] * x,
        coefy[0] + coefy[1] * x,
        coefz[0] + coefz[1] * x,
    ];
    let mut second_t = 1.0;
    let mut last_t = x;
    let mut next_t;

    // The derivative of the first kind is defined by the recurrence relation:
    // d T_i / dx = i * U_{i-1}
    let mut second_u = 1.0;
    let mut last_u = x2;
    let mut next_u;

    let mut der_val = [coefx[1], coefy[1], coefz[1]];

    for (idx, ((x, y), z)) in coefx.iter().zip(coefy).zip(coefz).enumerate().skip(2) {
        next_t = x2 * last_t - second_t;
        val[0] += x * next_t;
        val[1] += y * next_t;
        val[2] += z * next_t;

        second_t = last_t;
        last_t = next_t;

        next_u = x2 * last_u - second_u;
        der_val[0] += x * last_u * (idx as f64);
        der_val[1] += y * last_u * (idx as f64);
        der_val[2] += z * last_u * (idx as f64);

        second_u = last_u;
        last_u = next_u;
    }

    Ok((val, der_val))
}

/// Interpolate using Hermite interpolation.
///
/// # Arguments
///
/// * `times` - Times where `x` and `d`x are evaluated at.
/// * `x` - The values of the function `f` evaluated at the specified times.
/// * `dx` - The values of the derivative of the function `f`.
/// * `eval_time` - Time at which to evaluate the interpolation function.
pub fn hermite_interpolation(times: &[f64], x: &[f64], dx: &[f64], eval_time: f64) -> (f64, f64) {
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
