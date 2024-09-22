//! # Lambert's Problem
//!
//! Lambert's problem is one of solving an orbit from two positional measurements
//! and the time passage between the two points. This assumes strict two body
//! motion around a central mass, with infinitesimal point mass in orbit.

use std::f64::consts::PI;


pub fn lambert_time(M, Q:f64, X:f64) {
 
    let QSQ = Q*Q;

    let QSQFM1 = 1.0 - QSQ;
    let XSQ = X *X;
    let u = (1.0 - X) * (1.0 + X);

    let y = u.abs().sqrt();

    let z = (QSQFM1 + QSQ * XSQ).sqrt();
    let QX = Q *X;
    if QX.is_sign_negative(){
        let
    }

}

pub fn lambert_time_direct(t:f64, x:f64, q:f64, lm1:bool) -> (f64, f64, f64, f64){
    let q_squared = q.powi(2);
    let x_squared = x.powi(2);
    let u = (1.0 - x) * (1.0 + x);
    let m = 0.0;

    let mut ans: (f64, f64, f64, f64) = (t, 0.0, 0.0, 0.0);
    let y = u.abs().sqrt();
    let z = 1.0 - q_squared + q_squared * x_squared;
    let qx = q * x;

    let mut a = 0.0;
    let mut b = 0.0;
    let mut aa = 0.0;
    let mut bb = 0.0;
    if qx.is_sign_negative(){
        a = z - qx;
        b = q* z - x;

        if lm1{
            aa = (1.0 - q_squared) / a;
            bb = (1.0 - q_squared) * ( q_squared * u - x_squared) / b;
        }
    } 
    if (qx == 0.0 && lm1) || qx.is_sign_positive() {
        aa = z + qx;
        bb = q * z + x;
    }
    if qx > 0.0{
        a = (1.0 - q_squared) / aa;
        b = (1.0 - q_squared) *(q_squared * u - x_squared) / bb;
    }
    if !lm1{
        let g = if qx * u >= 0.0{
            x*z + q*u
        } else{
            (x_squared - q_squared * u) / (x*z - q * u)
        };
        let f = a * y;
        if x <= 1.0{
            ans.0 = m * PI + f.atan2(g);
        } else{
            let fg1 = f / (g + 1.0);
            
        }
    }



    ans

}
