//! # Universal Constants

/// Standard Gravitational Constants of the Sun
/// AU^3 / (Day^2 * Solar Mass)
pub const GMS: f64 = 0.00029591220828411956;

/// Gaussian gravitational constant, equivalent to sqrt of GMS.
/// AU^(3/2) per (Day sqrt(Solar Mass))
pub const GMS_SQRT: f64 = 0.01720209894846;

/// effective surface temperature of the sun in K.
pub const SUN_TEMP: f64 = 5778.0;

/// Sun diameter in km
pub const SUN_DIAMETER: f64 = 1.3914e6;

/// Sun J2 Parameter
///
/// This paper below a source, however there are several papers which all put
/// the Sun's J2 at 2.2e-7.
///
/// "Prospects of Dynamical Determination of General Relativity Parameter β and Solar
/// Quadrupole Moment J2 with Asteroid Radar Astronomy"
/// The Astrophysical Journal, 845:166 (5pp), 2017 August 20
pub const SUN_J2: f64 = 2.2e-7;

/// W/m^2 measured at 1 AU.
pub const SOLAR_FLUX: f64 = 1360.8;

/// Speed of light in AU / Day
pub const C_AU_PER_DAY: f64 = 173.14463267424034;

/// Speed of light in meters / second (Definition)
pub const C_M_PER_S: f64 = 299792458.0;

/// Inverse of the Speed of light in Day / AU
pub const C_AU_PER_DAY_INV: f64 = 0.005775518331436995;

/// Square of the inverse of the speed of light in Day^2 / AU^2
pub const C_AU_PER_DAY_INV_SQUARED: f64 = 3.33566119967647e-05;

/// Km per AU (Definition)
pub const AU_KM: f64 = 149597870.7;

/// Stefan Boltzmann constant in W / (m^2 * K^4)
pub const STEFAN_BOLTZMANN: f64 = 5.670374419e-8;

/// Golden Ratio, unit-less
pub const GOLDEN_RATIO: f64 = 1.618033988749894;

/// Zero point V band magnitude in Jy. This is approximately equivalent to Vega's flux.
/// <https://coolwiki.ipac.caltech.edu/index.php/Central_wavelengths_and_zero_points>
pub const V_MAG_ZERO: f64 = 3597.28;

/// The ID value for the massive objects which needs to be numerically evaluated for
/// the effects of gravity.
///
/// Recorded values are:
/// (ID, GM of the object in solar masses, Radius of the object in AU)
pub const MASSIVE_OBJECTS: &[(i64, f64, f64)] = &[
    (10, 1.0, 0.004654758765894654),                    // Sun
    (1, 1.66012082548908e-07, 1.63139354095098e-05),    // Mercury Barycenter
    (2, 2.44783828779694e-06, 4.04537843465442e-05),    // Venus Barycenter
    (399, 3.00348961546514e-06, 4.33036684926559e-05),  // Earth, not 3 because earth != barycenter
    (301, 3.69430335010988e-08, 1.16138016662292e-05),  // Moon
    (4, 3.2271560829139e-07, 2.27021279387769e-05),     // Mars Barycenter
    (5, 0.0009547919099414246, 0.00047789450254521576), // Jupiter Barycenter
    (6, 0.00028588567002459455, 0.0004028666966848747), // Saturn Barycenter
    (7, 4.36624961322212e-05, 0.0001708513622580592),   // Uranus Barycenter
    (8, 5.15138377265457e-05, 0.0001655371154958558),   // Neptune Barycenter
];

/// Recorded values are:
/// (ID, GM of the object in solar masses, Radius of the object in AU)
pub const MASSIVE_OBJECTS_EXTENDED: &[(i64, f64, f64)] = &[
    (10, 1.0, 0.004654758765894654),                    // Sun
    (1, 1.66012082548908e-07, 1.63139354095098e-05),    // Mercury Barycenter
    (2, 2.44783828779694e-06, 4.04537843465442e-05),    // Venus Barycenter
    (399, 3.00348961546514e-06, 4.33036684926559e-05),  // Earth, not 3 because earth != barycenter
    (301, 3.69430335010988e-08, 1.16138016662292e-05),  // Moon
    (4, 3.2271560829139e-07, 2.27021279387769e-05),     // Mars Barycenter
    (5, 0.0009547919099414246, 0.00047789450254521576), // Jupiter Barycenter
    (6, 0.00028588567002459455, 0.0004028666966848747), // Saturn Barycenter
    (7, 4.36624961322212e-05, 0.0001708513622580592),   // Uranus Barycenter
    (8, 5.15138377265457e-05, 0.0001655371154958558),   // Neptune Barycenter
    (20000001, 4.7191659919436e-10, 6.276827e-6),       // Ceres
    (20000002, 1.0259143964958e-10, 3.415824e-6),       // Pallas
    (20000004, 1.3026452498655e-10, 3.512082e-6),       // Vesta
    (20000010, 4.3953391300849e-11, 2.894426e-6),       // Hygiea
    (20000704, 1.911e-11, 2.219283e-6),                 // Interamnia
];

/// V-band constant for the relationship between D, H_V, and p_v, in km.
pub const C_V: f64 = 1329.0;
