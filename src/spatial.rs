pub const DEGRA: f64 = std::f64::consts::PI / 180.0;

const RGE: [[f64; 3]; 3] = [
    [-0.054875539, -0.873437105, -0.483834992],
    [0.494109454, -0.444829594, 0.746982249],
    [-0.867666136, -0.198076390, 0.455983795],
];

/// Convert RA and DEC to Galactic coordinates
///
/// # Arguments
///
/// * `ra` - Right Ascension in degrees
/// * `dec` - Declination in degrees
///
/// # Returns
///
/// * `(f64, f64)` - Tuple containing Galactic longitude and latitude in degrees
///
/// # Examples
///
/// ```
/// use flare::spatial::radec2lb;
///
/// let ra = 45.0;
/// let dec = 45.0;
/// let (l, b) = radec2lb(ra, dec);
/// assert_eq!(l, 145.56299769769032);
/// assert_eq!(b, -12.148257544681918);
/// println!("Galactic longitude: {}, Galactic latitude: {}", l, b);
/// ```
pub fn radec2lb(ra: f64, dec: f64) -> (f64, f64) {
    let ra_rad = ra.to_radians();
    let dec_rad = dec.to_radians();
    let u = vec![
        ra_rad.cos() * dec_rad.cos(),
        ra_rad.sin() * dec_rad.cos(),
        dec_rad.sin(),
    ];
    // next do a dot product of RGE and u
    let ug = vec![
        RGE[0][0] * u[0] + RGE[0][1] * u[1] + RGE[0][2] * u[2],
        RGE[1][0] * u[0] + RGE[1][1] * u[1] + RGE[1][2] * u[2],
        RGE[2][0] * u[0] + RGE[2][1] * u[1] + RGE[2][2] * u[2],
    ];
    let x = ug[0];
    let y = ug[1];
    let z = ug[2];
    let galactic_l = y.atan2(x);
    let galactic_b = z.atan2((x * x + y * y).sqrt());
    (galactic_l.to_degrees(), galactic_b.to_degrees())
}

/// Convert degrees to hours, minutes, and seconds
///
/// # Arguments
///
/// * `deg` - Angle in degrees
///
/// # Returns
///
/// * `String` - String representation of the angle in hours, minutes, and seconds
///
/// # Examples
///
/// ```
/// use flare::spatial::deg2hms;
///
/// let deg = 45.0;
/// let hms = deg2hms(deg);
/// assert_eq!(hms, "03:00:00.0000");
/// println!("{}", hms);
/// ```
pub fn deg2hms(deg: f64) -> String {
    if deg <= 0.0 || deg > 360.0 {
        panic!("Invalid RA input: {}", deg);
    }

    let h = deg * 12.0 / 180.0;
    let hours = h.floor() as i32;
    let m = h * 60.0 - hours as f64 * 60.0;
    let minutes = m.floor() as i32;
    let seconds = (m - minutes as f64) * 60.0;
    let hms = format!("{:02.0}:{:02.0}:{:07.4}", hours, minutes, seconds);
    hms
}

/// Convert degrees to degrees, minutes, and seconds
///
/// # Arguments
///
/// * `deg` - Angle in degrees
///
/// # Returns
///
/// * `String` - String representation of the angle in degrees, minutes, and seconds
///
/// # Examples
///
/// ```
/// use flare::spatial::deg2dms;
///
/// let deg = 45.0;
/// let dms = deg2dms(deg);
/// assert_eq!(dms, "45:00:00.000");
/// println!("{}", dms);
/// ```
pub fn deg2dms(deg: f64) -> String {
    if deg <= -90.0 || deg >= 90.0 {
        panic!("Invalid DEC input: {}", deg);
    }

    let degrees = deg.signum() * deg.abs().floor();
    let m = (deg - degrees).abs() * 60.0;
    let minutes = m.floor();
    let seconds = (m - minutes).abs() * 60.0;
    let dms = format!("{:02.0}:{:02.0}:{:06.3}", degrees, minutes, seconds);
    dms
}

/// Calculate the great circle distance between two points on the celestial sphere
///     using the Haversine formula
///
/// # Arguments
///
/// * `ra1_deg` - Right Ascension of the first point in degrees
/// * `dec1_deg` - Declination of the first point in degrees
/// * `ra2_deg` - Right Ascension of the second point in degrees
/// * `dec2_deg` - Declination of the second point in degrees
///
/// # Returns
///
/// * `f64` - Great circle distance in degrees
///
/// # Examples
///
/// ```
/// use flare::spatial::great_circle_distance;
///
/// let ra1 = 45.0;
/// let dec1 = 45.0;
/// let ra2 = 46.0;
/// let dec2 = 46.0;
/// let distance = great_circle_distance(ra1, dec1, ra2, dec2);
/// assert_eq!((distance - 1.221153650840359).abs() < 1e-6, true);
/// println!("{}", distance);
/// ```
pub fn great_circle_distance(ra1_deg: f64, dec1_deg: f64, ra2_deg: f64, dec2_deg: f64) -> f64 {
    let ra1 = ra1_deg * DEGRA;
    let dec1 = dec1_deg * DEGRA;
    let ra2 = ra2_deg * DEGRA;
    let dec2 = dec2_deg * DEGRA;
    let delta_ra = (ra2 - ra1).abs();
    let mut distance = (dec2.cos() * delta_ra.sin()).powi(2)
        + (dec1.cos() * dec2.sin() - dec1.sin() * dec2.cos() * delta_ra.cos()).powi(2);
    distance = distance
        .sqrt()
        .atan2(dec1.sin() * dec2.sin() + dec1.cos() * dec2.cos() * delta_ra.cos());
    distance * 180.0 / std::f64::consts::PI
}

/// Determine if a point is within an ellipse
///
/// # Arguments
///
/// * `alpha` - Right Ascension of the point in degrees
/// * `delta0` - Declination of the point in degrees
/// * `alpha1` - Right Ascension of the center of the ellipse in degrees
/// * `delta01` - Declination of the center of the ellipse in degrees
/// * `d0` - Distance from the center of the ellipse to the edge in degrees
/// * `axis_ratio` - Ratio of the minor axis to the major axis
/// * `pao` - Position angle of the minor axis in degrees
///
/// # Returns
///
/// * `bool` - True if the point is within the ellipse, false otherwise
///
/// # Examples
///
/// ```
/// use flare::spatial::in_ellipse;
///
/// let alpha = 45.0;
/// let delta0 = 45.0;
/// let alpha1 = 46.0;
/// let delta01 = 46.0;
/// let d0 = 1.0;
/// let axis_ratio = 1.0 as f64;
/// let pao = 0 as f64;
/// let inside = in_ellipse(alpha, delta0, alpha1, delta01, d0, axis_ratio, pao);
/// assert_eq!(inside, false);
///
/// let d0 = 1.23;
/// let inside = in_ellipse(alpha, delta0, alpha1, delta01, d0, axis_ratio, pao);
/// assert_eq!(inside, true);
/// ```
pub fn in_ellipse(
    alpha: f64,
    delta0: f64,
    alpha1: f64,
    delta01: f64,
    d0: f64,
    axis_ratio: f64,
    pao: f64,
) -> bool {
    let d_alpha = (alpha1 - alpha) * DEGRA;
    let delta1 = delta01 * DEGRA;
    let delta = delta0 * DEGRA;
    let pa = pao * DEGRA;
    let d = d0 * DEGRA;

    let e = (1.0 - axis_ratio.powi(2)).sqrt();

    let t1 = d_alpha.cos();
    let t22 = d_alpha.sin();
    let t3 = delta1.cos();
    let t32 = delta1.sin();
    let t6 = delta.cos();
    let t26 = delta.sin();
    let t9 = d.cos();
    let t55 = d.sin();

    if t3 * t6 * t1 + t32 * t26 < 0.0 {
        return false;
    }

    let t2 = t1 * t1;
    let t4 = t3 * t3;
    let t5 = t2 * t4;
    let t7 = t6 * t6;
    let t8 = t5 * t7;
    let t10 = t9 * t9;
    let t11 = t7 * t10;
    let t13 = pa.cos();
    let t14 = t13 * t13;
    let t15 = t14 * t10;
    let t18 = t7 * t14;
    let t19 = t18 * t10;

    let t24 = pa.sin();

    let t31 = t1 * t3;

    let t36 = 2.0 * t31 * t32 * t26 * t6;
    let t37 = t31 * t32;
    let t38 = t26 * t6;
    let t45 = t4 * t10;

    let t56 = t55 * t55;
    let t57 = t4 * t7;

    let t60 = -t8 + t5 * t11 + 2.0 * t5 * t15
        - t5 * t19
        - 2.0 * t1 * t4 * t22 * t10 * t24 * t13 * t26
        - t36
        + 2.0 * t37 * t38 * t10
        - 2.0 * t37 * t38 * t15
        - t45 * t14
        - t45 * t2
        + 2.0 * t22 * t3 * t32 * t6 * t24 * t10 * t13
        - t56
        + t7
        - t11
        + t4
        - t57
        + t57 * t10
        + t19
        - t18 * t45;

    let t61 = e * e;
    let t63 = t60 * t61 + t8 + t57 - t4 - t7 + t56 + t36;

    let inside = t63 > 0.0;
    inside
}
