#[inline]
pub fn is_finite_positive(x: f32) -> bool {
    x.is_finite() && x > 0.0
}

#[inline]
pub fn closest_index(xs: &[f64], v: f64) -> usize {
    if xs.is_empty() {
        return 0;
    }
    let (mut lo, mut hi) = (0usize, xs.len());
    while lo < hi {
        let mid = (lo + hi) / 2;
        if xs[mid] < v {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    if lo == 0 {
        0
    } else if lo >= xs.len() {
        xs.len() - 1
    } else if (v - xs[lo - 1]).abs() <= (xs[lo] - v).abs() {
        lo - 1
    } else {
        lo
    }
}

pub fn xy_integration(x: &[f64], y: &[f64]) -> (f64, f64) {
    let n = x.len();
    if n == 0 || n != y.len() {
        return (0.0, f64::NEG_INFINITY);
    }
    if n == 1 {
        return (0.0, y[0] as f64);
    }
    let mut s = 0.0f64;
    let mut m = y[0];
    for i in 0..(n - 1) {
        let dx: f64 = (x[i + 1] - x[i]) * 60.0;
        let yi = y[i] as f64;
        let yj = y[i + 1] as f64;
        s += dx * (yi + yj) * 0.5;
        if y[i + 1] > m {
            m = y[i + 1];
        }
    }
    (s, m as f64)
}

#[inline]
/// ## Error function
///  https://en.wikipedia.org/wiki/Error_function
/// erfc(x) Approximation with max error of 1.2e-7.
/// ### Function
/// - `t = 1 / (1 + 0.5 * |x|)`
/// - `τ = t * exp( -x^2 - 1.26551223 + t * P(t) )`
/// - `P(t) = c_1 + c_2 t + c_3 * t^2 + … + c_9 * t8`,
/// ## Coefficients:
/// `c = [1.00002368, 0.37409196, 0.09678418, -0.18628806, 0.27886807,
///      -1.13520398, 1.48851587, -0.82215223, 0.17087277]`.
pub fn erfc_approx(x: f64) -> f64 {
    const A: [f64; 9] = [
        1.000_023_68,
        0.374_091_96,
        0.096_784_18,
        -0.186_288_06,
        0.278_868_07,
        -1.135_203_98,
        1.488_515_87,
        -0.822_152_23,
        0.170_872_77,
    ];
    let z = x.abs();
    if z > 26.0 {
        return if x.is_sign_positive() { 0.0 } else { 2.0 };
    }
    let t = 1.0 / (1.0 + 0.5 * z);
    let mut poly = 0.0;
    for &a in A.iter().rev() {
        poly = poly * t + a;
    }
    let ans = t * (-z * z - 1.265_512_23 + t * poly).exp();
    if x >= 0.0 { ans } else { 2.0 - ans }
}
