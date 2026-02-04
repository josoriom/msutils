//! Savitzky–Golay (Generalized) — smoothing & derivatives (SGG)
//!
//! Rust implementation of JavaScript version from:
//! https://github.com/mljs/savitzky-golay-generalized
//!
//! # References
//! * Peter A. Gorry, (1990).
//! “General Least-Squares Smoothing and Differentiation by the Convolution (Savitzky–Golay) Method.”
//! * MLJS `savitzky-golay-generalized` repository.

#[derive(Clone, Copy, Debug)]
pub struct SggOptions {
    pub window_size: usize,
    pub derivative: usize,
    pub polynomial: usize,
}

impl Default for SggOptions {
    fn default() -> Self {
        Self {
            window_size: 9,
            derivative: 0,
            polynomial: 3,
        }
    }
}

pub fn sgg(ys: &[f64], xs: &[f64], opts: SggOptions) -> Vec<f64> {
    let window_size = opts.window_size;
    let derivative = opts.derivative;
    let polynomial = opts.polynomial;

    if window_size % 2 == 0 || window_size < 5 {
        panic!("Invalid window size (should be odd and at least 5 integer number)");
    }
    if ys.is_empty() {
        panic!("Y values must be an array");
    }
    if xs.is_empty() {
        panic!("X must be defined");
    }
    if window_size > ys.len() {
        panic!(
            "Window size is higher than the data length {}>{}",
            window_size,
            ys.len()
        );
    }
    if polynomial < 1 {
        panic!("Polynomial should be a positive integer");
    }
    if polynomial >= 6 {
        eprintln!(
            "You should not use polynomial grade higher than 5 if you are not sure that your data arises from such a model. Possible polynomial oscillation problems"
        );
    }

    let half = window_size / 2;
    let np = ys.len();
    let mut ans = vec![0.0f64; np];
    let weights = full_weights(window_size, polynomial, derivative);

    let constant_h = xs.len() == 1;
    let mut hs_val = 1.0f64;
    if constant_h {
        hs_val = xs[0].powi(derivative as i32);
    }

    for i in 0..half {
        let wg1 = &weights[half - i - 1];
        let wg2 = &weights[half + i + 1];
        let mut d1 = 0.0f64;
        let mut d2 = 0.0f64;
        for l in 0..window_size {
            d1 += wg1[l] * ys[l];
            d2 += wg2[l] * ys[np - window_size + l];
        }
        if constant_h {
            ans[half - i - 1] = d1 / hs_val;
            ans[np - half + i] = d2 / hs_val;
        } else {
            let hs_left = get_hs(xs, half - i - 1, half, derivative);
            ans[half - i - 1] = d1 / hs_left;
            let hs_right = get_hs(xs, np - half + i, half, derivative);
            ans[np - half + i] = d2 / hs_right;
        }
    }

    let wg = &weights[half];
    for i in window_size..=np {
        let mut d = 0.0f64;
        for l in 0..window_size {
            d += wg[l] * ys[l + i - window_size];
        }
        let denom = if constant_h {
            hs_val
        } else {
            get_hs(xs, i - half - 1, half, derivative)
        };
        ans[i - half - 1] = d / denom;
    }

    ans
}

fn get_hs(h: &[f64], center: usize, half: usize, derivative: usize) -> f64 {
    if derivative == 0 {
        return 1.0;
    }
    let mut hs = 0.0f64;
    let mut count = 0usize;
    let start = center as isize - half as isize;
    let end = center as isize + half as isize;
    let last = h.len().saturating_sub(1) as isize;
    let mut i = start;
    while i < end {
        if i >= 0 && i < last {
            let iu = i as usize;
            hs += h[iu + 1] - h[iu];
            count += 1;
        }
        i += 1;
    }
    if count == 0 {
        1.0
    } else {
        (hs / count as f64).powi(derivative as i32)
    }
}

fn gen_fact(a: i32, b: i32) -> f64 {
    if a >= b {
        let mut gf = 1.0f64;
        let start = a - b + 1;
        let mut j = start;
        while j <= a {
            gf *= j as f64;
            j += 1;
        }
        gf
    } else {
        1.0
    }
}

fn full_weights(m: usize, n: usize, s: usize) -> Vec<Vec<f64>> {
    let half_window = (m / 2) as i32;
    let polynomial_degree = n as i32;
    let derivative_order = s;

    let mut weights = vec![vec![0.0f64; m]; m];

    let mut coefficient_by_order = vec![0.0f64; n + 1];
    for polynomial_order in 0..=polynomial_degree {
        let numerator = gen_fact(2 * half_window, polynomial_order);
        let denominator = gen_fact(2 * half_window + polynomial_order + 1, polynomial_order + 1);
        coefficient_by_order[polynomial_order as usize] =
            (2 * polynomial_order + 1) as f64 * (numerator / denominator);
    }

    let mut gram_tables: Vec<Vec<Vec<f64>>> = Vec::with_capacity(m);
    let mut relative_position = -half_window;
    while relative_position <= half_window {
        gram_tables.push(gram_table(
            relative_position,
            half_window,
            polynomial_degree,
            derivative_order as i32,
        ));
        relative_position += 1;
    }

    let mut row_position = -half_window;
    while row_position <= half_window {
        let row_index = (row_position + half_window) as usize;
        let row_gram = &gram_tables[row_index];

        let mut column_position = -half_window;
        while column_position <= half_window {
            let column_index = (column_position + half_window) as usize;
            let column_gram = &gram_tables[column_index];

            let mut accumulator = 0.0f64;
            for polynomial_order in 0..=polynomial_degree {
                let order_index = polynomial_order as usize;
                accumulator += coefficient_by_order[order_index]
                    * column_gram[order_index][0]
                    * row_gram[order_index][derivative_order];
            }
            weights[row_index][column_index] = accumulator;

            column_position += 1;
        }

        row_position += 1;
    }

    weights
}

fn gram_table(i: i32, m: i32, n: i32, s: i32) -> Vec<Vec<f64>> {
    let max_polynomial_order = n as usize;
    let max_derivative_order = s.max(0) as usize;

    let mut gram = vec![vec![0.0f64; max_derivative_order + 1]; max_polynomial_order + 1];
    gram[0][0] = 1.0;

    let position = i as f64;

    for polynomial_order in 1..=max_polynomial_order {
        let k = polynomial_order as i32;
        let denominator = (k * (2 * m - k + 1)) as f64;

        let a = (4 * k - 2) as f64 / denominator;
        let b = ((k - 1) * (2 * m + k)) as f64 / denominator;

        for derivative_order in 0..=max_derivative_order {
            let mut mixed_term = position * gram[polynomial_order - 1][derivative_order];
            if derivative_order > 0 {
                mixed_term +=
                    (derivative_order as f64) * gram[polynomial_order - 1][derivative_order - 1];
            }

            let two_orders_back = if polynomial_order >= 2 {
                gram[polynomial_order - 2][derivative_order]
            } else {
                0.0
            };

            gram[polynomial_order][derivative_order] = a * mixed_term - b * two_orders_back;
        }
    }

    gram
}
