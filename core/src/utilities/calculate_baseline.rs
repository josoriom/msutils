use crate::utilities::cheminfo::air_pls::{AirPlsOptions, air_pls};

#[derive(Clone, Copy, Debug)]
pub struct BaselineOptions {
    pub baseline_window: Option<f64>,
    pub baseline_window_factor: Option<usize>,
    pub level: Option<u8>,
}

impl Default for BaselineOptions {
    fn default() -> Self {
        Self {
            baseline_window: Some(20.0),
            baseline_window_factor: Some(4),
            level: Some(5),
        }
    }
}
const EDGE_SLOPE_POINTS: usize = 4;

pub fn calculate_baseline(y: &[f64], options: BaselineOptions) -> Vec<f64> {
    let n = y.len();
    if n == 0 {
        return Vec::new();
    }

    let x: Vec<f64> = (0..n).map(|i| i as f64).collect();

    let defaults = BaselineOptions::default();
    let window = options
        .baseline_window
        .unwrap_or(defaults.baseline_window.unwrap_or(10.0));
    let factor = options
        .baseline_window_factor
        .unwrap_or(defaults.baseline_window_factor.unwrap_or(100));
    let level = options.level.unwrap_or(defaults.level.unwrap_or(3));

    let (ys, weights) = curate_edges(y, level, EDGE_SLOPE_POINTS);

    let result = air_pls(
        &x,
        &ys,
        AirPlsOptions {
            lambda: window,
            max_iterations: factor,
            weights: Some(weights),
            ..Default::default()
        },
    );
    result.baseline
}

fn median(v: &[f64]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    let mut t = v.to_vec();
    t.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = t.len();
    if n % 2 == 1 {
        t[n / 2]
    } else {
        (t[n / 2 - 1] + t[n / 2]) / 2.0
    }
}

fn slope_threshold(y: &[f64], level: u8) -> f64 {
    let n = y.len();
    if n < 3 {
        return 0.0;
    }
    let diffs: Vec<f64> = y.windows(2).map(|w| (w[1] - w[0]).abs()).collect();
    let base = median(&diffs).max(1e-12);
    match level.min(3).max(1) {
        1 => 3.0 * base,
        2 => 2.5 * base,
        _ => 2.0 * base,
    }
}

fn starts_with_tail_k(y: &[f64], k: usize, level: u8) -> bool {
    if y.len() < k {
        return false;
    }
    let s = (y[k - 1] - y[0]) / ((k - 1) as f64);
    s < -slope_threshold(y, level)
}

fn ends_with_tail_k(y: &[f64], k: usize, level: u8) -> bool {
    let n = y.len();
    if n < k {
        return false;
    }
    let s = (y[n - 1] - y[n - k]) / ((k - 1) as f64);
    s > slope_threshold(y, level)
}

fn left_tail_span(y: &[f64], k: usize) -> usize {
    let n = y.len();
    if n < 2 {
        return 0;
    }
    let mut len = 1usize;
    for i in 0..(n - 1) {
        if y[i + 1] - y[i] < 0.0 {
            len += 1;
        } else {
            break;
        }
    }
    len = len.max(k);
    let cap = (n as f64 * 0.25).round() as usize;
    len.min(cap).min(n - 1)
}

fn right_tail_span(y: &[f64], k: usize) -> usize {
    let n = y.len();
    if n < 2 {
        return 0;
    }
    let mut len = 1usize;
    for i in (1..n).rev() {
        if y[i] - y[i - 1] > 0.0 {
            len += 1;
        } else {
            break;
        }
    }
    len = len.max(k);
    let cap = (n as f64 * 0.25).round() as usize;
    len.min(cap).min(n - 1)
}

fn curate_edges(y: &[f64], level: u8, k: usize) -> (Vec<f64>, Vec<f64>) {
    let n = y.len();
    if n == 0 {
        return (Vec::new(), Vec::new());
    }
    let mut ys = y.to_vec();
    let mut wv = vec![1.0; n];

    if starts_with_tail_k(y, k, level) {
        let l = left_tail_span(y, k);
        if l > 0 {
            let v = y[l.min(n - 1)];
            for i in 0..l {
                ys[i] = v;
                wv[i] = 0.0;
            }
        }
    }

    if ends_with_tail_k(y, k, level) {
        let r = right_tail_span(y, k);
        if r > 0 {
            let v = y[n - 1 - r];
            for i in (n - r)..n {
                ys[i] = v;
                wv[i] = 0.0;
            }
        }
    }

    (ys, wv)
}
