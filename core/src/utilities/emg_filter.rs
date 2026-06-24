use crate::utilities::cheminfo::lm::{Data2D, LevenbergMarquardtOptions, Weights, lm};
use crate::utilities::cheminfo::sgg::{SggOptions, sgg};
use crate::utilities::functions::{emg_fn, sigma_from_fwhm};

#[derive(Clone, Debug)]
pub struct EmgFilterOptions {
    pub smooth_window: usize,
    pub smooth_polynomial: usize,
    pub min_fit_points: usize,
    pub weight_floor: f64,
    pub fit_iterations: usize,
    pub fit_timeout_seconds: f64,
    pub central_difference: bool,
    pub initial_tails: Vec<f64>,
    pub good_fit: f64,
}

impl Default for EmgFilterOptions {
    fn default() -> Self {
        Self {
            smooth_window: 5,
            smooth_polynomial: 3,
            min_fit_points: 5,
            weight_floor: 0.02,
            fit_iterations: 20,
            fit_timeout_seconds: 0.5,
            central_difference: true,
            initial_tails: vec![0.3, 1.5, 4.0],
            good_fit: 0.9,
        }
    }
}

pub trait Candidate {
    fn bounds(&self) -> (usize, usize);
    fn set_score(&mut self, r2: f64);
}

pub struct EmgFilter {
    min_r2: f64,
    options: EmgFilterOptions,
}

impl EmgFilter {
    pub fn new(min_r2: f64) -> Self {
        Self {
            min_r2,
            options: EmgFilterOptions::default(),
        }
    }

    pub fn with_options(min_r2: f64, options: EmgFilterOptions) -> Self {
        Self { min_r2, options }
    }

    pub fn is_off(&self) -> bool {
        !(self.min_r2 > 0.0)
    }

    pub fn score(&self, x: &[f64], y: &[f64]) -> Option<f64> {
        let options = &self.options;
        if x.len() < options.min_fit_points {
            return None;
        }
        let smoothed = smooth(x, y, options.smooth_window, options.smooth_polynomial);
        let (apex, height) = find_peak_top(&smoothed)?;
        let center = x[apex];
        let width = estimate_fwhm(x, &smoothed, apex, height);
        if !width.is_finite() || width <= 0.0 {
            return None;
        }
        let scaled_x: Vec<f64> = x.iter().map(|value| (value - center) / width).collect();
        let scaled_y: Vec<f64> = smoothed.iter().map(|value| value / height).collect();
        let weight: Vec<f64> = scaled_y
            .iter()
            .map(|value| value.max(options.weight_floor).sqrt())
            .collect();
        let params = fit_emg(&scaled_x, &scaled_y, &weight, options)?;
        let curve = build_curve(&params);
        let r2 = weighted_r_squared(&scaled_x, &scaled_y, &weight, curve.as_ref());
        r2.is_finite().then_some(r2)
    }

    pub fn keeps_peak(&self, x: &[f64], y: &[f64]) -> bool {
        if self.is_off() {
            return true;
        }
        match self.score(x, y) {
            Some(r2) => r2 >= self.min_r2,
            None => true,
        }
    }

    pub fn retain<C: Candidate>(&self, candidates: &mut Vec<C>, x: &[f64], y: &[f64], baseline: &[f64]) {
        if self.is_off() {
            for candidate in candidates.iter_mut() {
                candidate.set_score(f64::NAN);
            }
            return;
        }
        candidates.retain_mut(|candidate| self.keep(candidate, x, y, baseline));
    }

    fn keep<C: Candidate>(&self, candidate: &mut C, x: &[f64], y: &[f64], baseline: &[f64]) -> bool {
        let (start, end) = candidate.bounds();
        if start >= end || end >= y.len() || baseline.len() != y.len() {
            candidate.set_score(f64::NAN);
            return true;
        }
        let peak = subtract_baseline(&y[start..=end], &baseline[start..=end]);
        match self.score(&x[start..=end], &peak) {
            Some(r2) => {
                candidate.set_score(r2);
                r2 >= self.min_r2
            }
            None => {
                candidate.set_score(f64::NAN);
                true
            }
        }
    }
}

fn subtract_baseline(y: &[f64], baseline: &[f64]) -> Vec<f64> {
    y.iter()
        .zip(baseline.iter())
        .map(|(value, base)| (value - base).max(0.0))
        .collect()
}

fn smooth(x: &[f64], y: &[f64], window: usize, polynomial: usize) -> Vec<f64> {
    if y.len() < window {
        return y.to_vec();
    }
    sgg(
        y,
        x,
        SggOptions {
            window_size: window,
            derivative: 0,
            polynomial,
        },
    )
    .iter()
    .map(|value| value.max(0.0))
    .collect()
}

fn find_peak_top(y: &[f64]) -> Option<(usize, f64)> {
    if y.len() < 3 {
        return None;
    }
    let mut apex = 0;
    let mut height = y[0];
    for (index, &value) in y.iter().enumerate() {
        if value > height {
            height = value;
            apex = index;
        }
    }
    (height.is_finite() && height > 0.0).then_some((apex, height))
}

fn estimate_fwhm(x: &[f64], y: &[f64], apex: usize, height: f64) -> f64 {
    let half = 0.5 * height;
    let left = cross_before(x, y, apex, half);
    let right = cross_after(x, y, apex, half);
    match (left, right) {
        (Some(left), Some(right)) => (right - left).abs(),
        (Some(left), None) => 2.0 * (x[apex] - left).abs(),
        (None, Some(right)) => 2.0 * (right - x[apex]).abs(),
        (None, None) => (x[x.len() - 1] - x[0]).abs() / 6.0,
    }
}

fn cross_before(x: &[f64], y: &[f64], apex: usize, half: f64) -> Option<f64> {
    let mut index = apex;
    while index > 0 {
        if y[index] >= half && y[index - 1] < half {
            return Some(crossing(x[index - 1], y[index - 1], x[index], y[index], half));
        }
        index -= 1;
    }
    None
}

fn cross_after(x: &[f64], y: &[f64], apex: usize, half: f64) -> Option<f64> {
    let mut index = apex;
    while index + 1 < y.len() {
        if y[index] >= half && y[index + 1] < half {
            return Some(crossing(x[index], y[index], x[index + 1], y[index + 1], half));
        }
        index += 1;
    }
    None
}

fn crossing(x1: f64, y1: f64, x2: f64, y2: f64, target: f64) -> f64 {
    let slope = y2 - y1;
    if slope.abs() < 1e-18 {
        return x1;
    }
    x1 + (target - y1) / slope * (x2 - x1)
}

fn fit_emg(x: &[f64], y: &[f64], weight: &[f64], options: &EmgFilterOptions) -> Option<Vec<f64>> {
    let mut best: Option<(f64, Vec<f64>)> = None;
    for &initial_tail in &options.initial_tails {
        if let Some((r2, params)) = fit_once(x, y, weight, initial_tail, options) {
            if best.as_ref().is_none_or(|(score, _)| r2 > *score) {
                best = Some((r2, params));
            }
            if r2 >= options.good_fit {
                break;
            }
        }
    }
    best.map(|(_, params)| params)
}

fn fit_once(
    x: &[f64],
    y: &[f64],
    weight: &[f64],
    initial_tail: f64,
    options: &EmgFilterOptions,
) -> Option<(f64, Vec<f64>)> {
    let center = -0.4;
    let height = 1.0 / emg_fn(0.0, 1.0, center, sigma_from_fwhm(1.0), initial_tail).max(1e-6);
    let lm_options = LevenbergMarquardtOptions {
        initial_values: vec![height, center, 1.0, initial_tail],
        min_values: Some(vec![0.0, -4.0, 0.1, 1e-4]),
        max_values: Some(vec![20.0 * height, 0.3, 8.0, 16.0]),
        max_iterations: Some(options.fit_iterations),
        central_difference: Some(options.central_difference),
        timeout: Some(options.fit_timeout_seconds),
        weights: Some(Weights::Array(weight.iter().map(|value| value.sqrt()).collect())),
        ..Default::default()
    };
    let data = Data2D {
        x: x.to_vec(),
        y: y.to_vec(),
    };
    let fit = lm(&data, &build_curve, &lm_options).ok()?;
    let curve = build_curve(&fit.parameter_values);
    Some((weighted_r_squared(x, y, weight, curve.as_ref()), fit.parameter_values))
}

fn build_curve(params: &[f64]) -> Box<dyn Fn(f64) -> f64> {
    let height = params[0];
    let center = params[1];
    let width = params[2].abs().max(1e-9);
    let tail = params[3].abs().max(1e-9);
    let sigma = sigma_from_fwhm(width);
    Box::new(move |rt| emg_fn(rt, height, center, sigma, tail))
}

fn weighted_r_squared(x: &[f64], y: &[f64], weight: &[f64], curve: &dyn Fn(f64) -> f64) -> f64 {
    if y.len() < 3 {
        return 0.0;
    }
    let mean = y.iter().sum::<f64>() / y.len() as f64;
    let mut residual = 0.0;
    let mut total = 0.0;
    for index in 0..y.len() {
        let expected = curve(x[index]);
        residual += weight[index] * (y[index] - expected).powi(2);
        total += weight[index] * (y[index] - mean).powi(2);
    }
    if total <= 1e-18 {
        return 0.0;
    }
    1.0 - residual / total
}
