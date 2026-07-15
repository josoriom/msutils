use serde::Serialize;

use crate::utilities::{
    cheminfo::lm::{Data2D, LevenbergMarquardtOptions, lm},
    functions::{emg_fn, gaussian_fn, sigma_from_fwhm},
    structs::{DataXY, Peak},
};

const FIT_ITERATIONS: usize = 40;
const FIT_TIMEOUT_SECONDS: f64 = 1.0;
const MIN_POINTS: usize = 5;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum PeakShape {
    Gaussian,
    EMG,
}

#[derive(Clone, Copy, Debug)]
pub struct PeakSeed {
    pub rt: f64,
    pub intensity: f64,
}

impl From<&Peak> for PeakSeed {
    fn from(peak: &Peak) -> Self {
        Self {
            rt: peak.rt,
            intensity: peak.intensity,
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct PeakParameters {
    pub shape: PeakShape,
    pub height: f64,
    pub center: f64,
    pub fwhm: f64,
    pub tail: f64,
    pub r2: f64,
}

pub fn fit_peak(data: &DataXY, seed: &PeakSeed, shape: PeakShape) -> Option<PeakParameters> {
    if data.x.len() < MIN_POINTS || data.x.len() != data.y.len() || !(seed.intensity > 0.0) {
        return None;
    }
    let width = estimate_fwhm(&data.x, &data.y);
    if !width.is_finite() || width <= 0.0 {
        return None;
    }
    let scaled_x: Vec<f64> = data
        .x
        .iter()
        .map(|value| (value - seed.rt) / width)
        .collect();
    let scaled_y: Vec<f64> = data.y.iter().map(|value| value / seed.intensity).collect();

    let mut params = match shape {
        PeakShape::Gaussian => {
            let fwhm = fit_gaussian(&scaled_x, &scaled_y)?;
            PeakParameters {
                shape: PeakShape::Gaussian,
                height: seed.intensity,
                center: seed.rt,
                fwhm: fwhm * width,
                tail: 0.0,
                r2: 0.0,
            }
        }
        PeakShape::EMG => {
            let fit = fit_emg(&scaled_x, &scaled_y)?;
            PeakParameters {
                shape: PeakShape::EMG,
                height: fit[0] * seed.intensity,
                center: seed.rt + fit[1] * width,
                fwhm: fit[2].abs().max(1e-9) * width,
                tail: fit[3].abs().max(1e-9) * width,
                r2: 0.0,
            }
        }
    };

    let curve = build_curve(&params);
    params.r2 = r_squared(&data.x, &data.y, curve.as_ref());
    Some(params)
}

pub fn draw_peak(data: &DataXY, params: &PeakParameters) -> DataXY {
    let curve = build_curve(params);
    let y = data.x.iter().map(|&value| curve(value)).collect();
    DataXY {
        x: data.x.clone(),
        y,
    }
}

fn build_curve(params: &PeakParameters) -> Box<dyn Fn(f64) -> f64> {
    let height = params.height;
    let center = params.center;
    let sigma = sigma_from_fwhm(params.fwhm.abs().max(1e-9));
    let tail = params.tail.abs().max(1e-9);
    match params.shape {
        PeakShape::Gaussian => Box::new(move |x| gaussian_fn(x, height, center, sigma)),
        PeakShape::EMG => Box::new(move |x| emg_fn(x, height, center, sigma, tail)),
    }
}

fn fit_gaussian(x: &[f64], y: &[f64]) -> Option<f64> {
    let options = LevenbergMarquardtOptions {
        initial_values: vec![1.0],
        min_values: Some(vec![0.05]),
        max_values: Some(vec![20.0]),
        max_iterations: Some(FIT_ITERATIONS),
        central_difference: Some(true),
        timeout: Some(FIT_TIMEOUT_SECONDS),
        ..Default::default()
    };
    let data = Data2D {
        x: x.to_vec(),
        y: y.to_vec(),
    };
    lm(&data, &gaussian_unit, &options)
        .ok()
        .map(|fit| fit.parameter_values[0].abs().max(1e-9))
}

fn fit_emg(x: &[f64], y: &[f64]) -> Option<Vec<f64>> {
    let center = -0.4;
    let tail = 0.3;
    let height = 1.0 / emg_fn(0.0, 1.0, center, sigma_from_fwhm(1.0), tail).max(1e-6);
    let options = LevenbergMarquardtOptions {
        initial_values: vec![height, center, 1.0, tail],
        min_values: Some(vec![0.0, -4.0, 0.1, 1e-4]),
        max_values: Some(vec![20.0 * height, 0.3, 8.0, 8.0]),
        max_iterations: Some(FIT_ITERATIONS),
        central_difference: Some(true),
        timeout: Some(FIT_TIMEOUT_SECONDS),
        ..Default::default()
    };
    let data = Data2D {
        x: x.to_vec(),
        y: y.to_vec(),
    };
    lm(&data, &emg_unit, &options)
        .ok()
        .map(|fit| fit.parameter_values)
}

fn gaussian_unit(params: &[f64]) -> Box<dyn Fn(f64) -> f64> {
    let sigma = sigma_from_fwhm(params[0].abs().max(1e-9));
    Box::new(move |x| gaussian_fn(x, 1.0, 0.0, sigma))
}

fn emg_unit(params: &[f64]) -> Box<dyn Fn(f64) -> f64> {
    let height = params[0];
    let center = params[1];
    let sigma = sigma_from_fwhm(params[2].abs().max(1e-9));
    let tail = params[3].abs().max(1e-9);
    Box::new(move |x| emg_fn(x, height, center, sigma, tail))
}

fn estimate_fwhm(x: &[f64], y: &[f64]) -> f64 {
    let mut apex = 0;
    let mut height = y[0];
    for (index, &value) in y.iter().enumerate() {
        if value > height {
            height = value;
            apex = index;
        }
    }
    if !(height > 0.0) {
        return 0.0;
    }
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
            return Some(crossing(
                x[index - 1],
                y[index - 1],
                x[index],
                y[index],
                half,
            ));
        }
        index -= 1;
    }
    None
}

fn cross_after(x: &[f64], y: &[f64], apex: usize, half: f64) -> Option<f64> {
    let mut index = apex;
    while index + 1 < y.len() {
        if y[index] >= half && y[index + 1] < half {
            return Some(crossing(
                x[index],
                y[index],
                x[index + 1],
                y[index + 1],
                half,
            ));
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

fn r_squared(x: &[f64], y: &[f64], curve: &dyn Fn(f64) -> f64) -> f64 {
    let count = y.len();
    if count < 3 {
        return 0.0;
    }
    let mean = y.iter().sum::<f64>() / count as f64;
    let mut residual = 0.0;
    let mut total = 0.0;
    for index in 0..count {
        let expected = curve(x[index]);
        residual += (y[index] - expected).powi(2);
        total += (y[index] - mean).powi(2);
    }
    if total <= 1e-18 {
        return 0.0;
    }
    1.0 - residual / total
}
