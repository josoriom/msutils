pub use crate::utilities::fit_peak::PeakShape;
use crate::utilities::{
    fit_peak::{PeakSeed, fit_peak},
    structs::DataXY,
};

const MIN_FIT_POINTS: usize = 5;
const FIT_HALF_WIDTH_FWHM: f64 = 1.15;
const MIN_FLAT_TOP_POINTS: usize = 4;
const MIN_FLAT_TOP_SHARE: f64 = 0.2;

fn longest_run_at_maximum(y: &[f64], maximum: f64) -> usize {
    let tolerance = maximum.abs() * 1e-9;
    let mut longest = 0usize;
    let mut run = 0usize;
    for value in y {
        if (maximum - value).abs() <= tolerance {
            run += 1;
            if run > longest {
                longest = run;
            }
        } else {
            run = 0;
        }
    }
    longest
}

fn is_flat_topped(y: &[f64]) -> bool {
    let maximum = y.iter().copied().fold(f64::MIN, f64::max);
    if !(maximum > 0.0) {
        return false;
    }
    let plateau = longest_run_at_maximum(y, maximum);
    if plateau < MIN_FLAT_TOP_POINTS {
        return false;
    }
    let above_half = y.iter().filter(|value| **value >= 0.5 * maximum).count();
    plateau as f64 >= MIN_FLAT_TOP_SHARE * above_half.max(1) as f64
}

pub trait Candidate {
    fn bounds(&self) -> (usize, usize);
    fn set_score(&mut self, r2: f64);
}

pub struct ShapeFilter {
    min_r2: f64,
    shape: PeakShape,
}

impl ShapeFilter {
    pub fn new(min_r2: f64, shape: PeakShape) -> Self {
        Self { min_r2, shape }
    }

    pub fn is_off(&self) -> bool {
        !(self.min_r2 > 0.0)
    }

    pub fn score(&self, x: &[f64], y: &[f64]) -> Option<f64> {
        if x.len() < MIN_FIT_POINTS || x.len() != y.len() {
            return None;
        }
        let finite_y: Vec<f64> = y
            .iter()
            .map(|value| if value.is_finite() { *value } else { 0.0 })
            .collect();
        let (apex, height) = find_apex(&finite_y)?;
        if !(height > 0.0) || !x[apex].is_finite() {
            return None;
        }
        let data = DataXY {
            x: x.to_vec(),
            y: finite_y,
        };
        let seed = PeakSeed {
            rt: x[apex],
            intensity: height,
        };
        match fit_peak(&data, &seed, self.shape) {
            Some(params) if params.r2.is_finite() => Some(params.r2),
            _ => Some(0.0),
        }
    }

    pub fn keeps_peak(&self, x: &[f64], y: &[f64]) -> bool {
        if self.is_off() {
            return true;
        }
        if is_flat_topped(y) {
            return false;
        }
        match self.score(x, y) {
            Some(r2) => r2 >= self.min_r2,
            None => true,
        }
    }

    pub fn retain<C: Candidate>(
        &self,
        candidates: &mut Vec<C>,
        x: &[f64],
        y: &[f64],
        baseline: &[f64],
    ) {
        if self.is_off() {
            for candidate in candidates.iter_mut() {
                candidate.set_score(f64::NAN);
            }
            return;
        }
        candidates.retain_mut(|candidate| self.keep(candidate, x, y, baseline));
    }

    fn keep<C: Candidate>(
        &self,
        candidate: &mut C,
        x: &[f64],
        y: &[f64],
        baseline: &[f64],
    ) -> bool {
        let (start, end) = candidate.bounds();
        if start >= end || end >= y.len() || baseline.len() != y.len() {
            candidate.set_score(f64::NAN);
            return true;
        }
        let (from, to) = fixed_window(x, y, baseline, start, end);
        let peak = subtract_baseline(&y[from..=to], &baseline[from..=to]);
        match self.score(&x[from..=to], &peak) {
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

fn fixed_window(
    x: &[f64],
    y: &[f64],
    baseline: &[f64],
    start: usize,
    end: usize,
) -> (usize, usize) {
    let apex = apex_in_range(y, baseline, start, end);
    let Some(width) = full_width_at_half_max(x, y, baseline, apex) else {
        return (start, end);
    };
    if !(width > 0.0) {
        return (start, end);
    }
    let reach = FIT_HALF_WIDTH_FWHM * width;
    let from = x.partition_point(|&value| value < x[apex] - reach);
    let to = x
        .partition_point(|&value| value <= x[apex] + reach)
        .saturating_sub(1)
        .clamp(from, x.len() - 1);
    if to - from + 1 >= MIN_FIT_POINTS {
        (from, to)
    } else {
        (start, end)
    }
}

fn apex_in_range(y: &[f64], baseline: &[f64], start: usize, end: usize) -> usize {
    let mut apex = start;
    let mut height = y[start] - baseline[start];
    for index in start..=end {
        let value = y[index] - baseline[index];
        if value > height {
            height = value;
            apex = index;
        }
    }
    apex
}

fn full_width_at_half_max(x: &[f64], y: &[f64], baseline: &[f64], apex: usize) -> Option<f64> {
    let height = y[apex] - baseline[apex];
    if !(height > 0.0) {
        return None;
    }
    let half = 0.5 * height;
    let signal = |index: usize| y[index] - baseline[index];

    let mut left = apex;
    while left > 0 && signal(left) >= half {
        left -= 1;
    }
    let mut right = apex;
    while right + 1 < x.len() && signal(right) >= half {
        right += 1;
    }

    let crossed_left = signal(left) < half;
    let crossed_right = signal(right) < half;
    match (crossed_left, crossed_right) {
        (true, true) => Some(x[right] - x[left]),
        (true, false) => Some(2.0 * (x[apex] - x[left])),
        (false, true) => Some(2.0 * (x[right] - x[apex])),
        (false, false) => None,
    }
}

fn find_apex(y: &[f64]) -> Option<(usize, f64)> {
    let mut apex = 0;
    let mut height = f64::NEG_INFINITY;
    for (index, &value) in y.iter().enumerate() {
        if value > height {
            height = value;
            apex = index;
        }
    }
    height.is_finite().then_some((apex, height))
}

fn subtract_baseline(y: &[f64], baseline: &[f64]) -> Vec<f64> {
    y.iter()
        .zip(baseline.iter())
        .map(|(value, base)| (value - base).max(0.0))
        .collect()
}
