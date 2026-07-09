use std::cmp::Ordering;

use crate::utilities::cheminfo::sgg::{SggOptions, sgg};
use crate::utilities::{closest_index, structs::DataXY};

#[derive(Clone, Copy, Debug)]
pub struct Boundary {
    pub index: Option<usize>,
    pub value: Option<f64>,
}

#[derive(Clone, Copy, Debug)]
pub struct Boundaries {
    pub from: Boundary,
    pub to: Boundary,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum BoundaryMethod {
    #[default]
    Walk,
    Derivative,
}

#[derive(Clone, Copy, Debug)]
pub struct BoundariesOptions {
    pub method: BoundaryMethod,
    pub min_slope_step: f64,
    pub noise: f64,
    pub min_ascending_steps: usize,
    pub min_below_noise_run: usize,
    pub smooth_window: usize,
    pub smooth_polynomial: usize,
}

impl Default for BoundariesOptions {
    fn default() -> Self {
        Self {
            method: BoundaryMethod::default(),
            min_slope_step: 1e-5,
            noise: 0.0,
            min_ascending_steps: 3,
            min_below_noise_run: 2,
            smooth_window: 7,
            smooth_polynomial: 3,
        }
    }
}

const DESCEND_FRACTION: f64 = 0.02;
const RISE_FRACTION: f64 = 0.07;
const RISE_NOISE_MULTIPLE: f64 = 3.0;
const RECENTER_WINDOW: usize = 4;
const CONFIRM_STEPS: usize = 5;
const BASELINE_PERCENTILE: f64 = 0.10;

pub fn get_boundaries(data: &DataXY, peak_x: f64, options: Option<BoundariesOptions>) -> Boundaries {
    let n = data.x.len();
    if n < 2 || n != data.y.len() {
        return Boundaries {
            from: Boundary {
                index: None,
                value: None,
            },
            to: Boundary {
                index: None,
                value: None,
            },
        };
    }

    let opts = options.unwrap_or_default();
    let apex_index = closest_index(&data.x, peak_x);
    let lowest_value = min_value(&data.y);

    match opts.method {
        BoundaryMethod::Walk => boundaries_walk(data, apex_index, &opts),
        BoundaryMethod::Derivative => boundaries_derivative(data, apex_index, &opts, lowest_value),
    }
}

fn boundaries_walk(data: &DataXY, apex_index: usize, opts: &BoundariesOptions) -> Boundaries {
    let smoothed = smooth(&data.y, &data.x, opts);
    let apex = recenter_apex(&smoothed, apex_index);
    let baseline = baseline_of(&smoothed);
    let span = (smoothed[apex] - baseline).max(opts.min_slope_step);
    let floor = baseline + DESCEND_FRACTION * span;
    let rise = (RISE_FRACTION * span).max(RISE_NOISE_MULTIPLE * opts.noise);

    Boundaries {
        from: descend(&data.x, &smoothed, apex, -1, floor, rise),
        to: descend(&data.x, &smoothed, apex, 1, floor, rise),
    }
}

fn smooth(y: &[f64], x: &[f64], opts: &BoundariesOptions) -> Vec<f64> {
    match usable_window(opts.smooth_window, y.len()) {
        Some(window) => sgg(
            y,
            x,
            SggOptions {
                window_size: window,
                derivative: 0,
                polynomial: opts.smooth_polynomial,
            },
        ),
        None => y.to_vec(),
    }
}

fn recenter_apex(smoothed: &[f64], seed: usize) -> usize {
    let start = seed.saturating_sub(RECENTER_WINDOW);
    let end = (seed + RECENTER_WINDOW + 1).min(smoothed.len());
    let mut apex = seed;
    for index in start..end {
        if smoothed[index] > smoothed[apex] {
            apex = index;
        }
    }
    apex
}

fn min_value(values: &[f64]) -> f64 {
    values.iter().copied().fold(f64::INFINITY, f64::min)
}

fn baseline_of(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted = values.to_vec();
    let index = ((sorted.len() - 1) as f64 * BASELINE_PERCENTILE) as usize;
    sorted.select_nth_unstable_by(index, |a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    sorted[index]
}

fn descend(
    x: &[f64],
    smoothed: &[f64],
    apex: usize,
    direction: isize,
    floor: f64,
    rise: f64,
) -> Boundary {
    let length = smoothed.len() as isize;
    let mut lowest = smoothed[apex];
    let mut lowest_index = apex;
    let mut current = apex as isize;

    while current + direction >= 0 && current + direction < length {
        let next = (current + direction) as usize;
        let value = smoothed[next];
        if value <= floor {
            return boundary_at(x, next);
        }
        if value < lowest {
            lowest = value;
            lowest_index = next;
        } else if value - lowest >= rise && stays_up(smoothed, next, direction, lowest + rise) {
            return boundary_at(x, lowest_index);
        }
        current += direction;
    }

    let edge = if direction > 0 { smoothed.len() - 1 } else { 0 };
    boundary_at(x, edge)
}

fn stays_up(smoothed: &[f64], from: usize, direction: isize, threshold: f64) -> bool {
    let length = smoothed.len() as isize;
    let mut index = from as isize;
    for _ in 0..CONFIRM_STEPS {
        if index < 0 || index >= length || smoothed[index as usize] < threshold {
            return false;
        }
        index += direction;
    }
    true
}

fn usable_window(requested: usize, length: usize) -> Option<usize> {
    if length < 5 {
        return None;
    }
    let mut window = requested.min(length);
    if window.is_multiple_of(2) {
        window -= 1;
    }
    (window >= 5).then_some(window)
}

fn boundaries_derivative(
    data: &DataXY,
    apex_index: usize,
    opts: &BoundariesOptions,
    lowest_value: f64,
) -> Boundaries {
    let Some(window) = usable_window(opts.smooth_window, data.y.len()) else {
        return boundaries_walk(data, apex_index, opts);
    };

    let smoothed = sgg(
        &data.y,
        &data.x,
        SggOptions {
            window_size: window,
            derivative: 0,
            polynomial: opts.smooth_polynomial,
        },
    );
    let curvature = sgg(
        &data.y,
        &data.x,
        SggOptions {
            window_size: window,
            derivative: 2,
            polynomial: opts.smooth_polynomial,
        },
    );

    let noise = (lowest_value + opts.min_slope_step).max(opts.noise);

    let from = find_edge(&data.x, &smoothed, &curvature, apex_index, -1, noise);
    let to = find_edge(&data.x, &smoothed, &curvature, apex_index, 1, noise);

    Boundaries { from, to }
}

fn find_edge(
    x: &[f64],
    smoothed: &[f64],
    curvature: &[f64],
    apex_index: usize,
    direction: isize,
    noise: f64,
) -> Boundary {
    let length = smoothed.len() as isize;
    let mut current = apex_index as isize;
    let mut passed_inflection = false;

    while current + direction >= 0 && current + direction < length {
        let next = current + direction;
        let current_index = current as usize;
        let next_index = next as usize;

        if !passed_inflection && curvature[next_index] > 0.0 {
            passed_inflection = true;
        }

        if smoothed[next_index] <= noise {
            return boundary_at(x, next_index);
        }

        if passed_inflection && smoothed[next_index] >= smoothed[current_index] {
            return boundary_at(x, current_index);
        }

        current = next;
    }

    let edge = if direction > 0 { (length - 1) as usize } else { 0 };
    boundary_at(x, edge)
}

fn boundary_at(x: &[f64], index: usize) -> Boundary {
    Boundary {
        index: Some(index),
        value: Some(x[index]),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn grid(from: f64, to: f64, points: usize) -> Vec<f64> {
        (0..points)
            .map(|index| from + (to - from) * index as f64 / (points - 1) as f64)
            .collect()
    }

    fn bell(x: f64, center: f64, height: f64, fwhm: f64) -> f64 {
        let sigma = fwhm / 2.354_820_045;
        height * (-0.5 * ((x - center) / sigma).powi(2)).exp()
    }

    fn wobble(index: usize) -> f64 {
        let value = (index as f64 * 12.9898).sin() * 43758.547;
        2.0 * (value - value.floor()) - 1.0
    }

    fn derivative_options() -> BoundariesOptions {
        BoundariesOptions {
            method: BoundaryMethod::Derivative,
            ..Default::default()
        }
    }

    #[test]
    fn derivative_brackets_the_apex() {
        let x = grid(4.0, 6.0, 400);
        let y: Vec<f64> = x.iter().map(|&v| 0.05 + bell(v, 5.0, 1.0, 0.2)).collect();
        let data = DataXY { x, y };

        let edges = get_boundaries(&data, 5.0, Some(derivative_options()));
        let from = edges.from.value.unwrap();
        let to = edges.to.value.unwrap();

        assert!(from < 5.0 && to > 5.0, "bracket [{from}, {to}] must contain the apex");
        assert!(from < 4.9 && to > 5.1, "bracket [{from}, {to}] must reach past the flanks");
    }

    #[test]
    fn derivative_stops_at_the_valley_between_two_peaks() {
        let x = grid(4.0, 6.5, 500);
        let y: Vec<f64> = x
            .iter()
            .map(|&v| 0.05 + bell(v, 5.0, 1.0, 0.2) + bell(v, 5.35, 0.8, 0.2))
            .collect();
        let data = DataXY { x, y };

        let to = get_boundaries(&data, 5.0, Some(derivative_options()))
            .to
            .value
            .unwrap();

        assert!(to > 5.0 && to < 5.35, "right edge {to} must land in the valley, not the next peak");
    }

    #[test]
    fn derivative_window_is_stable_under_noise() {
        let x = grid(4.0, 6.5, 500);
        let clean: Vec<f64> = x
            .iter()
            .map(|&v| 0.05 + bell(v, 5.0, 1.0, 0.2) + bell(v, 5.35, 0.8, 0.2))
            .collect();
        let noisy: Vec<f64> = clean
            .iter()
            .enumerate()
            .map(|(index, &value)| value + 0.02 * wobble(index))
            .collect();

        let clean_to = get_boundaries(&DataXY { x: x.clone(), y: clean }, 5.0, Some(derivative_options()))
            .to
            .value
            .unwrap();
        let noisy_to = get_boundaries(&DataXY { x, y: noisy }, 5.0, Some(derivative_options()))
            .to
            .value
            .unwrap();

        assert!((clean_to - noisy_to).abs() <= 0.05, "right edge moved {clean_to} -> {noisy_to} under noise");
    }

    #[test]
    fn derivative_keeps_left_edge_when_right_side_is_truncated() {
        let full_x = grid(4.0, 6.0, 400);
        let full_y: Vec<f64> = full_x.iter().map(|&v| 0.05 + bell(v, 5.0, 1.0, 0.2)).collect();
        let full_from = get_boundaries(&DataXY { x: full_x.clone(), y: full_y.clone() }, 5.0, Some(derivative_options()))
            .from
            .value
            .unwrap();

        let cut = full_x.partition_point(|&v| v <= 5.25);
        let cut_x = full_x[..cut].to_vec();
        let cut_y = full_y[..cut].to_vec();
        let cut_edges = get_boundaries(&DataXY { x: cut_x, y: cut_y }, 5.0, Some(derivative_options()));
        let cut_from = cut_edges.from.value.unwrap();
        let cut_to = cut_edges.to.value.unwrap();

        assert!((full_from - cut_from).abs() <= 0.02, "left edge moved {full_from} -> {cut_from} after truncation");
        assert!(cut_to >= 5.2, "right edge {cut_to} should follow the data to its truncated end");
    }

    #[test]
    fn walk_method_brackets_the_apex() {
        let x = grid(4.0, 6.0, 400);
        let y: Vec<f64> = x.iter().map(|&v| 0.05 + bell(v, 5.0, 1.0, 0.2)).collect();
        let options = BoundariesOptions {
            method: BoundaryMethod::Walk,
            ..Default::default()
        };

        let edges = get_boundaries(&DataXY { x, y }, 5.0, Some(options));

        assert!(edges.from.value.unwrap() < 5.0 && edges.to.value.unwrap() > 5.0);
    }

    #[test]
    fn short_input_does_not_panic() {
        let x = grid(0.0, 1.0, 4);
        let y = vec![0.0, 1.0, 0.5, 0.2];
        let edges = get_boundaries(&DataXY { x, y }, 0.3, Some(derivative_options()));
        assert!(edges.from.index.is_some() && edges.to.index.is_some());
    }
}
