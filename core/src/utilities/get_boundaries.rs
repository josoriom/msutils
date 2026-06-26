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

pub fn get_boundaries(
    data: &DataXY,
    peak_x: f64,
    options: Option<BoundariesOptions>,
) -> Boundaries {
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

    let mut lowest_value = f64::INFINITY;
    for &value in &data.y {
        if value < lowest_value {
            lowest_value = value;
        }
    }

    match opts.method {
        BoundaryMethod::Walk => boundaries_walk(data, apex_index, &opts, lowest_value),
        BoundaryMethod::Derivative => boundaries_derivative(data, apex_index, &opts, lowest_value),
    }
}

fn boundaries_walk(
    data: &DataXY,
    apex_index: usize,
    opts: &BoundariesOptions,
    lowest_value: f64,
) -> Boundaries {
    let config = WalkConfig {
        epsilon: opts.min_slope_step,
        noise_value: opts.noise,
        n_steps: opts.min_ascending_steps,
        baseline_run: opts.min_below_noise_run,
        global_min: lowest_value,
    };

    let from = walk(&data.x, &data.y, apex_index, -1, config);
    let to = walk(&data.x, &data.y, apex_index, 1, config);

    Boundaries { from, to }
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
        return boundaries_walk(data, apex_index, opts, lowest_value);
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

#[derive(Clone, Copy)]
struct WalkConfig {
    epsilon: f64,
    noise_value: f64,
    n_steps: usize,
    baseline_run: usize,
    global_min: f64,
}

fn walk(x: &[f64], y: &[f64], start: usize, direction: isize, config: WalkConfig) -> Boundary {
    let n = x.len() as isize;
    if n < 2 {
        return Boundary {
            index: None,
            value: None,
        };
    }

    let dir = if direction > 0 { 1 } else { -1 };
    let noise = (config.global_min + config.epsilon).max(config.noise_value);

    let mut i = start as isize;

    let mut checking = false;
    let mut steps_up: usize = 0;
    let mut has_risen: bool = false;
    let mut valley_idx: isize = start as isize;
    let mut valley_val: f64 = y[start];
    let mut below_noise: bool = false;

    let mut below_noise_run_len: usize = 0;
    let mut below_noise_start: isize = 0;

    while i >= 0 && i + dir >= 0 && i + dir < n {
        let j = i + dir;
        let iu = i as usize;
        let ju = j as usize;

        let y_i = y[iu];
        let y_j = y[ju];

        if running_below_noise(
            y_j,
            noise,
            config.baseline_run,
            j,
            &mut below_noise_run_len,
            &mut below_noise_start,
        ) {
            let k = below_noise_start.clamp(0, n - 1) as usize;
            return Boundary {
                index: Some(k),
                value: Some(x[k]),
            };
        }

        let slope = compute_ratio_and_slope(x, y, iu, ju, dir, config.epsilon);

        if is_asc_or_flat(slope) {
            if !checking {
                checking = true;
                steps_up = 1;
                valley_idx = i;
                valley_val = y_i;
                below_noise = y_j <= noise;
            } else {
                steps_up += 1;
                if y_j <= noise {
                    below_noise = true;
                }
            }

            if steps_up >= config.n_steps {
                let end_val = y_j;
                let rise = end_val - valley_val;
                if !allow_rise(below_noise, end_val, noise, rise) {
                    let k = valley_idx.clamp(0, n - 1) as usize;
                    return Boundary {
                        index: Some(k),
                        value: Some(x[k]),
                    };
                }

                reset_state(
                    &mut checking,
                    &mut steps_up,
                    &mut has_risen,
                    &mut below_noise,
                );
            }
        } else {
            reset_state(
                &mut checking,
                &mut steps_up,
                &mut has_risen,
                &mut below_noise,
            );
        }

        i = j;
    }

    let edge = if dir > 0 { (n - 1) as usize } else { 0usize };
    Boundary {
        index: Some(edge),
        value: Some(x[edge]),
    }
}

fn running_below_noise(
    y_j: f64,
    noise: f64,
    baseline_run: usize,
    j: isize,
    below_noise_run_len: &mut usize,
    below_noise_start: &mut isize,
) -> bool {
    if y_j <= noise {
        if *below_noise_run_len == 0 {
            *below_noise_start = j;
        }
        *below_noise_run_len += 1;
        if *below_noise_run_len >= baseline_run {
            return true;
        }
    } else {
        *below_noise_run_len = 0;
    }
    false
}

fn reset_state(
    checking: &mut bool,
    steps_up: &mut usize,
    has_risen: &mut bool,
    below_noise: &mut bool,
) {
    *checking = false;
    *steps_up = 0;
    *has_risen = false;
    *below_noise = false;
}

fn is_asc_or_flat(slope: f64) -> bool {
    slope >= 0.0
}

#[inline]
fn allow_rise(below_noise: bool, end_val: f64, noise: f64, rise: f64) -> bool {
    !((below_noise && end_val <= noise) || (rise >= noise))
}

fn compute_ratio_and_slope(
    x: &[f64],
    y: &[f64],
    iu: usize,
    ju: usize,
    dir: isize,
    epsilon: f64,
) -> f64 {
    let dx = x[ju] - x[iu];
    let denom = if dx.abs() < epsilon {
        if dir > 0 { epsilon } else { -epsilon }
    } else {
        dx
    };
    let dy = (y[ju]) - (y[iu]);
    (dy / denom) * if dir > 0 { 1.0 } else { -1.0 }
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
