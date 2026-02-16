use core::f64;

use crate::utilities::cheminfo::{
    sgg::{SggOptions, sgg},
    {Point, kmeans},
};
use crate::utilities::closest_index;
use crate::utilities::structs::DataXY;

const DEFAULT_WINDOW_SIZES: &[usize] = &[5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27];
const EPSILON: f64 = 1e-5;

pub fn scan_for_peaks(data: &DataXY) -> Vec<f64> {
    let n = data.x.len();
    if n < 3 || n != data.y.len() {
        return Vec::new();
    }

    let mut window_sizes = Vec::with_capacity(DEFAULT_WINDOW_SIZES.len());
    window_sizes.extend(
        DEFAULT_WINDOW_SIZES
            .iter()
            .copied()
            .filter(|&window_size| window_size <= n),
    );
    if window_sizes.is_empty() {
        return Vec::new();
    }
    let mean_rt_step = mean_step(&data.x);
    let mut candidate_positions: Vec<f64> = Vec::new();
    let mut candidate_heights: Vec<f64> = Vec::new();
    for window_size in window_sizes.iter().copied() {
        for Peak { x, y } in scan_one_window(data, window_size, mean_rt_step, EPSILON) {
            candidate_positions.push(x);
            candidate_heights.push(y);
        }
    }
    if candidate_positions.is_empty() {
        return Vec::new();
    }

    let largest_window = *window_sizes.iter().max().unwrap() as f64;

    let mut grouping_tolerance = (8.0 * mean_rt_step).max(0.02);
    let window_tolerance_bump = 0.10 * largest_window * mean_rt_step;
    if window_tolerance_bump > grouping_tolerance {
        grouping_tolerance = window_tolerance_bump;
    }

    let groups = group_peak_positions(&candidate_positions, &candidate_heights, grouping_tolerance);
    if groups.is_empty() {
        return Vec::new();
    }

    let required_votes = if window_sizes.len() <= 3 {
        2
    } else {
        ((window_sizes.len() as f64) * 0.30).ceil() as usize
    };

    let mut global_best_height = 0.0f64;
    for g in &groups {
        if g.best_height > global_best_height {
            global_best_height = g.best_height;
        }
    }
    let strong_peak_gate = 0.60 * global_best_height;

    let mut kept_groups: Vec<(f64, f64)> = Vec::new();
    for g in &groups {
        if g.votes >= required_votes || g.best_height >= strong_peak_gate {
            kept_groups.push((g.center, g.best_height));
        }
    }
    if kept_groups.is_empty() {
        return Vec::new();
    }
    kept_groups.sort_by(|a, b| a.0.total_cmp(&b.0));

    let largest_window_size = *window_sizes.iter().max().unwrap();
    let signal_intensity: Vec<f64> = data.y.iter().map(|&v| v as f64).collect();

    let smoothed_intensity = if largest_window_size <= n {
        sgg(
            &signal_intensity,
            &data.x,
            SggOptions {
                window_size: largest_window_size,
                derivative: 0,
                polynomial: 3,
            },
        )
    } else {
        signal_intensity
    };

    let mut peak_groups: Vec<PeakGroup> = Vec::new();
    let mut group_index = 0usize;
    while group_index < kept_groups.len() {
        let (left_center, left_height) = kept_groups[group_index];
        if group_index + 1 >= kept_groups.len() {
            peak_groups.push(PeakGroup {
                center: left_center,
                height: left_height,
            });
            break;
        }
        let (right_center, right_height) = kept_groups[group_index + 1];

        let left_index = closest_index(&data.x, left_center);
        let right_index = closest_index(&data.x, right_center);
        let segment_left = left_index.min(right_index);
        let segment_right = left_index.max(right_index);

        let mut valley_value = f64::INFINITY;
        for j in segment_left..=segment_right {
            let v = smoothed_intensity[j];
            if v < valley_value {
                valley_value = v;
            }
        }

        let higher_peak_intensity = if left_height > right_height {
            left_height
        } else {
            right_height
        };

        let relative_valley_drop = if higher_peak_intensity > 0.0 {
            (higher_peak_intensity - valley_value) / higher_peak_intensity
        } else {
            0.0
        };

        if valley_value > higher_peak_intensity || relative_valley_drop < 0.08 {
            peak_groups.push(if left_height >= right_height {
                PeakGroup {
                    center: left_center,
                    height: left_height,
                }
            } else {
                PeakGroup {
                    center: right_center,
                    height: right_height,
                }
            });
            group_index += 2;
        } else {
            peak_groups.push(PeakGroup {
                center: left_center,
                height: left_height,
            });
            group_index += 1;
        }
    }
    if peak_groups.is_empty() {
        return Vec::new();
    }

    let mut peak_positions: Vec<f64> = Vec::new();
    let mut last_group = PeakGroup {
        center: f64::NEG_INFINITY,
        height: f64::NEG_INFINITY,
    };

    for PeakGroup {
        center: group_center_position,
        height: _group_height,
    } in peak_groups
    {
        let center_position_index = closest_index(&data.x, group_center_position);
        let apex_index = snap_to_local_max(&data.y, center_position_index);
        let apex_position = data.x[apex_index];
        let apex_height = data.y[apex_index] as f64;

        if !last_group.center.is_finite() || apex_position - last_group.center >= mean_rt_step {
            peak_positions.push(apex_position);
            last_group.center = apex_position;
            last_group.height = apex_height;
        } else if apex_height > last_group.height {
            if let Some(t) = peak_positions.last_mut() {
                *t = apex_position;
            }
            last_group.center = apex_position;
            last_group.height = apex_height;
        }
    }

    peak_positions
}

struct Peak {
    x: f64,
    y: f64,
}

struct PeakGroup {
    center: f64,
    height: f64,
}

fn scan_one_window(data: &DataXY, window_size: usize, mean_step: f64, epsilon: f64) -> Vec<Peak> {
    let DataXY { x, y } = data;
    let data_len = x.len();
    if data_len < 3 {
        return Vec::new();
    }

    let smooth_window = window_size.max(5);

    let smoothed_y = sgg(
        &y,
        &data.x,
        SggOptions {
            window_size: smooth_window,
            derivative: 0,
            polynomial: 3,
        },
    );
    let first_derivative = sgg(
        &y,
        &data.x,
        SggOptions {
            window_size: smooth_window,
            derivative: 1,
            polynomial: 3,
        },
    );
    let second_derivative = sgg(
        &y,
        &data.x,
        SggOptions {
            window_size: smooth_window,
            derivative: 2,
            polynomial: 3,
        },
    );

    let mut candidates: Vec<Peak> = Vec::with_capacity(data_len / 3);
    for i in 0..(data_len - 1) {
        if turns_down(first_derivative[i], first_derivative[i + 1], epsilon) {
            let idx = if smoothed_y[i] >= smoothed_y[i + 1] {
                i
            } else {
                i + 1
            };

            let apex_idx = snap_to_local_max(&smoothed_y, idx);
            if second_derivative[apex_idx] < 0.0 {
                candidates.push(Peak {
                    x: x[apex_idx],
                    y: smoothed_y[apex_idx],
                });
            }
        }
    }

    if candidates.is_empty() {
        return Vec::new();
    }

    candidates.sort_by(|a, b| a.x.total_cmp(&b.x));
    let min_candidate_separation = (1.2 * mean_step).max(0.01);

    let mut peaks: Vec<Peak> = Vec::with_capacity(candidates.len());
    let mut last_peak = Peak {
        x: f64::NEG_INFINITY,
        y: f64::NEG_INFINITY,
    };

    for Peak { x, y } in candidates {
        if !last_peak.x.is_finite() || x - last_peak.x >= min_candidate_separation {
            peaks.push(Peak { x, y });
            last_peak.x = x;
            last_peak.y = y;
        } else if y > last_peak.y {
            if let Some(peak) = peaks.last_mut() {
                peak.x = x;
                peak.y = y;
            }
            last_peak.x = x;
            last_peak.y = y;
        }
    }

    peaks
}

struct Group {
    center: f64,
    votes: usize,
    best_height: f64,
}

impl Default for Group {
    fn default() -> Self {
        Self {
            center: f64::NAN,
            votes: 0,
            best_height: f64::NEG_INFINITY,
        }
    }
}

fn group_peak_positions(positions: &[f64], heights: &[f64], tolerance: f64) -> Vec<Group> {
    let n = positions.len();
    if positions.is_empty() || n != heights.len() {
        return Vec::new();
    }
    let mut xmin = positions[0];
    let mut xmax = positions[0];
    for &position in positions {
        if position < xmin {
            xmin = position;
        }
        if position > xmax {
            xmax = position;
        }
    }

    let range: f64 = (xmax - xmin).abs();
    let number_of_clusters = if tolerance > 0.0 && range > 0.0 {
        (range / tolerance).ceil() as usize
    } else {
        1
    }
    .clamp(1, n);

    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|a, b| positions[*a].partial_cmp(&positions[*b]).unwrap());

    let mut initial_centroids = Vec::<Point>::with_capacity(number_of_clusters);
    let denom = number_of_clusters.saturating_sub(1);

    for c in 0..number_of_clusters {
        let ordered_index = if denom == 0 {
            order[n / 2]
        } else {
            let pos = c * (n - 1) / denom;
            order[pos]
        };

        initial_centroids.push(vec![positions[ordered_index]]);
    }

    let points: Vec<[f64; 1]> = positions.iter().map(|&x| [x]).collect();
    let centroids = kmeans(&points, initial_centroids);
    if centroids.is_empty() {
        return Vec::new();
    }

    let mut groups: Vec<Group> = Vec::with_capacity(centroids.len());
    for centroid in &centroids {
        groups.push(Group {
            center: centroid[0],
            votes: 0,
            best_height: f64::NEG_INFINITY,
        });
    }

    for i in 0..n {
        let position = positions[i];
        let height = heights[i];

        let mut idx = 0usize;
        let mut best = (position - groups[0].center).abs();
        for j in 1..groups.len() {
            let d = (position - groups[j].center).abs();
            if d < best {
                best = d;
                idx = j;
            }
        }

        groups[idx].votes += 1;
        if height > groups[idx].best_height {
            groups[idx].best_height = height;
        }
    }

    groups
}

fn snap_to_local_max(y: &[f64], mut idx: usize) -> usize {
    let n = y.len();
    if n == 0 {
        return idx;
    }
    if idx >= n {
        idx = n - 1;
    }

    loop {
        let mut best = idx;
        let mut best_y = y[idx];

        if idx > 0 && y[idx - 1] > best_y {
            best = idx - 1;
            best_y = y[idx - 1];
        }
        if idx + 1 < n && y[idx + 1] > best_y {
            best = idx + 1;
        }

        if best == idx {
            break;
        }
        idx = best;
    }
    idx
}

#[inline]
fn turns_down(left: f64, right: f64, epsilon: f64) -> bool {
    let left_up = left > epsilon;
    let left_down = left < -epsilon;
    let right_up = right > epsilon;
    let right_down = right < -epsilon;

    (left_up && !right_up) || (!left_down && right_down)
}

#[inline]
pub fn quad_peak(xs: &[f64], ys: &[f64], i: usize) -> f64 {
    let xm1 = xs[i - 1];
    let x0 = xs[i];
    let xp1 = xs[i + 1];
    let ym1 = ys[i - 1] as f64;
    let y0 = ys[i] as f64;
    let yp1 = ys[i + 1] as f64;

    let a0 = xm1 - x0;
    let a1 = xp1 - x0;
    let dy0 = ym1 - y0;
    let dy1 = yp1 - y0;
    let denom = a0 * a1 * (a0 - a1);
    if denom == 0.0 {
        return x0;
    }
    let a = (dy0 * a1 - dy1 * a0) / denom;
    if a == 0.0 {
        x0
    } else {
        let b = (dy1 * a0 * a0 - dy0 * a1 * a1) / denom;
        x0 - b / (2.0 * a)
    }
}

#[inline]
fn mean_step(xs: &[f64]) -> f64 {
    let mut sum = 0.0;
    let mut c = 0usize;
    for w in xs.windows(2) {
        let d = (w[1] - w[0]).abs();
        if d.is_finite() && d > 0.0 {
            sum += d;
            c += 1;
        }
    }
    if c == 0 {
        f64::EPSILON.max(0.01)
    } else {
        (sum / c as f64).max(f64::EPSILON)
    }
}
