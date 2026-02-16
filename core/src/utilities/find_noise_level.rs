use std::borrow::Cow;
use std::f64;

use crate::utilities::cheminfo::kmeans;

pub fn find_noise_level<'a, Y>(y: Y) -> f64
where
    Y: IntoF64Slice<'a>,
{
    let y = y.into_f64_slice();
    let segments = scan_segments(y.as_ref());
    if segments.is_empty() {
        return 0.0;
    }

    let params = ClusterParams::default();
    let labels = cluster_intensities(&segments, params);

    let mut counts = [0usize; 2];
    for &label in &labels {
        counts[label] += 1;
    }
    let noise_label = if counts[0] >= counts[1] { 0 } else { 1 };

    let mut max_in_noise = f64::NEG_INFINITY;
    for (
        i,
        &Segment {
            width: _,
            intensity,
        },
    ) in segments.iter().enumerate()
    {
        if labels[i] == noise_label && intensity.is_finite() && intensity > max_in_noise {
            max_in_noise = intensity;
        }
    }

    if max_in_noise.is_finite() {
        max_in_noise
    } else {
        0.0
    }
}

pub trait IntoF64Slice<'a> {
    fn into_f64_slice(self) -> Cow<'a, [f64]>;
}

impl<'a> IntoF64Slice<'a> for &'a [f64] {
    fn into_f64_slice(self) -> Cow<'a, [f64]> {
        Cow::Borrowed(self)
    }
}
impl<'a> IntoF64Slice<'a> for &'a Vec<f64> {
    fn into_f64_slice(self) -> Cow<'a, [f64]> {
        Cow::Borrowed(self.as_slice())
    }
}
impl<'a> IntoF64Slice<'a> for Vec<f64> {
    fn into_f64_slice(self) -> Cow<'a, [f64]> {
        Cow::Owned(self)
    }
}

impl<'a> IntoF64Slice<'a> for &'a [f32] {
    fn into_f64_slice(self) -> Cow<'a, [f64]> {
        let mut v = Vec::with_capacity(self.len());
        for &x in self {
            v.push(x as f64);
        }
        Cow::Owned(v)
    }
}
impl<'a> IntoF64Slice<'a> for &'a Vec<f32> {
    fn into_f64_slice(self) -> Cow<'a, [f64]> {
        let mut v = Vec::with_capacity(self.len());
        for &x in self {
            v.push(x as f64);
        }
        Cow::Owned(v)
    }
}
impl<'a> IntoF64Slice<'a> for Vec<f32> {
    fn into_f64_slice(self) -> Cow<'a, [f64]> {
        let mut v = Vec::with_capacity(self.len());
        for x in self {
            v.push(x as f64);
        }
        Cow::Owned(v)
    }
}

#[derive(Clone, Copy)]
struct ClusterParams {
    width_weight: f64,
    intensity_weight: f64,
}

impl Default for ClusterParams {
    fn default() -> Self {
        Self {
            width_weight: 2.0,
            intensity_weight: 1.0,
        }
    }
}

struct Segment {
    width: f64,
    intensity: f64,
}

fn scan_segments(y: &[f64]) -> Vec<Segment> {
    let mut segments: Vec<Segment> = Vec::new();
    if y.is_empty() {
        return segments;
    }

    let threshold = x_min_value(y);

    let mut in_segment = false;
    let mut segment_width: usize = 0;
    let mut segment_intensity = f64::NEG_INFINITY;

    for &intensity in y {
        let mask_on = intensity.is_finite() && intensity.abs() > threshold;

        if mask_on {
            if !in_segment {
                in_segment = true;
                segment_width = 1;
                segment_intensity = intensity;
            } else {
                segment_width += 1;
                if intensity > segment_intensity {
                    segment_intensity = intensity;
                }
            }
        } else if in_segment {
            segments.push(Segment {
                width: segment_width as f64,
                intensity: segment_intensity,
            });
            in_segment = false;
            segment_width = 0;
            segment_intensity = f64::NEG_INFINITY;
        }
    }

    if in_segment {
        segments.push(Segment {
            width: segment_width as f64,
            intensity: segment_intensity,
        });
    }

    segments
}

fn cluster_intensities(points: &[Segment], params: ClusterParams) -> Vec<usize> {
    if points.is_empty() {
        return Vec::new();
    }

    let mut widths: Vec<f64> = Vec::with_capacity(points.len());
    let mut intensities: Vec<f64> = Vec::with_capacity(points.len());
    for &Segment { width, intensity } in points.iter() {
        widths.push(log1p_safe(width));
        intensities.push(log1p_safe(intensity));
    }

    let meam_widths = mean(&widths);
    let mean_intensity = mean(&intensities);

    let width_std = match std_dev(&widths) {
        0.0 => 1.0,
        result => result,
    };
    let intensity_std = match std_dev(&intensities) {
        0.0 => 1.0,
        result => result,
    };

    let ClusterParams {
        width_weight,
        intensity_weight,
    } = params;

    let mut weighted_points: Vec<[f64; 2]> = Vec::with_capacity(points.len());
    for i in 0..points.len() {
        let a = (widths[i] - meam_widths) / width_std;
        let b = (intensities[i] - mean_intensity) / intensity_std;
        weighted_points.push([a * width_weight, b * intensity_weight]);
    }

    let seeds: Vec<Vec<f64>> = farthest_seeds(&weighted_points);
    let mut centroids = kmeans(&weighted_points, seeds);
    centroids.sort_by(|a, b| a[1].total_cmp(&b[1]));

    let mut result: Vec<usize> = Vec::with_capacity(weighted_points.len());
    for p in weighted_points.iter() {
        let mut idx = 0usize;
        let mut best = {
            let dx = p[0] - centroids[0][0];
            let dy = p[1] - centroids[0][1];
            dx * dx + dy * dy
        };
        for j in 1..centroids.len() {
            let dx = p[0] - centroids[j][0];
            let dy = p[1] - centroids[j][1];
            let d = dx * dx + dy * dy;
            if d < best {
                best = d;
                idx = j;
            }
        }
        result.push(idx);
    }

    result
}

fn farthest_seeds(weighted_points: &[[f64; 2]]) -> Vec<Vec<f64>> {
    let mut min = weighted_points[0];
    let mut max = weighted_points[0];
    for &point in weighted_points.iter() {
        if point[1] < min[1] {
            min = point;
        }
        if point[1] > max[1] {
            max = point;
        }
    }
    vec![vec![min[0], min[1]], vec![max[0], max[1]]]
}

fn x_min_value(array: &[f64]) -> f64 {
    let mut minimum = f64::INFINITY;
    for &value in array {
        if value.is_finite() {
            let absolute_value = value.abs();
            if absolute_value < minimum {
                minimum = absolute_value;
            }
        }
    }
    if minimum.is_finite() { minimum } else { 0.0 }
}

fn mean(array: &[f64]) -> f64 {
    if array.is_empty() {
        return 0.0;
    }
    let mut sum = 0.0f64;
    let mut n = 0usize;
    for &value in array {
        if value.is_finite() {
            sum += value;
            n += 1;
        }
    }
    if n == 0 { 0.0 } else { sum / n as f64 }
}

fn std_dev(array: &[f64]) -> f64 {
    if array.is_empty() {
        return 0.0;
    }
    let mean_value = mean(array);
    let mut sum_sq_diff = 0.0f64;
    let mut n = 0usize;
    for &value in array {
        if value.is_finite() {
            let diff = value - mean_value;
            sum_sq_diff += diff * diff;
            n += 1;
        }
    }
    if n == 0 {
        0.0
    } else {
        (sum_sq_diff / n as f64).sqrt()
    }
}

fn log1p_safe(v: f64) -> f64 {
    (0.0_f64.max(v)).ln_1p()
}
