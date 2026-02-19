use std::borrow::Cow;
use std::f64;

use crate::utilities::cheminfo::kmeans;

#[derive(Clone, Copy, Debug)]
pub struct Noise {
    pub width: usize,
    pub intensity: f64,
}

impl Default for Noise {
    fn default() -> Self {
        Self {
            width: 0,
            intensity: 0.0,
        }
    }
}

pub fn find_noise_level<'a, Y>(y: Y) -> Noise
where
    Y: IntoF64Slice<'a>,
{
    let y = y.into_f64_slice();
    let min_y = x_min_value(&y);
    let segments = scan_segments(y.as_ref(), min_y);
    if segments.is_empty() {
        return Noise::default();
    }

    let params = ClusterParams::default();
    let labels = cluster_intensities(&segments, params);

    let mut counts = [0usize; 2];
    for &label in &labels {
        counts[label] += 1;
    }
    let noise_label = if counts[0] >= counts[1] { 0 } else { 1 };

    let mut max_in_noise = f64::NEG_INFINITY;
    let mut max_noise_width = 0usize;
    for (i, &Segment { width, intensity }) in segments.iter().enumerate() {
        if labels[i] == noise_label {
            if intensity.is_finite() && intensity > max_in_noise {
                max_in_noise = intensity;
            }
            if width > max_noise_width {
                max_noise_width = width;
            }
        }
    }

    if max_in_noise.is_finite() {
        Noise {
            width: max_noise_width,
            intensity: max_in_noise,
        }
    } else {
        Noise::default()
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
    width: usize,
    intensity: f64,
}

fn scan_segments(y: &[f64], threshold: f64) -> Vec<Segment> {
    if y.is_empty() {
        return Vec::new();
    }

    let mut segments = Vec::with_capacity(32);

    let mut segment_width = 0usize;
    let mut segment_intensity = f64::NEG_INFINITY;

    for &val in y {
        if val > threshold && val.is_finite() {
            segment_width += 1;
            if val > segment_intensity {
                segment_intensity = val;
            }
        } else if segment_width > 0 {
            segments.push(Segment {
                width: segment_width,
                intensity: segment_intensity,
            });
            segment_width = 0;
            segment_intensity = f64::NEG_INFINITY;
        }
    }

    if segment_width > 0 {
        segments.push(Segment {
            width: segment_width,
            intensity: segment_intensity,
        });
    }
    segments
}

fn cluster_intensities(points: &[Segment], params: ClusterParams) -> Vec<usize> {
    if points.is_empty() {
        return Vec::new();
    }

    let mut w_sum = 0.0;
    let mut i_sum = 0.0;
    let mut w_sq_sum = 0.0;
    let mut i_sq_sum = 0.0;
    let n = points.len() as f64;

    for seg in points {
        let lw = log1p_safe(seg.width as f64);
        let li = log1p_safe(seg.intensity);
        w_sum += lw;
        i_sum += li;
        w_sq_sum += lw * lw;
        i_sq_sum += li * li;
    }

    let m_w = w_sum / n;
    let m_i = i_sum / n;

    // Variance = E[X^2] - (E[X])^2
    let s_w = ((w_sq_sum / n) - (m_w * m_w)).sqrt().max(1e-6);
    let s_i = ((i_sq_sum / n) - (m_i * m_i)).sqrt().max(1e-6);

    let weighted_points: Vec<[f64; 2]> = points
        .iter()
        .map(|seg| {
            let a = (log1p_safe(seg.width as f64) - m_w) / s_w;
            let b = (log1p_safe(seg.intensity) - m_i) / s_i;
            [a * params.width_weight, b * params.intensity_weight]
        })
        .collect();

    let seeds = farthest_seeds(&weighted_points);
    let mut centroids = kmeans(&weighted_points, seeds);
    centroids.sort_by(|a, b| a[1].total_cmp(&b[1]));

    weighted_points
        .iter()
        .map(|p| {
            let d0 = (p[0] - centroids[0][0]).powi(2) + (p[1] - centroids[0][1]).powi(2);
            let d1 = (p[0] - centroids[1][0]).powi(2) + (p[1] - centroids[1][1]).powi(2);
            if d0 < d1 { 0 } else { 1 }
        })
        .collect()
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

fn log1p_safe(v: f64) -> f64 {
    (0.0_f64.max(v)).ln_1p()
}
