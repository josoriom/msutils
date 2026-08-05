use std::{borrow::Cow, f64};

use crate::utilities::cheminfo::{
    kmeans,
    noise_san_plot::{NoiseSanPlotOptions, noise_san_plot},
};

const SPREAD_TO_DEVIATION: f64 = 1.4826;
const DEVIATIONS_ALLOWED: f64 = 12.0;

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
    let floor = min_value(&y);
    let segments = scan_segments(y.as_ref(), floor);
    if segments.is_empty() {
        return Noise::default();
    }
    if segments.len() == 1 {
        return Noise {
            width: 0,
            intensity: floor,
        };
    }

    let params = ClusterParams::default();
    let labels = cluster_intensities(&segments, params);

    let noise_label = 0usize;

    let mut noise_levels = Vec::with_capacity(segments.len());
    let mut noise_width = 0usize;
    for (i, &Segment { width, intensity }) in segments.iter().enumerate() {
        if labels[i] == noise_label && intensity.is_finite() {
            noise_levels.push(intensity);
            if width > noise_width {
                noise_width = width;
            }
        }
    }

    match noise_ceiling(&mut noise_levels) {
        Some(intensity) => Noise {
            width: noise_width,
            intensity,
        },
        None => Noise::default(),
    }
}

fn noise_ceiling(levels: &mut [f64]) -> Option<f64> {
    let middle = median(levels)?;
    let mut spreads: Vec<f64> = levels.iter().map(|level| (level - middle).abs()).collect();
    let spread = median(&mut spreads).unwrap_or(0.0);
    let ceiling = middle + spread * SPREAD_TO_DEVIATION * DEVIATIONS_ALLOWED;
    let highest = levels.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    Some(ceiling.min(highest))
}

fn median(values: &mut [f64]) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    values.sort_unstable_by(f64::total_cmp);
    let middle = values.len() / 2;
    if values.len() % 2 == 1 {
        Some(values[middle])
    } else {
        Some((values[middle - 1] + values[middle]) / 2.0)
    }
}

pub fn find_noise_level_san_plot<'a, Y>(y: Y) -> Noise
where
    Y: IntoF64Slice<'a>,
{
    let y = y.into_f64_slice();
    if y.is_empty() {
        return Noise::default();
    }

    let result = noise_san_plot(
        y.as_ref(),
        NoiseSanPlotOptions {
            cut_off: Some(0.5),
            refine: false,
            fix_offset: false,
            ..Default::default()
        },
    );
    let intensity = if result.positive.is_finite() && result.positive > 0.0 {
        result.positive
    } else {
        0.0
    };

    Noise {
        width: 0,
        intensity,
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

    let n = points.len() as f64;
    let mut w_sum = 0.0;
    let mut i_sum = 0.0;
    for seg in points {
        w_sum += log1p_safe(seg.width as f64);
        i_sum += log1p_safe(seg.intensity);
    }
    let m_w = w_sum / n;
    let m_i = i_sum / n;

    let mut w_var = 0.0;
    let mut i_var = 0.0;
    for seg in points {
        let dw = log1p_safe(seg.width as f64) - m_w;
        let di = log1p_safe(seg.intensity) - m_i;
        w_var += dw * dw;
        i_var += di * di;
    }
    let s_w = (w_var / n).sqrt().max(1e-6);
    let s_i = (i_var / n).sqrt().max(1e-6);

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

fn min_value(array: &[f64]) -> f64 {
    let mut minimum = f64::INFINITY;
    for &value in array {
        if value.is_finite() && value < minimum {
            minimum = value;
        }
    }
    if minimum.is_finite() { minimum } else { 0.0 }
}

fn log1p_safe(v: f64) -> f64 {
    (0.0_f64.max(v)).ln_1p()
}
