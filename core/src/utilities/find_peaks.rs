use std::cmp::Ordering;

use crate::utilities::calculate_baseline::{BaselineOptions, calculate_baseline};
use crate::utilities::cheminfo::sgg::{SggOptions, sgg};
use crate::utilities::closest_index;
use crate::utilities::find_noise_level::{find_noise_level, find_noise_level_san_plot};
use crate::utilities::fit_peak::PeakShape;
use crate::utilities::get_boundaries::{Boundaries, BoundariesOptions, get_boundaries};
use crate::utilities::math::xy_integration;
use crate::utilities::scan_for_peaks::scan_for_peaks;
use crate::utilities::shape_filter::{Candidate, ShapeFilter};
use crate::utilities::structs::{DataXY, Peak};

#[derive(Clone, Copy, Debug)]
pub struct ArtifactFilter {
    pub min_r2: f64,
    pub shape: PeakShape,
}

impl Default for ArtifactFilter {
    fn default() -> Self {
        Self {
            min_r2: 0.0,
            shape: PeakShape::EMG,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum NoiseMethod {
    #[default]
    FindNoiseLevel,
    SanPlot,
}

#[derive(Clone, Copy, Debug)]
pub struct PeakFilter {
    pub min_integral: Option<f64>,
    pub min_peak_width_points: Option<usize>,
    pub min_intensity: Option<f64>,
    pub noise: Option<f64>,
    pub auto_noise: Option<bool>,
    pub auto_baseline: Option<bool>,
    pub allow_overlap: Option<bool>,
    pub min_snr: Option<f64>,
    pub noise_method: Option<NoiseMethod>,
    pub kernel_size: Option<usize>,
}

impl Default for PeakFilter {
    fn default() -> Self {
        Self {
            min_integral: None,
            min_peak_width_points: Some(5),
            min_intensity: None,
            noise: None,
            auto_noise: Some(false),
            auto_baseline: Some(false),
            allow_overlap: Some(false),
            min_snr: Some(1.0),
            noise_method: None,
            kernel_size: None,
        }
    }
}

#[derive(Clone, Debug)]
struct PeakCandidate {
    from: f64,
    to: f64,
    from_idx: usize,
    to_idx: usize,
    rt: f64,
    integral: f64,
    intensity: f64,
    n_points: usize,
    r2: f64,
    noise: f64,
}

impl From<PeakCandidate> for Peak {
    fn from(c: PeakCandidate) -> Self {
        Peak {
            from: c.from,
            to: c.to,
            rt: c.rt,
            integral: c.integral,
            intensity: c.intensity,
            r2: if c.r2.is_finite() { Some(c.r2) } else { None },
            n_points: c.n_points,
            noise: c.noise,
        }
    }
}

impl Candidate for PeakCandidate {
    fn bounds(&self) -> (usize, usize) {
        (self.from_idx, self.to_idx)
    }

    fn set_score(&mut self, r2: f64) {
        self.r2 = r2;
    }
}

#[derive(Clone, Debug)]
pub struct FindPeaksOptions {
    pub boundaries: Option<BoundariesOptions>,
    pub filter: Option<PeakFilter>,
    pub baseline: Option<BaselineOptions>,
    pub artifact_filter: Option<ArtifactFilter>,
}

impl Default for FindPeaksOptions {
    fn default() -> Self {
        Self {
            boundaries: Some(BoundariesOptions::default()),
            filter: Some(PeakFilter::default()),
            baseline: Some(BaselineOptions::default()),
            artifact_filter: Some(ArtifactFilter::default()),
        }
    }
}

pub fn find_peaks(data: &DataXY, options: Option<FindPeaksOptions>) -> Vec<Peak> {
    let o = options.unwrap_or_default();
    let filter = o.filter.unwrap_or_default();
    let base_opts = o.baseline.unwrap_or_default();

    let n = data.y.len();
    if n == 0 {
        return Vec::new();
    }

    let auto_baseline = filter.auto_baseline.unwrap_or(false);
    let auto_noise = filter.auto_noise.unwrap_or(false);

    let baseline: Vec<f64> = if auto_baseline {
        let mut b = base_opts;
        b.edge_slope_level = Some(0);
        calculate_baseline(&data.y, b)
    } else {
        vec![0.0; n]
    };

    let y_center: Vec<f64> = (0..n).map(|i| (data.y[i] - baseline[i]).max(0.0)).collect();

    let normalized_data = DataXY {
        x: data.x.clone(),
        y: y_center,
    };

    let noise: f64 = if auto_noise {
        match filter.noise_method.unwrap_or_default() {
            NoiseMethod::FindNoiseLevel => find_noise_level(&normalized_data.y).intensity,
            NoiseMethod::SanPlot => find_noise_level_san_plot(&normalized_data.y).intensity,
        }
    } else {
        filter.noise.unwrap_or_default()
    };

    let positions = scan_for_peaks(&normalized_data);
    if positions.is_empty() {
        return Vec::new();
    }

    let mut bopt = o.boundaries.unwrap_or_default();
    bopt.noise = noise;

    let smoothed_signal = get_smoothed_signal(&normalized_data, filter.kernel_size);
    let boundary_input = smoothed_signal.as_ref().unwrap_or(&normalized_data);

    let mut candidates: Vec<PeakCandidate> = Vec::with_capacity(positions.len());
    for seed_rt in positions {
        let b = get_boundaries(boundary_input, seed_rt, Some(bopt));
        let seed_idx = closest_index(&normalized_data.x, seed_rt);
        let apex = apex_in_window(&normalized_data, &b);
        let (rt, apex_y) = if let Some(t) = apex {
            t
        } else {
            (normalized_data.x[seed_idx], normalized_data.y[seed_idx])
        };
        if apex_y <= noise {
            continue;
        }

        match (b.from.index, b.from.value, b.to.index, b.to.value) {
            (Some(fi), Some(fx), Some(ti), Some(tx)) if fi < ti => {
                let (integral, intensity) = xy_integration(&data.x[fi..=ti], &data.y[fi..=ti]);
                candidates.push(PeakCandidate {
                    from: fx,
                    to: tx,
                    from_idx: fi,
                    to_idx: ti,
                    rt,
                    integral,
                    intensity,
                    n_points: ti - fi + 1,
                    r2: f64::NAN,
                    noise,
                });
            }
            _ => {}
        }
    }
    if candidates.is_empty() {
        return Vec::new();
    }

    if let Some(artifact_opts) = o.artifact_filter {
        ShapeFilter::new(artifact_opts.min_r2, artifact_opts.shape).retain(
            &mut candidates,
            &data.x,
            &data.y,
            &baseline,
        );
    }

    if candidates.is_empty() {
        return Vec::new();
    }

    let mut peaks = filter_peak_candidates(candidates, filter);

    peaks.sort_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap_or(Ordering::Equal));
    peaks = dedupe_near_identical(peaks);

    if !peaks.is_empty() {
        let mut cutoff = 0.0_f64;
        if noise > 0.0 {
            let sn_mult = filter.min_snr.unwrap_or(1.0);
            cutoff = sn_mult * noise;
        }
        if let Some(user_int) = filter.min_intensity
            && user_int > cutoff
        {
            cutoff = user_int;
        }

        if cutoff > 0.0 {
            peaks.retain(|p| p.intensity > cutoff);
        }
    }

    if peaks.len() > 1 {
        peaks = suppress_contained_peaks(data, peaks);
    }
    peaks
}

fn get_smoothed_signal(input: &DataXY, kernel_size: Option<usize>) -> Option<DataXY> {
    let window = get_smoothing_window(kernel_size?, input.y.len())?;
    let smoothed = sgg(
        &input.y,
        &input.x,
        SggOptions {
            window_size: window,
            derivative: 0,
            polynomial: 3,
        },
    );
    Some(DataXY {
        x: input.x.clone(),
        y: smoothed.into_iter().map(keep_above_zero).collect(),
    })
}

fn get_smoothing_window(requested: usize, length: usize) -> Option<usize> {
    let mut window = requested.min(length);
    if window.is_multiple_of(2) {
        window -= 1;
    }
    (window >= 5).then_some(window)
}

fn keep_above_zero(value: f64) -> f64 {
    value.max(0.0)
}

pub fn apex_in_window(data: &DataXY, b: &Boundaries) -> Option<(f64, f64)> {
    const APEX_MIN_PAD: usize = 3;

    let n = data.y.len();
    if n == 0 {
        return None;
    }

    let mut l = b.from.index?;
    let mut r = b.to.index?;
    if l >= n || r >= n {
        return None;
    }
    if l > r {
        core::mem::swap(&mut l, &mut r);
    }

    let need = 2 * APEX_MIN_PAD + 1;
    let width = r - l + 1;
    if width < need {
        let c = (l + r) / 2;
        l = c.saturating_sub(APEX_MIN_PAD);
        r = (c + APEX_MIN_PAD).min(n - 1);
    }

    if l >= r {
        return None;
    }

    let mut max_off = 0usize;
    let mut max_y = f64::NEG_INFINITY;
    let mut k = 0usize;
    while k <= (r - l) {
        let y = data.y[l + k];
        if y > max_y {
            max_y = y;
            max_off = k;
        }
        k += 1;
    }
    let i = l + max_off;
    Some((data.x[i], max_y))
}

fn filter_peak_candidates(peaks: Vec<PeakCandidate>, opt: PeakFilter) -> Vec<Peak> {
    let mut out: Vec<Peak> = Vec::with_capacity(peaks.len());

    let min_intensity = opt.min_intensity;
    let min_width = opt.min_peak_width_points;

    for p in peaks {
        let mut pass = true;
        if let Some(mi) = min_intensity
            && p.intensity < mi
        {
            pass = false;
        }

        if pass
            && let Some(w) = min_width
            && p.n_points < w
        {
            pass = false;
        }

        if pass {
            out.push(Peak::from(p));
        }
    }
    out
}

fn dedupe_near_identical(peaks: Vec<Peak>) -> Vec<Peak> {
    let m = peaks.len();
    if m <= 1 {
        return peaks;
    }

    let eps_rt = 1e-6;
    let eps_w = 1e-6;

    let mut out: Vec<Peak> = Vec::with_capacity(m);
    let mut i = 0usize;

    while i < m {
        let p = &peaks[i];
        let mut best_idx = i;
        let mut j = i + 1;

        while j < m {
            let q = &peaks[j];
            let same = (p.from - q.from).abs() <= eps_w
                && (p.to - q.to).abs() <= eps_w
                && (p.rt - q.rt).abs() <= eps_rt;
            if same {
                if q.intensity > peaks[best_idx].intensity {
                    best_idx = j;
                }
                j += 1;
            } else {
                break;
            }
        }

        out.push(peaks[best_idx]);
        i = j;
    }

    out
}

fn suppress_contained_peaks(data: &DataXY, mut peaks: Vec<Peak>) -> Vec<Peak> {
    let m = peaks.len();
    if m <= 1 {
        return peaks;
    }

    const OVERLAP_FRAC_SHOULDER: f64 = 0.70;
    const INTENSITY_RATIO_SHOULDER: f64 = 0.10;
    const MERGE_SHOULDERS: bool = true;

    let mut order: Vec<usize> = (0..m).collect();
    order.sort_by(|&i, &j| {
        peaks[j]
            .intensity
            .partial_cmp(&peaks[i].intensity)
            .unwrap_or(Ordering::Equal)
    });

    let mut keep = vec![true; m];

    let mut a_rank = 0usize;
    while a_rank < order.len() {
        let ia = order[a_rank];
        if keep[ia] {
            let la = peaks[ia].from;
            let ra = peaks[ia].to;

            let mut b_rank = a_rank + 1;
            while b_rank < order.len() {
                let ib = order[b_rank];
                if keep[ib] {
                    let lb = peaks[ib].from;
                    let rb = peaks[ib].to;

                    let wb = (rb - lb).abs();
                    if wb <= 0.0 {
                        keep[ib] = false;
                        b_rank += 1;
                        continue;
                    }

                    let l = if la > lb { la } else { lb };
                    let r = if ra < rb { ra } else { rb };
                    let overlap = (r - l).max(0.0);
                    let frac = overlap / wb;

                    let apex_b_inside_a = peaks[ib].rt >= la && peaks[ib].rt <= ra;
                    let rel = peaks[ib].intensity / peaks[ia].intensity;

                    let enough_overlap = frac >= OVERLAP_FRAC_SHOULDER;
                    let inside_and_small = apex_b_inside_a && rel <= INTENSITY_RATIO_SHOULDER;
                    let is_shoulder = enough_overlap || inside_and_small;

                    if is_shoulder {
                        if MERGE_SHOULDERS {
                            let new_from = if la < lb { la } else { lb };
                            let new_to = if ra > rb { ra } else { rb };

                            let li = closest_index(&data.x, new_from);
                            let ri = closest_index(&data.x, new_to);

                            let lidx = if li <= ri { li } else { ri };
                            let ridx = if li <= ri { ri } else { li };

                            let (integral, _) =
                                xy_integration(&data.x[lidx..=ridx], &data.y[lidx..=ridx]);

                            peaks[ia].from = new_from;
                            peaks[ia].to = new_to;
                            peaks[ia].integral = integral;
                            peaks[ia].n_points = ridx - lidx + 1;
                        }

                        keep[ib] = false;
                    }
                }
                b_rank += 1;
            }
        }
        a_rank += 1;
    }

    let mut out = Vec::<Peak>::with_capacity(m);
    for i in 0..m {
        if keep[i] {
            out.push(peaks[i]);
        }
    }
    out.sort_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap_or(Ordering::Equal));
    out
}
