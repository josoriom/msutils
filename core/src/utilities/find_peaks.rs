use std::cmp::Ordering;

use crate::utilities::calculate_baseline::{BaselineOptions, calculate_baseline};
use crate::utilities::closest_index;
use crate::utilities::find_noise_level::find_noise_level;
use crate::utilities::get_boundaries::{Boundaries, BoundariesOptions, get_boundaries};
use crate::utilities::scan_for_peaks::{ScanPeaksOptions, scan_for_peaks_across_windows};
use crate::utilities::structs::{DataXY, Peak};
use crate::utilities::utilities::xy_integration;

#[derive(Clone, Copy, Debug)]
pub struct FilterPeaksOptions {
    pub integral_threshold: Option<f64>,
    pub width_threshold: Option<usize>,
    pub intensity_threshold: Option<f64>,
    pub noise: Option<f64>,
    pub auto_noise: Option<bool>,
    pub auto_baseline: Option<bool>,
    pub allow_overlap: Option<bool>,
    pub sn_ratio: Option<f64>,
}

impl Default for FilterPeaksOptions {
    fn default() -> Self {
        Self {
            integral_threshold: None,
            width_threshold: Some(5),
            intensity_threshold: None,
            noise: None,
            auto_noise: Some(false),
            auto_baseline: Some(false),
            allow_overlap: Some(false),
            sn_ratio: Some(1.0),
        }
    }
}

#[derive(Clone, Debug)]
struct PeakCandidate {
    from: f64,
    to: f64,
    rt: f64,
    integral: f64,
    intensity: f64,
    number_of_points: usize,
    ratio: f64,
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
            ratio: c.ratio,
            np: c.number_of_points,
            noise: c.noise,
        }
    }
}

#[derive(Clone, Debug)]
pub struct FindPeaksOptions {
    pub get_boundaries_options: Option<BoundariesOptions>,
    pub filter_peaks_options: Option<FilterPeaksOptions>,
    pub scan_peaks_options: Option<ScanPeaksOptions>,
    pub baseline_options: Option<BaselineOptions>,
}

impl Default for FindPeaksOptions {
    fn default() -> Self {
        Self {
            get_boundaries_options: Some(BoundariesOptions::default()),
            filter_peaks_options: Some(FilterPeaksOptions::default()),
            scan_peaks_options: Some(ScanPeaksOptions::default()),
            baseline_options: Some(BaselineOptions::default()),
        }
    }
}

pub fn find_peaks(data: &DataXY, options: Option<FindPeaksOptions>) -> Vec<Peak> {
    let o = options.unwrap_or_default();
    let filter_opts = o.filter_peaks_options.unwrap_or_default();
    let base_opts = o.baseline_options.unwrap_or_default();

    let n = data.y.len();
    if n == 0 {
        return Vec::new();
    }

    let auto_baseline = filter_opts.auto_baseline.unwrap_or(false);
    let auto_noise = filter_opts.auto_noise.unwrap_or(false);
    if auto_noise && filter_opts.noise.is_some() {
        panic!("auto_noise=true cannot be used with noise");
    }

    let mut y_center = Vec::with_capacity(n);
    if auto_baseline {
        let mut b = base_opts.clone();
        b.level = Some(0);
        let floor = calculate_baseline(&data.y, b);
        y_center.extend((0..n).map(|i| {
            let v = data.y[i] - floor[i];
            if v > 0.0 { v } else { 0.0 }
        }));
    } else {
        y_center.extend((0..n).map(|i| {
            let v = data.y[i];
            if v > 0.0 { v } else { 0.0 }
        }));
    }

    let normalized_data = DataXY {
        x: data.x.clone(),
        y: y_center,
    };

    let noise = if auto_noise {
        find_noise_level(&normalized_data.y)
    } else {
        filter_opts.noise.unwrap_or(0.0).max(0.0)
    };

    let positions = scan_for_peaks_across_windows(&normalized_data, o.scan_peaks_options);
    if positions.is_empty() {
        return Vec::new();
    }

    let mut bopt = o.get_boundaries_options.unwrap_or_default();
    bopt.noise = noise;

    let mut candidates: Vec<PeakCandidate> = Vec::with_capacity(positions.len());
    for seed_rt in positions {
        let b = get_boundaries(&normalized_data, seed_rt, Some(bopt));
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
                    rt,
                    integral,
                    intensity,
                    number_of_points: ti - fi + 1,
                    ratio: 0.0,
                    noise,
                });
            }
            _ => {}
        }
    }
    if candidates.is_empty() {
        return Vec::new();
    }

    // println!("{candidates:?}");

    let mut peaks = filter_peak_candidates(candidates, filter_opts);

    peaks.sort_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap_or(Ordering::Equal));
    peaks = dedupe_near_identical(peaks);

    if !peaks.is_empty() {
        let mut cutoff = 0.0_f64;
        if noise > 0.0 {
            let sn_mult = filter_opts.sn_ratio.unwrap_or(1.0) as f64;
            cutoff = sn_mult * noise;
        }
        if let Some(user_int) = filter_opts.intensity_threshold {
            if user_int > cutoff {
                cutoff = user_int;
            }
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
        let tmp = l;
        l = r;
        r = tmp;
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

fn filter_peak_candidates(peaks: Vec<PeakCandidate>, opt: FilterPeaksOptions) -> Vec<Peak> {
    let mut out: Vec<Peak> = Vec::with_capacity(peaks.len());

    let min_intensity = opt.intensity_threshold;
    let min_width = opt.width_threshold;

    for p in peaks {
        let mut pass = true;

        if let Some(mi) = min_intensity {
            if p.intensity < mi {
                pass = false;
            }
        }
        if pass {
            if let Some(w) = min_width {
                if p.number_of_points <= w {
                    pass = false;
                }
            }
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

        out.push(peaks[best_idx].clone());
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

                    let is_shoulder = (apex_b_inside_a && rel <= INTENSITY_RATIO_SHOULDER)
                        || (apex_b_inside_a && frac >= OVERLAP_FRAC_SHOULDER)
                        || (frac >= OVERLAP_FRAC_SHOULDER);

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
                            peaks[ia].np = ridx - lidx + 1;
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
            out.push(peaks[i].clone());
        }
    }
    out.sort_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap_or(Ordering::Equal));
    out
}
