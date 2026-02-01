use crate::utilities::calculate_eic::{
    CentroidScan, EicOptions, TimeUnit, collect_ms1_scans, compute_eic_for_mz, lower_bound,
    upper_bound, with_eic_apex_intensity,
};
use crate::utilities::find_peaks::{FindPeaksOptions, find_peaks};
use crate::utilities::structs::{DataXY, FromTo, Peak};
use octo::MzML;
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
use std::cmp::Ordering;
use std::mem;
use std::time::Instant;

#[derive(Clone, Debug)]
pub struct Feature {
    pub mz: f64,
    pub rt: f64,
    pub intensity: f64,
    pub from: f64,
    pub to: f64,
    pub np: usize,
}

pub struct MzScanGrid {
    pub mz_min: f64,
    pub mz_max: f64,
    pub step_size: f64,
}

impl Default for MzScanGrid {
    fn default() -> Self {
        Self {
            mz_min: 70.0,
            mz_max: 1000.0,
            step_size: 0.005,
        }
    }
}

pub struct FindFeaturesOptions {
    pub scan_eic_options: Option<EicOptions>,
    pub eic_options: Option<EicOptions>,
    pub find_peaks: Option<FindPeaksOptions>,
    pub mz_scan_grid: Option<MzScanGrid>,
    pub scan_width_threshold: Option<usize>,
}

impl Default for FindFeaturesOptions {
    fn default() -> Self {
        Self {
            scan_eic_options: Some(EicOptions {
                ppm_tolerance: 10.0,
                mz_tolerance: 0.003,
                ..Default::default()
            }),
            eic_options: Some(EicOptions {
                ppm_tolerance: 20.0,
                mz_tolerance: 0.005,
                ..Default::default()
            }),
            find_peaks: Some(FindPeaksOptions::default()),
            mz_scan_grid: Some(MzScanGrid::default()),
            scan_width_threshold: Some(5),
        }
    }
}

pub fn find_features(
    mzml: &MzML,
    time_window: FromTo,
    options: Option<FindFeaturesOptions>,
    cores: usize,
) -> Vec<Feature> {
    let t0 = Instant::now();

    let opts = options.unwrap_or_default();
    let scan_eic_options = opts.scan_eic_options.unwrap_or(EicOptions {
        ppm_tolerance: 10.0,
        mz_tolerance: 0.003,
        ..Default::default()
    });
    let eic_options = opts.eic_options.unwrap_or(EicOptions {
        ppm_tolerance: 20.0,
        mz_tolerance: 0.005,
        ..Default::default()
    });
    let find_peak_options = opts.find_peaks.unwrap_or_default();
    let mz_scan_grid_options = opts.mz_scan_grid.unwrap_or_default();
    let scan_width_threshold = opts.scan_width_threshold.unwrap_or(5);

    if !mz_scan_grid_options.step_size.is_finite() || mz_scan_grid_options.step_size <= 0.0 {
        panic!(
            "[panic] step_size must be > 0 in Da, got {}",
            mz_scan_grid_options.step_size
        );
    }

    let grid = build_mz_grid(
        mz_scan_grid_options.mz_min,
        mz_scan_grid_options.mz_max,
        mz_scan_grid_options.step_size,
    );
    if grid.is_empty() {
        panic!("[panic] empty grid");
    }
    if grid.len() > 2_000_000 {
        panic!("[panic] grid too large: {}", grid.len());
    }

    let (rts, scans) = collect_ms1_scans(mzml, time_window, TimeUnit::Minutes);
    if scans.is_empty() {
        panic!("[panic] no scans in time window");
    }
    let time = rts;

    let pool = ThreadPoolBuilder::new()
        .num_threads(cores.max(1))
        .thread_name(|i| format!("ff-{}", i))
        .build()
        .expect("failed to build rayon pool");

    pool.install(|| {
        let t1 = Instant::now();

        let mut coarse_find_opts = find_peak_options.clone();
        let mut coarse_filter = coarse_find_opts.filter_peaks_options.unwrap_or_default();
        coarse_filter.width_threshold = Some(scan_width_threshold);
        coarse_find_opts.filter_peaks_options = Some(coarse_filter);

        let masses: Vec<f64> = grid
            .par_iter()
            .map_init(
                || {
                    let bins_scratch: Vec<f64> = Vec::new();
                    bins_scratch
                },
                |bins_scratch, &mz_target| {
                    let y0 = compute_eic_for_mz(&scans, time.len(), mz_target, scan_eic_options);
                    let data = DataXY {
                        x: time.clone(),
                        y: y0,
                    };

                    let (rt_from, rt_to) = find_peaks(&data, Some(coarse_find_opts.clone()))
                        .into_iter()
                        .max_by(|a, b| {
                            a.intensity
                                .partial_cmp(&b.intensity)
                                .unwrap_or(Ordering::Equal)
                        })
                        .map(|p| (p.from, p.to))
                        .unwrap_or((time_window.from, time_window.to));

                    let refined_opt = refine_mz_for_peak(
                        &scans,
                        &time,
                        mz_target,
                        rt_from,
                        rt_to,
                        eic_options,
                        bins_scratch,
                    );

                    match refined_opt {
                        Some(refined_mz) => weighted_centroid(
                            &scans,
                            &time,
                            rt_from,
                            rt_to,
                            refined_mz,
                            eic_options,
                        ),
                        None => None,
                    }
                },
            )
            .filter_map(|v| v)
            .collect();

        if masses.is_empty() {
            panic!("[panic] refine_mz_for_peak returned empty list");
        }

        let unique_masses = dedup_masses_dynamic(masses, eic_options);

        let t2 = Instant::now();

        let final_w = find_peak_options
            .filter_peaks_options
            .as_ref()
            .and_then(|o| o.width_threshold)
            .unwrap_or(scan_width_threshold)
            .max(0);

        let features_raw: Vec<Feature> = unique_masses
            .par_iter()
            .map(|&mz| {
                let y = compute_eic_for_mz(&scans, time.len(), mz, eic_options);
                let data_local = DataXY { x: time.clone(), y };

                let mut peaks = find_peaks(&data_local, Some(find_peak_options.clone()));
                if peaks.is_empty() {
                    return Vec::<Feature>::new();
                }
                for i in 0..peaks.len() {
                    let p = mem::take(&mut peaks[i]);
                    peaks[i] = with_eic_apex_intensity(&data_local.x, &data_local.y, p);
                }
                sort_peaks_desc(&mut peaks);
                let mut out = Vec::with_capacity(peaks.len());
                for p in peaks {
                    if final_w > 0 && p.np < final_w {
                        continue;
                    }
                    out.push(Feature {
                        mz,
                        rt: p.rt,
                        intensity: p.intensity,
                        from: p.from,
                        to: p.to,
                        np: p.np,
                    });
                }
                out
            })
            .reduce(Vec::new, |mut a, mut b| {
                if !b.is_empty() {
                    a.append(&mut b);
                }
                a
            });

        let t3 = Instant::now();

        let mut features = dedup_features_dynamic_ppm(features_raw, eic_options, 0.80);

        features.sort_unstable_by(|a, b| {
            a.rt.partial_cmp(&b.rt)
                .unwrap_or(Ordering::Equal)
                .then_with(|| {
                    b.intensity
                        .partial_cmp(&a.intensity)
                        .unwrap_or(Ordering::Equal)
                })
                .then_with(|| a.mz.partial_cmp(&b.mz).unwrap_or(Ordering::Equal))
        });

        let _ = (t0, t1, t2, t3);
        features
    })
}

fn dedup_masses_dynamic(mut masses: Vec<f64>, opts: EicOptions) -> Vec<f64> {
    if masses.is_empty() {
        return masses;
    }
    masses.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let mut out = Vec::with_capacity(masses.len());
    let mut last = masses[0];
    out.push(last);
    for &m in masses.iter().skip(1) {
        if !mass_close_dynamic(m, last, opts) {
            out.push(m);
            last = m;
        }
    }
    out
}

pub fn weighted_centroid(
    scans: &[CentroidScan],
    rt: &[f64],
    rt_from: f64,
    rt_to: f64,
    approx_mz: f64,
    opts: EicOptions,
) -> Option<f64> {
    let n = rt.len();
    if scans.len() != n || n == 0 {
        return None;
    }
    let ppm = if opts.ppm_tolerance > 0.0 {
        opts.ppm_tolerance * 1e-6
    } else {
        0.0
    };
    let tol = opts.mz_tolerance.max(approx_mz * ppm);
    if !(tol > 0.0) || !tol.is_finite() {
        return None;
    }
    let lo_mz = approx_mz - tol;
    let hi_mz = approx_mz + tol;

    let mut start = lower_bound(rt, rt_from);
    let mut end = upper_bound(rt, rt_to);
    if start > n {
        start = n;
    }
    if end > n {
        end = n;
    }
    if start >= end {
        return None;
    }

    let mut w_sum = 0.0f64;
    let mut mw_sum = 0.0f64;

    for s in start..end {
        let mzs = &scans[s].mz;
        let its = &scans[s].intensity;
        if mzs.is_empty() || mzs.len() != its.len() {
            continue;
        }
        let mut j = lower_bound(mzs, lo_mz);
        let mlen = mzs.len();
        while j < mlen {
            let m = mzs[j];
            if m > hi_mz {
                break;
            }
            let it = its[j];
            if it.is_finite() && it > 0.0 {
                w_sum += it;
                mw_sum += m * it;
            }
            j += 1;
        }
    }

    if w_sum > 0.0 {
        Some(mw_sum / w_sum)
    } else {
        None
    }
}

fn mass_close_dynamic(a: f64, b: f64, opts: EicOptions) -> bool {
    let d = (a - b).abs();
    let c = 0.5 * (a + b).abs();
    let tol_ppm = if opts.ppm_tolerance > 0.0 {
        (opts.ppm_tolerance * 1e-6) * c
    } else {
        0.0
    };
    let tol = tol_ppm.max(opts.mz_tolerance.max(0.0));
    d <= tol
}

pub fn refine_mz_for_peak(
    scans: &[CentroidScan],
    rt: &[f64],
    approx: f64,
    rt_from: f64,
    rt_to: f64,
    opts: EicOptions,
    bins_scratch: &mut Vec<f64>,
) -> Option<f64> {
    let i0 = lower_bound(rt, rt_from);
    let i1 = upper_bound(rt, rt_to).min(scans.len());
    if i0 >= i1 {
        return None;
    }

    let tol_ppm = if opts.ppm_tolerance > 0.0 {
        (opts.ppm_tolerance * 1e-6) * approx
    } else {
        0.0
    };
    let tol = tol_ppm.max(opts.mz_tolerance.max(0.0));
    if !tol.is_finite() || tol <= 0.0 {
        return None;
    }

    let lo = approx - tol;
    let hi = approx + tol;
    let span = hi - lo;
    if !span.is_finite() || span <= 0.0 {
        return None;
    }

    let bin_da = (span / 400.0).max(1e-9);
    if !bin_da.is_finite() || bin_da <= 0.0 {
        return None;
    }

    let n_bins = ((span / bin_da).ceil() as usize).saturating_add(1);
    if n_bins == 0 || n_bins > 20_000_000 {
        return None;
    }

    if bins_scratch.len() < n_bins {
        bins_scratch.resize(n_bins, 0.0);
    } else {
        for v in &mut bins_scratch[..n_bins] {
            *v = 0.0;
        }
    }
    let bins = &mut bins_scratch[..n_bins];

    for s in i0..i1 {
        let mzs = &scans[s].mz;
        let ints = &scans[s].intensity;
        let mut j = lower_bound(mzs, lo);
        let mlen = mzs.len();
        while j < mlen {
            let m = mzs[j];
            if m > hi {
                break;
            }
            let it = ints[j];
            if it.is_finite() && it > 0.0 && m.is_finite() {
                let idx_f = (m - lo) / bin_da;
                if idx_f.is_finite() {
                    let idx = idx_f.floor() as isize;
                    if idx >= 0 {
                        let k = idx as usize;
                        if k < n_bins {
                            bins[k] += it;
                        }
                    }
                }
            }
            j += 1;
        }
    }

    let mut any = false;
    for &v in bins.iter() {
        if v > 0.0 {
            any = true;
            break;
        }
    }
    if !any {
        return Some(approx);
    }

    let mut w_bins = (tol / bin_da).round() as isize;
    if w_bins < 1 {
        w_bins = 1;
    }
    let w_bins = w_bins as usize;

    let mut best_i: usize = 0;
    if n_bins >= w_bins {
        let mut window_sum = 0.0;
        for i in 0..w_bins {
            window_sum += bins[i];
        }
        let mut best_sum = window_sum;
        best_i = 0;
        let last = n_bins - w_bins;
        for i in 1..=last {
            window_sum += bins[i + w_bins - 1];
            window_sum -= bins[i - 1];
            if window_sum > best_sum {
                best_sum = window_sum;
                best_i = i;
            }
        }
    }

    let start = best_i;
    let end = (best_i + w_bins).min(n_bins);
    let mut max_v = -1.0f64;
    let mut max_k = start;
    for k in start..end {
        if bins[k] > max_v {
            max_v = bins[k];
            max_k = k;
        }
    }

    let mz = lo + (max_k as f64 + 0.5) * bin_da;
    if !mz.is_finite() {
        return Some(approx);
    }
    Some(mz)
}

fn build_mz_grid(start: f64, end: f64, step_da: f64) -> Vec<f64> {
    let (lo, hi) = if start <= end {
        (start, end)
    } else {
        (end, start)
    };
    if !lo.is_finite() || !hi.is_finite() || hi <= lo {
        return Vec::new();
    }
    if step_da <= 0.0 || !step_da.is_finite() {
        return vec![lo, hi];
    }
    let approx_n = (((hi - lo) / step_da).floor() as usize).saturating_add(2);
    let mut xs = Vec::with_capacity(approx_n);
    let mut m = lo;
    while m <= hi {
        xs.push(m);
        m += step_da;
    }
    const EPS: f64 = 1e-9;
    if let Some(last) = xs.last_mut() {
        if (hi - *last).abs() > EPS {
            xs.push(hi);
        } else {
            *last = hi;
        }
    } else {
        xs.push(hi);
    }
    xs
}

fn sort_peaks_desc(xs: &mut Vec<Peak>) {
    xs.sort_unstable_by(|a, b| {
        b.intensity
            .partial_cmp(&a.intensity)
            .unwrap_or(Ordering::Equal)
            .then_with(|| {
                b.integral
                    .partial_cmp(&a.integral)
                    .unwrap_or(Ordering::Equal)
            })
    });
}

fn nearly_eq(a: f64, b: f64, eps: f64) -> bool {
    (a - b).abs() <= eps
}

fn mass_close_for_dedup(a: f64, b: f64, eic: EicOptions) -> bool {
    let c = 0.5 * (a + b).abs();
    let tol_ppm = if eic.ppm_tolerance > 0.0 {
        (eic.ppm_tolerance * 1e-6) * c
    } else {
        0.0
    };
    let base = tol_ppm.max(eic.mz_tolerance.max(0.0));
    let tol = base * 1.2;
    (a - b).abs() <= tol
}

fn rt_overlap_fraction(a_from: f64, a_to: f64, b_from: f64, b_to: f64) -> f64 {
    let l = a_from.max(b_from);
    let r = a_to.min(b_to);
    let overlap = (r - l).max(0.0);
    let wa = (a_to - a_from).max(0.0);
    let wb = (b_to - b_from).max(0.0);
    let base = wa.max(wb).max(f64::EPSILON);
    overlap / base
}

fn rt_overlap_fraction_min(a_from: f64, a_to: f64, b_from: f64, b_to: f64) -> f64 {
    let l = a_from.max(b_from);
    let r = a_to.min(b_to);
    let overlap = (r - l).max(0.0);
    let wa = (a_to - a_from).max(0.0);
    let wb = (b_to - b_from).max(0.0);
    let base = wa.min(wb).max(f64::EPSILON);
    overlap / base
}

fn better(a: &Feature, b: &Feature) -> bool {
    if a.np != b.np {
        return a.np > b.np;
    }
    if a.intensity != b.intensity {
        return a.intensity > b.intensity;
    }
    let wa = (a.to - a.from).abs();
    let wb = (b.to - b.from).abs();
    if wa != wb {
        return wa < wb;
    }
    a.mz <= b.mz
}

fn dedup_features_dynamic_ppm(
    mut xs: Vec<Feature>,
    eic: EicOptions,
    min_rt_overlap: f64,
) -> Vec<Feature> {
    if xs.is_empty() {
        return xs;
    }

    xs.sort_unstable_by(|a, b| {
        a.rt.partial_cmp(&b.rt)
            .unwrap_or(Ordering::Equal)
            .then_with(|| a.from.partial_cmp(&b.from).unwrap_or(Ordering::Equal))
            .then_with(|| a.to.partial_cmp(&b.to).unwrap_or(Ordering::Equal))
            .then_with(|| a.mz.partial_cmp(&b.mz).unwrap_or(Ordering::Equal))
    });

    let eps_rt = 1e-6;
    let eps_w = 1e-6;

    let mut out: Vec<Feature> = Vec::with_capacity(xs.len());
    let mut cluster: Vec<Feature> = Vec::new();

    for f in xs.into_iter() {
        if cluster.is_empty() {
            cluster.push(f);
            continue;
        }
        let last = cluster.last().unwrap();

        let same_window = nearly_eq(f.from, last.from, eps_w)
            && nearly_eq(f.to, last.to, eps_w)
            && nearly_eq(f.rt, last.rt, eps_rt);

        let ovl_max = rt_overlap_fraction(last.from, last.to, f.from, f.to);
        let ovl_min = rt_overlap_fraction_min(last.from, last.to, f.from, f.to);

        let mass_close = mass_close_for_dedup(f.mz, last.mz, eic);
        let same_apex = nearly_eq(f.rt, last.rt, eps_rt);

        let close_in_mass_and_time =
            mass_close && (same_apex || ovl_max >= min_rt_overlap || ovl_min >= 0.95);

        if same_window || close_in_mass_and_time {
            cluster.push(f);
        } else {
            let mut best = cluster[0].clone();
            for c in &cluster[1..] {
                if better(c, &best) {
                    best = c.clone();
                }
            }
            out.push(best);
            cluster.clear();
            cluster.push(f);
        }
    }

    if !cluster.is_empty() {
        let mut best = cluster[0].clone();
        for c in &cluster[1..] {
            if better(c, &best) {
                best = c.clone();
            }
        }
        out.push(best);
    }

    out.sort_unstable_by(|a, b| {
        a.mz.partial_cmp(&b.mz)
            .unwrap_or(Ordering::Equal)
            .then_with(|| a.rt.partial_cmp(&b.rt).unwrap_or(Ordering::Equal))
    });

    let mut final_out: Vec<Feature> = Vec::with_capacity(out.len());
    for f in out.into_iter() {
        if let Some(g) = final_out.last_mut() {
            let same_window = nearly_eq(f.from, g.from, eps_w)
                && nearly_eq(f.to, g.to, eps_w)
                && nearly_eq(f.rt, g.rt, eps_rt);

            let ovl_max = rt_overlap_fraction(g.from, g.to, f.from, f.to);
            let ovl_min = rt_overlap_fraction_min(g.from, g.to, f.from, f.to);

            let mass_close = mass_close_for_dedup(f.mz, g.mz, eic);
            let same_apex = nearly_eq(f.rt, g.rt, eps_rt);

            let close_in_mass_and_time =
                mass_close && (same_apex || ovl_max >= min_rt_overlap || ovl_min >= 0.95);

            if same_window || close_in_mass_and_time {
                if better(&f, g) {
                    *g = f;
                }
                continue;
            }
        }
        final_out.push(f);
    }
    final_out
}
