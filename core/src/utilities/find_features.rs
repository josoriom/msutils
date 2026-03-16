use octo::MzML;
use serde::Serialize;
use std::{
    cmp::Ordering,
    collections::HashSet,
    error::Error,
    fmt::{Display, Formatter},
    mem,
};

use crate::utilities::{
    calculate_eic::{
        CentroidScan, EicOptions, TimeUnit, collect_scans, compute_eic_for_mz, lower_bound,
        with_eic_apex_intensity,
    },
    find_peaks::{FindPeaksOptions, find_peaks},
    structs::{DataXY, FromTo, Peak},
};

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use rayon::{ThreadPoolBuilder, prelude::*};

#[derive(Clone, Debug)]
pub struct MzTolerance {
    pub mz_abs: f64,
    pub ppm: f64,
}

#[derive(Debug)]
pub enum FeatureError {
    InvalidStepSize(f64),
    EmptyGrid,
    GridTooLarge(usize),
    NoScans,
    NoCandidateMasses,
}

impl Display for FeatureError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidStepSize(v) => write!(f, "invalid step size: {}", v),
            Self::EmptyGrid => write!(f, "empty grid"),
            Self::GridTooLarge(n) => write!(f, "grid too large: {}", n),
            Self::NoScans => write!(f, "no scans in time window"),
            Self::NoCandidateMasses => write!(f, "no candidate masses found"),
        }
    }
}

impl Error for FeatureError {}

#[derive(Clone, Debug, Serialize, Default)]
pub struct Feature {
    pub mz: f64,
    pub rt: f64,
    pub intensity: f64,
    pub from: f64,
    pub to: f64,
    pub np: usize,
    pub noise: f64,
    pub integral: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct FeaturePoint {
    pub rt: f64,
    pub mz: f64,
    pub intensity: f64,
}

impl MzTolerance {
    pub fn window_at(&self, mz: f64) -> (f64, f64) {
        let tol = self.tol_at(mz.abs());
        (mz - tol, mz + tol)
    }

    pub fn are_close(&self, a: f64, b: f64) -> bool {
        (a - b).abs() <= self.tol_at(0.5 * (a + b).abs())
    }

    pub fn are_close_to_ref(&self, a: f64, ref_mz: f64) -> bool {
        (a - ref_mz).abs() <= self.tol_at(ref_mz.abs())
    }

    fn tol_at(&self, center: f64) -> f64 {
        let tol_ppm = if self.ppm > 0.0 {
            (self.ppm * 1e-6) * center
        } else {
            0.0
        };
        tol_ppm.max(self.mz_abs)
    }
}

#[derive(Clone, Debug)]
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
            step_size: 0.006,
        }
    }
}

#[derive(Clone, Debug)]
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

pub struct FindFeaturesConfig {
    pub scan_eic_options: EicOptions,
    pub eic_options: EicOptions,
    pub find_peaks: FindPeaksOptions,
    pub mz_scan_grid: MzScanGrid,
    pub scan_width_threshold: usize,
}

impl From<FindFeaturesOptions> for FindFeaturesConfig {
    fn from(opts: FindFeaturesOptions) -> Self {
        Self {
            scan_eic_options: opts.scan_eic_options.unwrap_or(EicOptions {
                ppm_tolerance: 10.0,
                mz_tolerance: 0.003,
                ..Default::default()
            }),
            eic_options: opts.eic_options.unwrap_or(EicOptions {
                ppm_tolerance: 20.0,
                mz_tolerance: 0.005,
                ..Default::default()
            }),
            find_peaks: opts.find_peaks.unwrap_or_default(),
            mz_scan_grid: opts.mz_scan_grid.unwrap_or_default(),
            scan_width_threshold: opts.scan_width_threshold.unwrap_or(5),
        }
    }
}

pub fn find_features(
    mzml: &MzML,
    time_window: FromTo,
    options: Option<FindFeaturesOptions>,
    cores: usize,
) -> Result<Vec<Feature>, FeatureError> {
    let config = FindFeaturesConfig::from(options.unwrap_or_default());

    if !config.mz_scan_grid.step_size.is_finite() || config.mz_scan_grid.step_size <= 0.0 {
        return Err(FeatureError::InvalidStepSize(config.mz_scan_grid.step_size));
    }

    let grid = build_mz_grid(
        config.mz_scan_grid.mz_min,
        config.mz_scan_grid.mz_max,
        config.mz_scan_grid.step_size,
    );
    if grid.is_empty() {
        return Err(FeatureError::EmptyGrid);
    }
    if grid.len() > 2_000_000 {
        return Err(FeatureError::GridTooLarge(grid.len()));
    }

    let (time, scans) = collect_scans(mzml, time_window, TimeUnit::Minutes, 1, false);
    if scans.is_empty() {
        return Err(FeatureError::NoScans);
    }

    let masses = collect_candidate_masses(&scans, &time, &grid, &config, cores)?;
    let masses = deduplicate_masses(masses, config.eic_options);
    let features = extract_features_from_masses(&masses, &scans, &time, &config, cores);
    let features = deduplicate_features(features, config.eic_options, 0.05);
    Ok(sort_features(features))
}

fn build_coarse_opts(config: &FindFeaturesConfig) -> FindPeaksOptions {
    let mut coarse = config.find_peaks.clone();
    let mut filter = coarse.filter_peaks_options.unwrap_or_default();
    filter.width_threshold = Some(config.scan_width_threshold);
    coarse.filter_peaks_options = Some(filter);
    coarse
}

fn process_single_grid_point(
    scans: &[CentroidScan],
    time: &[f64],
    mz_target: f64,
    coarse_opts: &FindPeaksOptions,
    scan_eic_options: EicOptions,
) -> Option<f64> {
    let y = compute_eic_for_mz(scans, time.len(), mz_target, scan_eic_options);
    let data = DataXY {
        x: time.to_vec(),
        y,
    };
    let apex_rt = find_peaks(&data, Some(coarse_opts.clone()))
        .into_iter()
        .max_by(|a, b| {
            a.intensity
                .partial_cmp(&b.intensity)
                .unwrap_or(Ordering::Equal)
        })
        .map(|p| p.rt)?;
    max_intensity_centroid(scans, time, apex_rt, mz_target, scan_eic_options)
}

fn collect_candidate_masses(
    scans: &[CentroidScan],
    time: &[f64],
    grid: &[f64],
    config: &FindFeaturesConfig,
    cores: usize,
) -> Result<Vec<f64>, FeatureError> {
    let coarse_opts = build_coarse_opts(config);

    #[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
    let masses: Vec<f64> = grid
        .iter()
        .filter_map(|&mz| {
            process_single_grid_point(scans, time, mz, &coarse_opts, config.scan_eic_options)
        })
        .collect();

    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    let masses: Vec<f64> = ThreadPoolBuilder::new()
        .num_threads(cores.max(1))
        .thread_name(|i| format!("ff-{}", i))
        .build()
        .expect("failed to build rayon pool")
        .install(|| {
            grid.par_iter()
                .filter_map(|&mz| {
                    process_single_grid_point(
                        scans,
                        time,
                        mz,
                        &coarse_opts,
                        config.scan_eic_options,
                    )
                })
                .collect()
        });

    if masses.is_empty() {
        return Err(FeatureError::NoCandidateMasses);
    }
    Ok(masses)
}

pub(crate) fn deduplicate_masses(mut masses: Vec<f64>, opts: EicOptions) -> Vec<f64> {
    if masses.is_empty() {
        return masses;
    }
    masses.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));

    let tolerance = MzTolerance {
        mz_abs: opts.mz_tolerance.max(0.0),
        ppm: opts.ppm_tolerance,
    };

    let mut out: Vec<f64> = Vec::with_capacity(masses.len());
    let mut cluster = vec![masses[0]];

    for &m in masses.iter().skip(1) {
        let median = cluster[cluster.len() / 2];
        if tolerance.are_close(m, median) {
            cluster.push(m);
        } else {
            out.push(cluster[cluster.len() / 2]);
            cluster.clear();
            cluster.push(m);
        }
    }
    out.push(cluster[cluster.len() / 2]);
    out
}

fn extract_peaks_for_mass(
    mz: f64,
    scans: &[CentroidScan],
    time: &[f64],
    config: &FindFeaturesConfig,
    final_w: usize,
) -> Vec<Feature> {
    let y = compute_eic_for_mz(scans, time.len(), mz, config.eic_options);
    let data = DataXY {
        x: time.to_vec(),
        y,
    };
    let mut peaks = find_peaks(&data, Some(config.find_peaks.clone()));
    if peaks.is_empty() {
        return Vec::new();
    }
    for i in 0..peaks.len() {
        let p = mem::take(&mut peaks[i]);
        peaks[i] = with_eic_apex_intensity(&data.x, &data.y, p);
    }
    sort_peaks_desc(&mut peaks);
    peaks
        .into_iter()
        .filter(|p| final_w == 0 || p.np >= final_w)
        .map(|p| Feature {
            mz,
            rt: p.rt,
            intensity: p.intensity,
            from: p.from,
            to: p.to,
            np: p.np,
            integral: p.integral,
            noise: p.noise,
        })
        .collect()
}

fn extract_features_from_masses(
    masses: &[f64],
    scans: &[CentroidScan],
    time: &[f64],
    config: &FindFeaturesConfig,
    cores: usize,
) -> Vec<Feature> {
    let final_w = config
        .find_peaks
        .filter_peaks_options
        .as_ref()
        .and_then(|o| o.width_threshold)
        .unwrap_or(config.scan_width_threshold)
        .max(0);

    #[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
    let features: Vec<Feature> = masses
        .iter()
        .flat_map(|&mz| extract_peaks_for_mass(mz, scans, time, config, final_w))
        .collect();

    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    let features: Vec<Feature> = ThreadPoolBuilder::new()
        .num_threads(cores.max(1))
        .thread_name(|i| format!("ff-{}", i))
        .build()
        .expect("failed to build rayon pool")
        .install(|| {
            masses
                .par_iter()
                .map(|&mz| extract_peaks_for_mass(mz, scans, time, config, final_w))
                .reduce(Vec::new, |mut a, mut b| {
                    if !b.is_empty() {
                        a.append(&mut b);
                    }
                    a
                })
        });

    features
}

fn deduplicate_features(xs: Vec<Feature>, eic: EicOptions, rt_tolerance: f64) -> Vec<Feature> {
    if xs.is_empty() {
        return xs;
    }
    let points: Vec<(f64, f64, f64)> = xs.iter().map(|f| (f.rt, f.mz, f.intensity)).collect();
    let survivors: HashSet<(u64, u64, u64)> =
        dedup_points(points, eic.mz_tolerance, eic.ppm_tolerance, rt_tolerance)
            .into_iter()
            .map(|(rt, mz, intensity)| (rt.to_bits(), mz.to_bits(), intensity.to_bits()))
            .collect();
    xs.into_iter()
        .filter(|f| survivors.contains(&(f.rt.to_bits(), f.mz.to_bits(), f.intensity.to_bits())))
        .collect()
}

pub(crate) fn sort_features(mut features: Vec<Feature>) -> Vec<Feature> {
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
    features
}

pub fn max_intensity_centroid(
    scans: &[CentroidScan],
    rt: &[f64],
    apex_rt: f64,
    approx_mz: f64,
    opts: EicOptions,
) -> Option<f64> {
    let n = rt.len();
    if scans.len() != n || n == 0 {
        return None;
    }

    let tolerance = MzTolerance {
        mz_abs: opts.mz_tolerance,
        ppm: opts.ppm_tolerance,
    };
    let (lo_mz, hi_mz) = tolerance.window_at(approx_mz);
    if lo_mz >= hi_mz || !lo_mz.is_finite() || !hi_mz.is_finite() {
        return None;
    }

    let closest_idx = rt
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            (*a - apex_rt)
                .abs()
                .partial_cmp(&(*b - apex_rt).abs())
                .unwrap_or(Ordering::Equal)
        })
        .map(|(i, _)| i)?;

    let scan = &scans[closest_idx];
    if scan.mz.is_empty() || scan.mz.len() != scan.intensity.len() {
        return None;
    }

    let mut best_mz: Option<f64> = None;
    let mut best_intensity = 0.0f64;

    let mut j = lower_bound(&scan.mz, lo_mz);
    while j < scan.mz.len() {
        let m = scan.mz[j];
        if m > hi_mz {
            break;
        }
        let it = scan.intensity[j];
        if it.is_finite() && it > best_intensity {
            best_intensity = it;
            best_mz = Some(m);
        }
        j += 1;
    }

    best_mz
}

pub(crate) fn group_by_mz<F: Fn(f64, f64) -> bool>(
    points: Vec<FeaturePoint>,
    mz_close: F,
) -> Vec<Vec<FeaturePoint>> {
    let mut groups: Vec<Vec<FeaturePoint>> = Vec::new();
    let mut current: Vec<FeaturePoint> = Vec::new();

    for p in points {
        if current.is_empty() {
            current.push(p);
        } else {
            let pivot_mz = current[current.len() / 2].mz;
            if mz_close(p.mz, pivot_mz) {
                current.push(p);
            } else {
                groups.push(current);
                current = vec![p];
            }
        }
    }
    if !current.is_empty() {
        groups.push(current);
    }
    groups
}

pub(crate) fn pick_best_per_rt_cluster(
    group: Vec<FeaturePoint>,
    rt_tolerance: f64,
) -> Vec<FeaturePoint> {
    let mut rt_sorted = group;
    rt_sorted.sort_unstable_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap_or(Ordering::Equal));

    let mut out: Vec<FeaturePoint> = Vec::new();
    let mut cluster: Vec<FeaturePoint> = Vec::new();

    for p in rt_sorted {
        if cluster.is_empty() {
            cluster.push(p);
        } else {
            let pivot_rt = cluster[cluster.len() / 2].rt;
            if (p.rt - pivot_rt).abs() <= rt_tolerance {
                cluster.push(p);
            } else {
                let best = cluster
                    .iter()
                    .max_by(|a, b| {
                        a.intensity
                            .partial_cmp(&b.intensity)
                            .unwrap_or(Ordering::Equal)
                    })
                    .copied()
                    .unwrap();
                out.push(best);
                cluster = vec![p];
            }
        }
    }
    if !cluster.is_empty() {
        let best = cluster
            .iter()
            .max_by(|a, b| {
                a.intensity
                    .partial_cmp(&b.intensity)
                    .unwrap_or(Ordering::Equal)
            })
            .copied()
            .unwrap();
        out.push(best);
    }
    out
}

pub(crate) fn dedup_points(
    mut points: Vec<(f64, f64, f64)>,
    mz_tolerance: f64,
    ppm_tolerance: f64,
    rt_tolerance: f64,
) -> Vec<(f64, f64, f64)> {
    if points.is_empty() {
        return points;
    }

    points.sort_unstable_by(|a, b| {
        a.1.partial_cmp(&b.1)
            .unwrap_or(Ordering::Equal)
            .then_with(|| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal))
    });

    let mz_close = |a: f64, b: f64| -> bool {
        let c = 0.5 * (a + b).abs();
        let tol_ppm = if ppm_tolerance > 0.0 {
            (ppm_tolerance * 1e-6) * c
        } else {
            0.0
        };
        let tol = tol_ppm.max(mz_tolerance) * 1.2;
        (a - b).abs() <= tol
    };

    let feature_points: Vec<FeaturePoint> = points
        .into_iter()
        .map(|(rt, mz, intensity)| FeaturePoint { rt, mz, intensity })
        .collect();

    group_by_mz(feature_points, mz_close)
        .into_iter()
        .flat_map(|group| pick_best_per_rt_cluster(group, rt_tolerance))
        .map(|p| (p.rt, p.mz, p.intensity))
        .collect()
}

pub(crate) fn build_mz_grid(start: f64, end: f64, step_da: f64) -> Vec<f64> {
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
