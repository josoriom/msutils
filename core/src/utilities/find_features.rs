use serde::Serialize;

use crate::utilities::structs::ser_finite_f64;
#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use std::sync::Arc;
use std::{
    cmp::Ordering,
    collections::HashSet,
    error::Error,
    fmt::{Display, Formatter},
};

use crate::utilities::{
    calculate_eic::{
        CentroidScan, EicOptions, EicReader, FastError, MS1_LEVEL, get_scan_times, lower_bound,
        mz_tolerance_for, read_mz_window, summed_intensity_in_window, with_eic_apex_intensity,
    },
    find_peaks::{FindPeaksOptions, find_peaks},
    parallel::run_with_cores,
    structs::{DataXY, FromTo},
};

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use crate::utilities::gpu::{
    GpuContext,
    processor::{FlattenedScans, GpuBatchOptions, GpuGridProcessor},
};

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use rayon::prelude::*;

#[derive(Clone, Debug)]
pub struct MzTolerance {
    pub mz_absolute: f64,
    pub ppm: f64,
}

impl MzTolerance {
    pub(crate) fn window_at(&self, mz: f64) -> (f64, f64) {
        let tol = self.tol_at(mz.abs());
        (mz - tol, mz + tol)
    }

    pub(crate) fn are_close(&self, a: f64, b: f64) -> bool {
        (a - b).abs() <= self.tol_at(0.5 * (a + b).abs())
    }

    pub(crate) fn are_close_to_ref(&self, a: f64, ref_mz: f64) -> bool {
        (a - ref_mz).abs() <= self.tol_at(ref_mz.abs())
    }

    fn tol_at(&self, center: f64) -> f64 {
        let tol_ppm = if self.ppm > 0.0 {
            (self.ppm * 1e-6) * center
        } else {
            0.0
        };
        tol_ppm.max(self.mz_absolute)
    }
}

#[derive(Debug)]
pub enum FeatureError {
    InvalidStepSize(f64),
    EmptyGrid,
    GridTooLarge(usize),
    NoScans,
    NoCandidateMasses,
    ReadFailed(FastError),
    TooManyCandidateMasses(usize),
}

impl Display for FeatureError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidStepSize(v) => write!(f, "invalid step size: {}", v),
            Self::EmptyGrid => write!(f, "empty grid"),
            Self::GridTooLarge(n) => write!(f, "grid too large: {}", n),
            Self::NoScans => write!(f, "no scans in time window"),
            Self::NoCandidateMasses => write!(f, "no candidate masses found"),
            Self::ReadFailed(e) => write!(f, "read failed: {}", e),
            Self::TooManyCandidateMasses(n) => write!(f, "too many candidate masses: {}", n),
        }
    }
}

impl Error for FeatureError {}

impl From<FastError> for FeatureError {
    fn from(e: FastError) -> Self {
        Self::ReadFailed(e)
    }
}

#[derive(Clone, Debug, Serialize, Default)]
pub struct Feature {
    #[serde(serialize_with = "ser_finite_f64")]
    pub mz: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    pub rt: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    pub from: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    pub to: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    pub intensity: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    pub integral: f64,
    pub n_points: usize,
    #[serde(skip)]
    pub noise: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct FeaturePoint {
    pub rt: f64,
    pub mz: f64,
    pub intensity: f64,
}

#[derive(Clone, Debug)]
pub struct MzScanGrid {
    pub min_mz: f64,
    pub max_mz: f64,
    pub step: f64,
}

impl Default for MzScanGrid {
    fn default() -> Self {
        Self {
            min_mz: 40.0,
            max_mz: 1000.0,
            step: 0.006,
        }
    }
}

#[derive(Clone, Debug)]
pub struct FindFeaturesOptions {
    pub seed_eic_options: EicOptions,
    pub final_eic_options: EicOptions,
    pub peak_options: FindPeaksOptions,
    pub mz_scan_grid: MzScanGrid,
    pub min_seed_width_points: usize,
    pub memory_budget_bytes: usize,
    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    pub gpu_context: Option<Arc<GpuContext>>,
    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    pub use_gpu: bool,
    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    pub batch_size: Option<usize>,
}

impl Default for FindFeaturesOptions {
    fn default() -> Self {
        Self {
            seed_eic_options: EicOptions {
                ppm_tolerance: 10.0,
                mz_tolerance: 0.003,
                ..Default::default()
            },
            final_eic_options: EicOptions {
                ppm_tolerance: 20.0,
                mz_tolerance: 0.005,
                ..Default::default()
            },
            peak_options: FindPeaksOptions::default(),
            mz_scan_grid: MzScanGrid::default(),
            min_seed_width_points: 5,
            memory_budget_bytes: 256 * 1024 * 1024,
            #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
            gpu_context: None,
            #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
            use_gpu: false,
            #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
            batch_size: None,
        }
    }
}

pub fn find_features(
    reader: &mut EicReader,
    time_window: FromTo,
    options: Option<FindFeaturesOptions>,
    cores: usize,
) -> Result<Vec<Feature>, FeatureError> {
    let opts = options.unwrap_or_default();

    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    let gpu_context = if opts.use_gpu {
        opts.gpu_context.clone().or_else(|| {
            let context = GpuContext::try_init().map(Arc::new);
            if context.is_none() {
                eprintln!(
                    "[find_features] GPU requested but initialization failed, using CPU"
                );
            }
            context
        })
    } else {
        None
    };

    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    let gpu_batch = GpuBatchOptions {
        batch_size: opts.batch_size.unwrap_or_else(|| GpuBatchOptions::default().batch_size),
        vram_override: None,
    };

    if !opts.mz_scan_grid.step.is_finite() || opts.mz_scan_grid.step <= 0.0 {
        return Err(FeatureError::InvalidStepSize(opts.mz_scan_grid.step));
    }

    let grid = build_mz_grid(
        opts.mz_scan_grid.min_mz,
        opts.mz_scan_grid.max_mz,
        opts.mz_scan_grid.step,
    );
    if grid.is_empty() {
        return Err(FeatureError::EmptyGrid);
    }
    if grid.len() > 2_000_000 {
        return Err(FeatureError::GridTooLarge(grid.len()));
    }

    let scan_times = get_scan_times(reader, time_window.from, time_window.to, MS1_LEVEL);
    if scan_times.is_empty() {
        return Err(FeatureError::NoScans);
    }

    let time: Vec<f64> = scan_times.iter().map(|scan| scan.rt).collect();

    let scans = read_grid_scans(reader, &scan_times, &grid, &opts)?;

    run_with_cores(cores, || {
        #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
        let gpu = gpu_context.as_deref().map(|context| (context, &gpu_batch));
        #[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
        let gpu: GpuArgs = None;

        let candidates = collect_candidates(&scans, &grid, &time, &opts, gpu);
        if candidates.is_empty() {
            return Err(FeatureError::NoCandidateMasses);
        }

        let masses = deduplicate_masses(candidates, opts.final_eic_options);
        if masses.len() > 100_000 {
            return Err(FeatureError::TooManyCandidateMasses(masses.len()));
        }

        let features = extract_features(&scans, &masses, &time, &opts);
        if features.is_empty() {
            return Err(FeatureError::NoCandidateMasses);
        }

        let features = deduplicate_features(features, opts.final_eic_options, 0.05);
        Ok(sort_features(features))
    })
}

type Scans = Vec<(Vec<f64>, Vec<f64>)>;

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
type GpuArgs<'a> = Option<(&'a GpuContext, &'a GpuBatchOptions)>;

#[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
type GpuArgs<'a> = Option<std::convert::Infallible>;

fn grid_read_window(grid: &[f64], opts: &FindFeaturesOptions) -> (f64, f64) {
    let first_mz = grid[0];
    let last_mz = grid[grid.len() - 1];
    let low = first_mz - mz_tolerance_for(first_mz, opts.seed_eic_options);
    let high = last_mz + mz_tolerance_for(last_mz, opts.seed_eic_options);
    (low, high)
}

fn read_grid_scans(
    reader: &mut EicReader,
    scan_times: &[crate::utilities::calculate_eic::ScanTime],
    grid: &[f64],
    opts: &FindFeaturesOptions,
) -> Result<Scans, FeatureError> {
    let (mz_low, mz_high) = grid_read_window(grid, opts);
    let mut scans = Vec::with_capacity(scan_times.len());

    for scan_time in scan_times {
        let mut mz = Vec::new();
        let mut intensity = Vec::new();
        read_mz_window(reader, scan_time.index, mz_low, mz_high, &mut mz, &mut intensity)?;
        scans.push((mz, intensity));
    }

    Ok(scans)
}

fn seed_intensity_threshold(opts: &FindFeaturesOptions) -> f64 {
    opts.peak_options
        .filter
        .as_ref()
        .and_then(|filter| filter.min_intensity)
        .unwrap_or(0.0)
}

fn eic_row_for_mass(scans: &Scans, target_mz: f64, options: EicOptions) -> Vec<f64> {
    let tolerance = mz_tolerance_for(target_mz, options);
    let mz_low = target_mz - tolerance;
    let mz_high = target_mz + tolerance;
    scans
        .iter()
        .map(|(mz, intensity)| summed_intensity_in_window(mz, intensity, mz_low, mz_high))
        .collect()
}

fn row_has_peak(
    eic_row: &[f64],
    time: &[f64],
    intensity_threshold: f64,
    coarse_opts: &FindPeaksOptions,
) -> bool {
    let max_intensity = eic_row.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if max_intensity <= intensity_threshold {
        return false;
    }
    let data = DataXY {
        x: time.to_vec(),
        y: eic_row.to_vec(),
    };
    find_peaks(&data, Some(coarse_opts.clone()))
        .iter()
        .any(|peak| peak.intensity >= intensity_threshold)
}

fn collect_candidates(
    scans: &Scans,
    grid: &[f64],
    time: &[f64],
    opts: &FindFeaturesOptions,
    gpu: GpuArgs<'_>,
) -> Vec<f64> {
    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    if let Some((context, batch)) = gpu {
        return collect_candidates_gpu(context, batch, scans, grid, time, opts);
    }

    #[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
    let _ = gpu;

    collect_candidates_cpu(scans, grid, time, opts)
}

fn collect_candidates_cpu(
    scans: &Scans,
    grid: &[f64],
    time: &[f64],
    opts: &FindFeaturesOptions,
) -> Vec<f64> {
    let coarse_opts = build_coarse_opts(opts);
    let intensity_threshold = seed_intensity_threshold(opts);

    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    let candidate_iter = grid.par_iter();
    #[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
    let candidate_iter = grid.iter();

    candidate_iter
        .filter_map(|&target_mz| {
            let eic_row = eic_row_for_mass(scans, target_mz, opts.seed_eic_options);
            if row_has_peak(&eic_row, time, intensity_threshold, &coarse_opts) {
                Some(target_mz)
            } else {
                None
            }
        })
        .collect()
}

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
fn collect_candidates_gpu(
    context: &GpuContext,
    batch: &GpuBatchOptions,
    scans: &Scans,
    grid: &[f64],
    time: &[f64],
    opts: &FindFeaturesOptions,
) -> Vec<f64> {
    let flattened = FlattenedScans::from_windows(scans);
    let mut processor = GpuGridProcessor::new(context, batch.clone());

    let survivors = match processor.process(&flattened, grid, opts) {
        Ok(masses) => masses,
        Err(error) => {
            eprintln!("[find_features] GPU failed: {error}, using CPU");
            return collect_candidates_cpu(scans, grid, time, opts);
        }
    };

    let coarse_opts = build_coarse_opts(opts);
    let intensity_threshold = seed_intensity_threshold(opts);

    survivors
        .into_par_iter()
        .filter(|&target_mz| {
            let eic_row = eic_row_for_mass(scans, target_mz, opts.seed_eic_options);
            row_has_peak(&eic_row, time, intensity_threshold, &coarse_opts)
        })
        .collect()
}

fn extract_features(
    scans: &Scans,
    masses: &[f64],
    time: &[f64],
    opts: &FindFeaturesOptions,
) -> Vec<Feature> {
    let min_width = opts
        .peak_options
        .filter
        .as_ref()
        .and_then(|filter| filter.min_peak_width_points)
        .unwrap_or(opts.min_seed_width_points);

    let extract_one = |&target_mz: &f64| -> Vec<Feature> {
        let eic = eic_row_for_mass(scans, target_mz, opts.final_eic_options);
        let data = DataXY {
            x: time.to_vec(),
            y: eic,
        };
        find_peaks(&data, Some(opts.peak_options.clone()))
            .into_iter()
            .filter(|peak| min_width == 0 || peak.n_points >= min_width)
            .map(|peak| {
                let peak = with_eic_apex_intensity(&data.x, &data.y, peak);
                Feature {
                    mz: target_mz,
                    rt: peak.rt,
                    intensity: peak.intensity,
                    from: peak.from,
                    to: peak.to,
                    n_points: peak.n_points,
                    integral: peak.integral,
                    noise: peak.noise,
                }
            })
            .collect()
    };

    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    {
        masses.par_iter().flat_map_iter(extract_one).collect()
    }
    #[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
    {
        masses.iter().flat_map(extract_one).collect()
    }
}

fn build_coarse_opts(config: &FindFeaturesOptions) -> FindPeaksOptions {
    let mut coarse = config.peak_options.clone();
    let mut filter = coarse.filter.unwrap_or_default();
    filter.min_peak_width_points = Some(config.min_seed_width_points);
    coarse.filter = Some(filter);
    coarse
}

pub(crate) fn deduplicate_masses(mut masses: Vec<f64>, opts: EicOptions) -> Vec<f64> {
    if masses.is_empty() {
        return masses;
    }
    masses.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));

    let tolerance = MzTolerance {
        mz_absolute: opts.mz_tolerance.max(0.0),
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

pub(crate) fn max_intensity_centroid(
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
        mz_absolute: opts.mz_tolerance,
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
