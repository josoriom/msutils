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
        CentroidScan, EicOptions, EicReader, FastError, MS1_LEVEL, SpectrumKind, get_scan_times,
        get_spectrum_kind, lower_bound, mz_tolerance_for, read_mz_window,
        summed_intensity_in_window, with_eic_apex_intensity,
    },
    find_masses::find_masses,
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

    #[cfg(test)]
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
    pub ppm: f64,
}

impl Default for MzScanGrid {
    fn default() -> Self {
        Self {
            min_mz: 40.0,
            max_mz: 1000.0,
            step: 0.006,
            ppm: 20.0,
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

    if !opts.mz_scan_grid.step.is_finite() || opts.mz_scan_grid.step <= 0.0 {
        return Err(FeatureError::InvalidStepSize(opts.mz_scan_grid.step));
    }

    let grid = build_mz_grid(
        opts.mz_scan_grid.min_mz,
        opts.mz_scan_grid.max_mz,
        opts.mz_scan_grid.step,
        opts.mz_scan_grid.ppm,
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
    let kind = get_spectrum_kind(reader);
    let scans = read_grid_scans(reader, &scan_times, &grid, &opts)?;
    let grid = keep_masses_with_signal(
        &grid,
        &scans,
        seed_intensity_threshold(&opts),
        opts.seed_eic_options,
        cores,
    );

    detect_features(&scans, &grid, &time, kind, &opts, cores)
}

pub fn detect_features(
    scans: &[Scan],
    grid: &[f64],
    time: &[f64],
    kind: SpectrumKind,
    opts: &FindFeaturesOptions,
    cores: usize,
) -> Result<Vec<Feature>, FeatureError> {
    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    let gpu_context = if opts.use_gpu {
        opts.gpu_context.clone().or_else(|| {
            let context = GpuContext::try_init().map(Arc::new);
            if context.is_none() {
                eprintln!("[find_features] GPU requested but initialization failed, using CPU");
            }
            context
        })
    } else {
        None
    };

    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    let gpu_batch = GpuBatchOptions {
        batch_size: opts
            .batch_size
            .unwrap_or_else(|| GpuBatchOptions::default().batch_size),
        vram_override: None,
    };

    run_with_cores(cores, || {
        #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
        let gpu = gpu_context.as_deref().map(|context| (context, &gpu_batch));
        #[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
        let gpu: GpuArgs = None;

        let candidates = match kind {
            SpectrumKind::Centroid => grid.to_vec(),
            SpectrumKind::Profile => collect_candidates(scans, grid, opts, gpu),
        };
        if candidates.is_empty() {
            return Err(FeatureError::NoCandidateMasses);
        }

        let masses = deduplicate_masses(candidates, opts.final_eic_options);
        if masses.len() > 100_000 {
            return Err(FeatureError::TooManyCandidateMasses(masses.len()));
        }

        let features = extract_features(scans, &masses, time, opts, kind);
        if features.is_empty() {
            return Err(FeatureError::NoCandidateMasses);
        }

        let features = deduplicate_features(features, opts.final_eic_options, 0.05);
        Ok(sort_features(features))
    })
}

#[derive(Clone, Debug, Default)]
pub struct Scan {
    pub mz: Vec<f64>,
    pub intensity: Vec<f64>,
}

type Scans = Vec<Scan>;

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
        read_mz_window(
            reader,
            scan_time.index,
            mz_low,
            mz_high,
            &mut mz,
            &mut intensity,
        )?;
        scans.push(Scan { mz, intensity });
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

pub(crate) fn keep_masses_with_signal(
    grid: &[f64],
    scans: &[Scan],
    min_intensity: f64,
    eic_options: EicOptions,
    cores: usize,
) -> Vec<f64> {
    let low_bounds: Vec<f64> = grid
        .iter()
        .map(|&mass| mass - mz_tolerance_for(mass, eic_options))
        .collect();
    let high_bounds: Vec<f64> = grid
        .iter()
        .map(|&mass| mass + mz_tolerance_for(mass, eic_options))
        .collect();

    let has_signal = find_masses_with_signal(&low_bounds, &high_bounds, scans, min_intensity, cores);

    grid.iter()
        .copied()
        .zip(has_signal)
        .filter_map(|(mass, found)| found.then_some(mass))
        .collect()
}

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
fn find_masses_with_signal(
    low_bounds: &[f64],
    high_bounds: &[f64],
    scans: &[Scan],
    min_intensity: f64,
    cores: usize,
) -> Vec<bool> {
    let width = low_bounds.len();
    run_with_cores(cores, || {
        scans
            .par_iter()
            .fold(
                || vec![false; width],
                |mut found, scan| {
                    mark_signal_for_scan(scan, low_bounds, high_bounds, min_intensity, &mut found);
                    found
                },
            )
            .reduce(
                || vec![false; width],
                |mut base, other| {
                    for (slot, value) in base.iter_mut().zip(other) {
                        *slot |= value;
                    }
                    base
                },
            )
    })
}

#[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
fn find_masses_with_signal(
    low_bounds: &[f64],
    high_bounds: &[f64],
    scans: &[Scan],
    min_intensity: f64,
    _cores: usize,
) -> Vec<bool> {
    let mut found = vec![false; low_bounds.len()];
    for scan in scans {
        mark_signal_for_scan(scan, low_bounds, high_bounds, min_intensity, &mut found);
    }
    found
}

fn mark_signal_for_scan(
    scan: &Scan,
    low_bounds: &[f64],
    high_bounds: &[f64],
    min_intensity: f64,
    found: &mut [bool],
) {
    let masses = &scan.mz;
    let intensities = &scan.intensity;
    if masses.len() != intensities.len() {
        return;
    }
    let mut point = 0usize;
    for index in 0..found.len() {
        while point < masses.len()
            && (masses[point] < low_bounds[index] || intensities[point] < min_intensity)
        {
            point += 1;
        }
        if point >= masses.len() {
            break;
        }
        if masses[point] <= high_bounds[index] {
            found[index] = true;
        }
    }
}

fn eic_row_for_mass(scans: &[Scan], target_mz: f64, options: EicOptions) -> Vec<f64> {
    let tolerance = mz_tolerance_for(target_mz, options);
    let mz_low = target_mz - tolerance;
    let mz_high = target_mz + tolerance;
    scans
        .iter()
        .map(|scan| summed_intensity_in_window(&scan.mz, &scan.intensity, mz_low, mz_high))
        .collect()
}

pub(crate) fn max_eic_per_mass(scans: &[Scan], masses: &[f64], options: EicOptions) -> Vec<f64> {
    if masses.is_empty() {
        return Vec::new();
    }
    let window_low: Vec<f64> = masses
        .iter()
        .map(|&mass| mass - mz_tolerance_for(mass, options))
        .collect();
    let window_high: Vec<f64> = masses
        .iter()
        .map(|&mass| mass + mz_tolerance_for(mass, options))
        .collect();
    let mut max_sum = vec![0.0f64; masses.len()];
    for scan in scans {
        let points = &scan.mz;
        let intensities = &scan.intensity;
        if points.len() != intensities.len() {
            continue;
        }
        let mut start = lower_bound(points, window_low[0]);
        let mut end = start;
        for index in 0..masses.len() {
            while start < points.len() && points[start] < window_low[index] {
                start += 1;
            }
            if end < start {
                end = start;
            }
            while end < points.len() && points[end] <= window_high[index] {
                end += 1;
            }
            let window_sum: f64 = intensities[start..end].iter().sum();
            if window_sum > max_sum[index] {
                max_sum[index] = window_sum;
            }
        }
    }
    max_sum
}

fn keep_reaching(scans: &[Scan], masses: &[f64], threshold: f64, options: EicOptions) -> Vec<f64> {
    let max_value = max_eic_per_mass(scans, masses, options);
    masses
        .iter()
        .zip(max_value)
        .filter_map(|(&mass, value)| if value > threshold { Some(mass) } else { None })
        .collect()
}

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
const CHUNKS_PER_THREAD: usize = 8;

fn filter_by_max(scans: &[Scan], masses: &[f64], threshold: f64, options: EicOptions) -> Vec<f64> {
    if masses.is_empty() {
        return Vec::new();
    }
    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    {
        let threads = rayon::current_num_threads().max(1);
        let chunk = masses.len().div_ceil(threads * CHUNKS_PER_THREAD).max(1);
        masses
            .par_chunks(chunk)
            .flat_map_iter(|part| keep_reaching(scans, part, threshold, options))
            .collect()
    }
    #[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
    {
        keep_reaching(scans, masses, threshold, options)
    }
}

fn collect_candidates(
    scans: &[Scan],
    grid: &[f64],
    opts: &FindFeaturesOptions,
    gpu: GpuArgs<'_>,
) -> Vec<f64> {
    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    if let Some((context, batch)) = gpu {
        return collect_candidates_gpu(context, batch, scans, grid, opts);
    }

    #[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
    let _ = gpu;

    collect_candidates_cpu(scans, grid, opts)
}

fn collect_candidates_cpu(scans: &[Scan], grid: &[f64], opts: &FindFeaturesOptions) -> Vec<f64> {
    let threshold = seed_intensity_threshold(opts);
    filter_by_max(scans, grid, threshold, opts.seed_eic_options)
}

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
fn collect_candidates_gpu(
    context: &GpuContext,
    batch: &GpuBatchOptions,
    scans: &[Scan],
    grid: &[f64],
    opts: &FindFeaturesOptions,
) -> Vec<f64> {
    let flattened = FlattenedScans::from_windows(scans);
    let mut processor = GpuGridProcessor::new(context, batch.clone());

    let survivors = match processor.process(&flattened, grid, opts) {
        Ok(masses) => masses,
        Err(error) => {
            eprintln!("[find_features] GPU failed: {error}, using CPU");
            return collect_candidates_cpu(scans, grid, opts);
        }
    };

    let threshold = seed_intensity_threshold(opts);
    let mut survivors = survivors;
    survivors.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    filter_by_max(scans, &survivors, threshold, opts.seed_eic_options)
}

fn extract_features(
    scans: &[Scan],
    masses: &[f64],
    time: &[f64],
    opts: &FindFeaturesOptions,
    kind: SpectrumKind,
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
                let measured_mz = match apex_scan_index(&data.x, &data.y, peak.from, peak.to) {
                    Some(index) => {
                        peak_apex_mz(&scans[index], target_mz, opts.final_eic_options, kind)
                    }
                    None => target_mz,
                };
                Feature {
                    mz: measured_mz,
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

fn apex_scan_index(time: &[f64], eic: &[f64], rt_from: f64, rt_to: f64) -> Option<usize> {
    let mut apex = None;
    let mut best = 0.0;
    for (index, &rt) in time.iter().enumerate() {
        if rt < rt_from || rt > rt_to {
            continue;
        }
        if eic[index] > best {
            best = eic[index];
            apex = Some(index);
        }
    }
    apex
}

fn peak_apex_mz(apex_scan: &Scan, center_mz: f64, options: EicOptions, kind: SpectrumKind) -> f64 {
    let tolerance = mz_tolerance_for(center_mz, options);
    let mz_lo = center_mz - tolerance;
    let mz_hi = center_mz + tolerance;
    find_masses(&apex_scan.mz, &apex_scan.intensity, mz_lo, mz_hi, kind)
        .into_iter()
        .min_by(|a, b| {
            (a.mz - center_mz)
                .abs()
                .partial_cmp(&(b.mz - center_mz).abs())
                .unwrap_or(Ordering::Equal)
        })
        .map(|mass| mass.mz)
        .unwrap_or(center_mz)
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
    let mut anchor = masses[0];
    let mut sum = masses[0];
    let mut count = 1usize;

    for &m in masses.iter().skip(1) {
        if tolerance.are_close_to_ref(m, anchor) {
            sum += m;
            count += 1;
        } else {
            out.push(sum / count as f64);
            anchor = m;
            sum = m;
            count = 1;
        }
    }
    out.push(sum / count as f64);
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

pub(crate) fn build_mz_grid(start: f64, end: f64, step_da: f64, step_ppm: f64) -> Vec<f64> {
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
        let step = step_da.max(step_ppm * m * 1e-6);
        m += step;
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
