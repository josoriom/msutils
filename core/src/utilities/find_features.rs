use ionic::SpectrumSource;
use serde::Serialize;
use std::{
    cmp::Ordering,
    collections::HashSet,
    error::Error,
    fmt::{Display, Formatter},
    mem,
};

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use crate::utilities::gpu::{GpuContext, processor::GpuBatchOptions};
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
    pub scan_eic_options: EicOptions,
    pub eic_options: EicOptions,
    pub find_peaks: FindPeaksOptions,
    pub mz_scan_grid: MzScanGrid,
    pub scan_width_threshold: usize,
    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    pub use_gpu: bool,
    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    pub batch_size: Option<usize>,
}

impl Default for FindFeaturesOptions {
    fn default() -> Self {
        Self {
            scan_eic_options: EicOptions {
                ppm_tolerance: 10.0,
                mz_tolerance: 0.003,
                ..Default::default()
            },
            eic_options: EicOptions {
                ppm_tolerance: 20.0,
                mz_tolerance: 0.005,
                ..Default::default()
            },
            find_peaks: FindPeaksOptions::default(),
            mz_scan_grid: MzScanGrid::default(),
            scan_width_threshold: 5,
            #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
            use_gpu: false,
            #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
            batch_size: None,
        }
    }
}

pub fn find_features(
    source: &mut impl SpectrumSource,
    time_window: FromTo,
    options: Option<FindFeaturesOptions>,
    cores: usize,
) -> Result<Vec<Feature>, FeatureError> {
    let opts = options.unwrap_or_default();

    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    let gpu_ctx = if opts.use_gpu {
        let ctx = GpuContext::try_init();
        if ctx.is_none() {
            eprintln!(
                "[find_features] GPU requested but initialization failed, falling back to CPU"
            );
        } else {
            eprintln!("[find_features] GPU initialized successfully");
        }
        ctx
    } else {
        eprintln!("[find_features] GPU disabled, using CPU");
        None
    };

    if !opts.mz_scan_grid.step_size.is_finite() || opts.mz_scan_grid.step_size <= 0.0 {
        return Err(FeatureError::InvalidStepSize(opts.mz_scan_grid.step_size));
    }

    let grid = build_mz_grid(
        opts.mz_scan_grid.mz_min,
        opts.mz_scan_grid.mz_max,
        opts.mz_scan_grid.step_size,
    );
    if grid.is_empty() {
        return Err(FeatureError::EmptyGrid);
    }
    if grid.len() > 2_000_000 {
        return Err(FeatureError::GridTooLarge(grid.len()));
    }

    let (time, scans) = collect_scans(source, time_window, TimeUnit::Minutes, 1);
    if scans.is_empty() {
        return Err(FeatureError::NoScans);
    }

    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    let gpu_batch = gpu_ctx.as_ref().map(|ctx| {
        safe_batch_options(
            ctx,
            time.len(),
            flatten_scans_size_estimate(&scans),
            opts.batch_size,
        )
    });

    run_with_cores(cores, || {
        #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
        let gpu_args: GpuArgs<'_> = gpu_ctx
            .as_ref()
            .zip(gpu_batch.as_ref())
            .map(|(ctx, batch)| (ctx, batch));

        #[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
        let gpu_args: GpuArgs<'_> = None;

        let masses = collect_candidate_masses(&scans, &time, &grid, &opts, gpu_args)?;

        let masses = deduplicate_masses(masses, opts.eic_options);

        let features = extract_features_from_masses(&masses, &scans, &time, &opts);

        let features = deduplicate_features(features, opts.eic_options, 0.05);

        Ok(sort_features(features))
    })
}

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
fn run_with_cores<F, T>(cores: usize, f: F) -> T
where
    F: FnOnce() -> T + Send,
    T: Send,
{
    if cores <= 1 {
        return f();
    }

    ThreadPoolBuilder::new()
        .num_threads(cores.max(1))
        .thread_name(|i| format!("ff-{}", i))
        .build()
        .expect("failed to build rayon pool")
        .install(f)
}

#[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
fn run_with_cores<F, T>(_cores: usize, f: F) -> T
where
    F: FnOnce() -> T,
{
    f()
}

fn build_coarse_opts(config: &FindFeaturesOptions) -> FindPeaksOptions {
    let mut coarse = config.find_peaks.clone();
    let mut filter = coarse.filter_peaks_options.unwrap_or_default();
    filter.width_threshold = Some(config.scan_width_threshold);
    coarse.filter_peaks_options = Some(filter);
    coarse
}

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
type GpuArgs<'a> = Option<(&'a GpuContext, &'a GpuBatchOptions)>;

#[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
type GpuArgs<'a> = Option<std::convert::Infallible>;

fn collect_candidate_masses(
    scans: &[CentroidScan],
    time: &[f64],
    grid: &[f64],
    config: &FindFeaturesOptions,
    gpu: GpuArgs<'_>,
) -> Result<Vec<f64>, FeatureError> {
    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    if let Some((ctx, batch_opts)) = gpu {
        return collect_candidates_gpu(ctx, batch_opts, scans, time, grid, config);
    }
    collect_candidates_cpu(scans, time, grid, config)
}

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
fn collect_candidates_gpu(
    ctx: &GpuContext,
    batch_opts: &GpuBatchOptions,
    scans: &[CentroidScan],
    time: &[f64],
    grid: &[f64],
    config: &FindFeaturesOptions,
) -> Result<Vec<f64>, FeatureError> {
    use crate::utilities::gpu::GpuGridProcessor;

    let coarse_opts = build_coarse_opts(config);
    let intensity_threshold = config
        .find_peaks
        .filter_peaks_options
        .as_ref()
        .and_then(|o| o.intensity_threshold)
        .unwrap_or(0.0);

    let mut processor = GpuGridProcessor::new(ctx, batch_opts.clone());
    match processor.process(scans, grid, config) {
        Ok(survivors) if !survivors.is_empty() => {
            let masses: Vec<f64> = survivors
                .par_iter()
                .filter_map(|&mz| {
                    let eic = compute_eic_for_mz(scans, time.len(), mz, config.scan_eic_options);
                    evaluate_eic_row(
                        &eic,
                        time,
                        scans,
                        mz,
                        intensity_threshold,
                        &coarse_opts,
                        config.scan_eic_options,
                    )
                })
                .collect();
            Ok(masses)
        }
        Ok(_) => {
            eprintln!("[gpu] falling back to CPU");
            collect_candidates_cpu(scans, time, grid, config)
        }
        Err(e) => {
            eprintln!("[gpu] failed: {e}, falling back to CPU");
            collect_candidates_cpu(scans, time, grid, config)
        }
    }
}

fn collect_candidates_cpu(
    scans: &[CentroidScan],
    time: &[f64],
    grid: &[f64],
    config: &FindFeaturesOptions,
) -> Result<Vec<f64>, FeatureError> {
    let coarse_opts = build_coarse_opts(config);
    let intensity_threshold = config
        .find_peaks
        .filter_peaks_options
        .as_ref()
        .and_then(|o| o.intensity_threshold)
        .unwrap_or(0.0);

    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    let masses: Vec<f64> = grid
        .par_iter()
        .filter_map(|&mz| {
            let eic = compute_eic_for_mz(scans, time.len(), mz, config.scan_eic_options);
            evaluate_eic_row(
                &eic,
                time,
                scans,
                mz,
                intensity_threshold,
                &coarse_opts,
                config.scan_eic_options,
            )
        })
        .collect();

    #[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
    let masses: Vec<f64> = grid
        .iter()
        .filter_map(|&mz| {
            let eic = compute_eic_for_mz(scans, time.len(), mz, config.scan_eic_options);
            evaluate_eic_row(
                &eic,
                time,
                scans,
                mz,
                intensity_threshold,
                &coarse_opts,
                config.scan_eic_options,
            )
        })
        .collect();

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
    config: &FindFeaturesOptions,
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
    config: &FindFeaturesOptions,
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
    let features: Vec<Feature> = masses
        .par_iter()
        .map(|&mz| extract_peaks_for_mass(mz, scans, time, config, final_w))
        .reduce(Vec::new, |mut a, mut b| {
            if !b.is_empty() {
                a.append(&mut b);
            }
            a
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

pub(crate) fn evaluate_eic_row(
    eic: &[f64],
    time: &[f64],
    scans: &[CentroidScan],
    mz_target: f64,
    intensity_threshold: f64,
    coarse_opts: &FindPeaksOptions,
    scan_eic_options: EicOptions,
) -> Option<f64> {
    let eic_f32: Vec<f64> = eic.iter().map(|&v| v as f32 as f64).collect();
    if eic_f32.iter().copied().fold(f64::NEG_INFINITY, f64::max) <= intensity_threshold {
        return None;
    }
    let data = DataXY {
        x: time.to_vec(),
        y: eic_f32,
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

fn flatten_scans_size_estimate(scans: &[CentroidScan]) -> u64 {
    let total: usize = scans.iter().map(|s| s.mz.len()).sum();
    ((total * 2 * size_of::<f32>()) + (scans.len() * 2 * size_of::<u32>())) as u64
}

pub fn safe_batch_options(
    ctx: &GpuContext,
    scan_count: usize,
    scan_vram_bytes: u64,
    requested: Option<usize>,
) -> GpuBatchOptions {
    let vram = ctx.vram_bytes.saturating_sub(scan_vram_bytes);
    let bytes_per_row = scan_count * size_of::<f32>();
    let max_from_vram = ((vram as f64 * 0.8) as usize / bytes_per_row).max(1);
    let max_from_buf = (ctx.device.limits().max_buffer_size as usize / bytes_per_row).max(1);
    let max_from_bind =
        (ctx.device.limits().max_storage_buffer_binding_size as usize / bytes_per_row).max(1);
    let max_safe = max_from_vram.min(max_from_buf).min(max_from_bind);
    let batch_size = requested.unwrap_or(max_safe);
    let clamped = batch_size.min(max_safe).max(1);

    if clamped < batch_size {
        eprintln!(
            "[gpu] batch_size clamped {} → {} (buf limit {} MB, bind limit {} MB, vram limit {} MB)",
            batch_size,
            clamped,
            ctx.device.limits().max_buffer_size / (1024 * 1024),
            ctx.device.limits().max_storage_buffer_binding_size / (1024 * 1024),
            vram / (1024 * 1024),
        );
    }

    GpuBatchOptions {
        batch_size: clamped,
        vram_override: None,
    }
}
