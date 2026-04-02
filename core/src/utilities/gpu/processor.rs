#![cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]

use std::mem::size_of;

use rayon::prelude::*;

use crate::utilities::{
    calculate_eic::CentroidScan,
    find_features::{FindFeaturesOptions, evaluate_eic_row},
    find_peaks::FindPeaksOptions,
};

use super::{context::GpuContext, kernel::EicKernel};
use std::fmt::{Display, Formatter};

#[derive(Clone, Debug)]
pub struct GpuBatchOptions {
    pub batch_size: usize,
    pub vram_override: Option<u64>,
}

impl Default for GpuBatchOptions {
    fn default() -> Self {
        Self {
            batch_size: 1024,
            vram_override: None,
        }
    }
}

pub struct GpuGridProcessor<'a> {
    ctx: &'a GpuContext,
    kernel: EicKernel,
    opts: GpuBatchOptions,
}

impl<'a> GpuGridProcessor<'a> {
    pub fn new(ctx: &'a GpuContext, opts: GpuBatchOptions) -> Self {
        Self {
            kernel: EicKernel::new(ctx),
            ctx,
            opts,
        }
    }

    pub fn process(
        &self,
        scans: &[CentroidScan],
        time: &[f64],
        grid: &[f64],
        config: &FindFeaturesOptions,
    ) -> Result<Vec<f64>, GpuError> {
        let flattened = flatten_scans(scans);
        let batch_size = effective_batch_size(self.ctx, time.len(), &self.opts);
        let coarse_opts = coarse_peak_opts(config);
        let ppm_tol = config.scan_eic_options.ppm_tolerance as f32;
        let mz_tol = config.scan_eic_options.mz_tolerance as f32;

        let mut candidates = Vec::new();

        for chunk in grid.chunks(batch_size) {
            let targets: Vec<f32> = chunk.iter().map(|&v| v as f32).collect();

            let eic_matrix = self.kernel.dispatch(
                self.ctx,
                &targets,
                &flattened,
                time.len(),
                ppm_tol,
                mz_tol,
            )?;

            let batch =
                pick_peaks_from_matrix(&eic_matrix, chunk, scans, time, config, &coarse_opts);

            candidates.extend(batch);
        }

        Ok(candidates)
    }
}

fn pick_peaks_from_matrix(
    eic_matrix: &[f32],
    chunk: &[f64],
    scans: &[CentroidScan],
    time: &[f64],
    config: &FindFeaturesOptions,
    coarse_opts: &FindPeaksOptions,
) -> Vec<f64> {
    let scan_count = time.len();
    let threshold = config
        .find_peaks
        .filter_peaks_options
        .as_ref()
        .and_then(|o| o.intensity_threshold)
        .unwrap_or(0.0);

    chunk
        .par_iter()
        .enumerate()
        .filter_map(|(i, &mz_target)| {
            let row: Vec<f64> = eic_matrix[i * scan_count..(i + 1) * scan_count]
                .iter()
                .map(|&v| v as f64)
                .collect();
            evaluate_eic_row(
                &row,
                time,
                scans,
                mz_target,
                threshold,
                coarse_opts,
                config.scan_eic_options,
            )
        })
        .collect()
}

fn effective_batch_size(ctx: &GpuContext, scan_count: usize, opts: &GpuBatchOptions) -> usize {
    let vram = opts.vram_override.unwrap_or(ctx.vram_bytes);
    let bytes_per_row = scan_count * size_of::<f32>();
    let max_from_vram = ((vram as f64 * 0.8) as usize / bytes_per_row).max(1);
    let max_from_buf = (ctx.device.limits().max_buffer_size as usize / bytes_per_row).max(1);
    let max_from_bind =
        (ctx.device.limits().max_storage_buffer_binding_size as usize / bytes_per_row).max(1);
    opts.batch_size
        .min(max_from_vram)
        .min(max_from_buf)
        .min(max_from_bind)
        .max(1)
}

fn coarse_peak_opts(config: &FindFeaturesOptions) -> FindPeaksOptions {
    let mut opts = config.find_peaks.clone();
    let mut filter = opts.filter_peaks_options.unwrap_or_default();
    filter.width_threshold = Some(config.scan_width_threshold);
    opts.filter_peaks_options = Some(filter);
    opts
}

pub(crate) struct FlattenedScans {
    pub mz: Vec<f32>,
    pub intensity: Vec<f32>,
    pub offsets: Vec<u32>,
    pub lengths: Vec<u32>,
}

fn flatten_scans(scans: &[CentroidScan]) -> FlattenedScans {
    let total: usize = scans.iter().map(|s| s.mz.len()).sum();

    let mut mz = Vec::with_capacity(total);
    let mut intensity = Vec::with_capacity(total);
    let mut offsets = Vec::with_capacity(scans.len());
    let mut lengths = Vec::with_capacity(scans.len());
    let mut cursor = 0u32;

    for scan in scans {
        let n = scan.mz.len() as u32;
        offsets.push(cursor);
        lengths.push(n);
        mz.extend(scan.mz.iter().map(|&v| v as f32));
        intensity.extend(scan.intensity.iter().map(|&v| v as f32));
        cursor += n;
    }

    FlattenedScans {
        mz,
        intensity,
        offsets,
        lengths,
    }
}

#[derive(Debug)]
pub enum GpuError {
    BufferOverflow { requested: u64, limit: u64 },
    MapFailed,
}

impl Display for GpuError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::BufferOverflow { requested, limit } => {
                write!(
                    f,
                    "GPU buffer overflow: requested {} bytes, limit {} bytes",
                    requested, limit
                )
            }
            Self::MapFailed => write!(f, "GPU staging buffer mapping failed"),
        }
    }
}

impl std::error::Error for GpuError {}

pub fn safe_batch_options(
    ctx: &GpuContext,
    scan_count: usize,
    requested: Option<usize>,
) -> GpuBatchOptions {
    let vram = ctx.vram_bytes;
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
