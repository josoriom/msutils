#![cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]

use std::mem::size_of;

use crate::utilities::{calculate_eic::CentroidScan, find_features::FindFeaturesOptions};

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
        &mut self,
        scans: &[CentroidScan],
        grid: &[f64],
        config: &FindFeaturesOptions,
    ) -> Result<Vec<f64>, GpuError> {
        let flattened = flatten_scans(scans);
        let scan_count = scans.len();

        self.kernel.load_scans(self.ctx, &flattened, scan_count);

        let scan_vram = self.kernel.scan_vram_bytes;
        let batch_size = effective_batch_size(self.ctx, scan_count, &self.opts, scan_vram);

        let ppm_tol = config.seed_eic_options.ppm_tolerance as f32;
        let mz_tol = config.seed_eic_options.mz_tolerance as f32;
        let intensity_threshold = config
            .peak_options
            .filter
            .as_ref()
            .and_then(|o| o.min_intensity)
            .unwrap_or(0.0);
        let coarse_threshold = relaxed_threshold(intensity_threshold);

        let mut candidates: Vec<f64> = Vec::new();

        for chunk in grid.chunks(batch_size) {
            let targets_gpu: Vec<f32> = chunk.iter().map(|&v| v as f32).collect();
            let survivor_indices =
                self.kernel
                    .dispatch(self.ctx, &targets_gpu, ppm_tol, mz_tol, coarse_threshold)?;
            let mut sorted_survivor_indices = survivor_indices;
            sorted_survivor_indices.sort_unstable();
            candidates.extend(
                sorted_survivor_indices
                    .into_iter()
                    .map(|target_index| chunk[target_index as usize]),
            );
        }

        self.kernel.unload_scans();
        Ok(candidates)
    }
}

fn effective_batch_size(
    ctx: &GpuContext,
    scan_count: usize,
    opts: &GpuBatchOptions,
    scan_vram_bytes: u64,
) -> usize {
    let vram = opts.vram_override.unwrap_or(ctx.vram_bytes);
    let available_vram = vram.saturating_sub(scan_vram_bytes);
    let bytes_per_row = scan_count * size_of::<f32>();
    let max_from_vram = ((available_vram as f64 * 0.8) as usize / bytes_per_row).max(1);
    let max_from_buf = (ctx.device.limits().max_buffer_size as usize / bytes_per_row).max(1);
    let max_from_bind =
        (ctx.device.limits().max_storage_buffer_binding_size as usize / bytes_per_row).max(1);
    opts.batch_size
        .min(max_from_vram)
        .min(max_from_buf)
        .min(max_from_bind)
        .max(1)
}

fn relaxed_threshold(intensity_threshold: f64) -> f32 {
    if !intensity_threshold.is_finite() || intensity_threshold <= 0.0 {
        0.0
    } else {
        let slack = intensity_threshold.abs() * 1e-3 + 1e-3;
        (intensity_threshold - slack).max(0.0) as f32
    }
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
    ScanBuffersNotLoaded,
    DispatchOverflow { x: u32, y: u32, limit: u32 },
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
            Self::ScanBuffersNotLoaded => write!(f, "GPU scan buffers not loaded"),
            Self::DispatchOverflow { x, y, limit } => {
                write!(
                    f,
                    "GPU dispatch overflow: ({}, {}) exceeds limit {}",
                    x, y, limit
                )
            }
        }
    }
}

impl std::error::Error for GpuError {}
