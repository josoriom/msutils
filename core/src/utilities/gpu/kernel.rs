use std::mem::size_of;
use std::sync::mpsc;

use bytemuck::{Pod, Zeroable, bytes_of, cast_slice};
use wgpu::util::DeviceExt;

use crate::utilities::gpu::{
    context::GpuContext,
    processor::{FlattenedScans, GpuError},
};

const SHADER: &str = include_str!("shader.wgsl");
const FILTER_SHADER: &str = include_str!("filter_shader.wgsl");
const WORKGROUP_X: u32 = 64;

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct Params {
    num_targets: u32,
    num_scans: u32,
    ppm_tol: f32,
    mz_tol: f32,
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct FilterParams {
    num_targets: u32,
    num_scans: u32,
    threshold: f32,
}

struct ScanBuffers {
    mz_buf: wgpu::Buffer,
    it_buf: wgpu::Buffer,
    off_buf: wgpu::Buffer,
    len_buf: wgpu::Buffer,
}

pub(crate) struct EicKernel {
    pipeline: wgpu::ComputePipeline,
    bind_group_layout: wgpu::BindGroupLayout,
    filter_pipeline: wgpu::ComputePipeline,
    filter_bind_group_layout: wgpu::BindGroupLayout,
    scan_buffers: Option<ScanBuffers>,
    pub(crate) scan_count: u32,
    pub(crate) scan_vram_bytes: u64,
}

impl EicKernel {
    pub(crate) fn new(ctx: &GpuContext) -> Self {
        let shader = ctx
            .device
            .create_shader_module(wgpu::ShaderModuleDescriptor {
                label: None,
                source: wgpu::ShaderSource::Wgsl(SHADER.into()),
            });

        let filter_shader = ctx
            .device
            .create_shader_module(wgpu::ShaderModuleDescriptor {
                label: None,
                source: wgpu::ShaderSource::Wgsl(FILTER_SHADER.into()),
            });

        let bind_group_layout =
            ctx.device
                .create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                    label: None,
                    entries: &[
                        uniform_entry(0),
                        storage_entry(1, true),
                        storage_entry(2, true),
                        storage_entry(3, true),
                        storage_entry(4, true),
                        storage_entry(5, true),
                        storage_entry(6, false),
                    ],
                });

        let filter_bind_group_layout =
            ctx.device
                .create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                    label: None,
                    entries: &[
                        uniform_entry(0),
                        storage_entry(1, true),
                        storage_entry(2, false),
                        storage_entry(3, false),
                    ],
                });

        let pipeline_layout = ctx
            .device
            .create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                label: None,
                bind_group_layouts: &[Some(&bind_group_layout)],
                immediate_size: 0,
            });

        let filter_pipeline_layout =
            ctx.device
                .create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                    label: None,
                    bind_group_layouts: &[Some(&filter_bind_group_layout)],
                    immediate_size: 0,
                });

        let pipeline = ctx
            .device
            .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: None,
                layout: Some(&pipeline_layout),
                module: &shader,
                entry_point: Some("main"),
                compilation_options: wgpu::PipelineCompilationOptions::default(),
                cache: None,
            });

        let filter_pipeline =
            ctx.device
                .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                    label: None,
                    layout: Some(&filter_pipeline_layout),
                    module: &filter_shader,
                    entry_point: Some("main"),
                    compilation_options: wgpu::PipelineCompilationOptions::default(),
                    cache: None,
                });

        Self {
            pipeline,
            bind_group_layout,
            filter_pipeline,
            filter_bind_group_layout,
            scan_buffers: None,
            scan_count: 0,
            scan_vram_bytes: 0,
        }
    }

    pub(crate) fn load_scans(
        &mut self,
        ctx: &GpuContext,
        scans: &FlattenedScans,
        scan_count: usize,
    ) {
        self.scan_vram_bytes = (scans.mz.len() * size_of::<f32>()
            + scans.intensity.len() * size_of::<f32>()
            + scans.offsets.len() * size_of::<u32>()
            + scans.lengths.len() * size_of::<u32>()) as u64;
        self.scan_count = scan_count as u32;
        self.scan_buffers = Some(ScanBuffers {
            mz_buf: storage_buf(&ctx.device, cast_slice(&scans.mz)),
            it_buf: storage_buf(&ctx.device, cast_slice(&scans.intensity)),
            off_buf: storage_buf(&ctx.device, cast_slice(&scans.offsets)),
            len_buf: storage_buf(&ctx.device, cast_slice(&scans.lengths)),
        });
    }

    pub(crate) fn unload_scans(&mut self) {
        self.scan_buffers = None;
        self.scan_vram_bytes = 0;
        self.scan_count = 0;
    }

    pub(crate) fn dispatch(
        &self,
        ctx: &GpuContext,
        targets: &[f32],
        ppm_tol: f32,
        mz_tol: f32,
        threshold: f32,
    ) -> Result<Vec<u32>, GpuError> {
        let bufs = self
            .scan_buffers
            .as_ref()
            .ok_or(GpuError::ScanBuffersNotLoaded)?;

        let num_targets = targets.len() as u32;
        let num_scans = self.scan_count;
        let output_bytes = (targets.len() * num_scans as usize * size_of::<f32>()) as u64;

        let max_buffer = ctx.device.limits().max_buffer_size;
        if output_bytes > max_buffer {
            return Err(GpuError::BufferOverflow {
                requested: output_bytes,
                limit: max_buffer,
            });
        }

        let max_wg = ctx.device.limits().max_compute_workgroups_per_dimension;
        let wg_x = num_targets.div_ceil(WORKGROUP_X);
        let wg_y = num_scans;
        if wg_x > max_wg || wg_y > max_wg {
            return Err(GpuError::DispatchOverflow {
                x: wg_x,
                y: wg_y,
                limit: max_wg,
            });
        }

        let params = Params {
            num_targets,
            num_scans,
            ppm_tol,
            mz_tol,
        };
        let filter_params = FilterParams {
            num_targets,
            num_scans,
            threshold,
        };

        let params_buf = uniform_buf(&ctx.device, bytes_of(&params));
        let filter_params_buf = uniform_buf(&ctx.device, bytes_of(&filter_params));
        let tgt_buf = storage_buf(&ctx.device, cast_slice(targets));

        let out_buf = ctx.device.create_buffer(&wgpu::BufferDescriptor {
            label: None,
            size: output_bytes,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });

        let survivor_index_bytes = (targets.len() * size_of::<u32>()) as u64;

        let survivor_count_buf = ctx
            .device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: None,
                contents: &0u32.to_le_bytes(),
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            });

        let survivor_indices_buf = ctx.device.create_buffer(&wgpu::BufferDescriptor {
            label: None,
            size: survivor_index_bytes,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });

        let staging_count = ctx.device.create_buffer(&wgpu::BufferDescriptor {
            label: None,
            size: 4,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let staging_indices = ctx.device.create_buffer(&wgpu::BufferDescriptor {
            label: None,
            size: survivor_index_bytes,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let eic_bind_group = ctx.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: None,
            layout: &self.bind_group_layout,
            entries: &[
                bind_entry(0, &params_buf),
                bind_entry(1, &bufs.mz_buf),
                bind_entry(2, &bufs.it_buf),
                bind_entry(3, &bufs.off_buf),
                bind_entry(4, &bufs.len_buf),
                bind_entry(5, &tgt_buf),
                bind_entry(6, &out_buf),
            ],
        });

        let filter_bind_group = ctx.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: None,
            layout: &self.filter_bind_group_layout,
            entries: &[
                bind_entry(0, &filter_params_buf),
                bind_entry(1, &out_buf),
                bind_entry(2, &survivor_count_buf),
                bind_entry(3, &survivor_indices_buf),
            ],
        });

        let mut encoder = ctx
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });

        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: None,
                timestamp_writes: None,
            });
            pass.set_pipeline(&self.pipeline);
            pass.set_bind_group(0, &eic_bind_group, &[]);
            pass.dispatch_workgroups(wg_x, wg_y, 1);
        }

        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: None,
                timestamp_writes: None,
            });
            pass.set_pipeline(&self.filter_pipeline);
            pass.set_bind_group(0, &filter_bind_group, &[]);
            pass.dispatch_workgroups(num_targets.div_ceil(WORKGROUP_X), 1, 1);
        }

        encoder.copy_buffer_to_buffer(&survivor_count_buf, 0, &staging_count, 0, 4);
        encoder.copy_buffer_to_buffer(
            &survivor_indices_buf,
            0,
            &staging_indices,
            0,
            survivor_index_bytes,
        );

        let _index = ctx.queue.submit(std::iter::once(encoder.finish()));

        let count_slice = staging_count.slice(..);
        let idx_slice = staging_indices.slice(..);

        let (tx1, rx1) = mpsc::channel();
        let (tx2, rx2) = mpsc::channel();
        count_slice.map_async(wgpu::MapMode::Read, move |r| {
            let _ = tx1.send(r);
        });
        idx_slice.map_async(wgpu::MapMode::Read, move |r| {
            let _ = tx2.send(r);
        });

        ctx.device
            .poll(wgpu::PollType::wait_indefinitely())
            .unwrap();
        rx1.recv().unwrap().map_err(|_| GpuError::MapFailed)?;
        rx2.recv().unwrap().map_err(|_| GpuError::MapFailed)?;

        let survivor_count: u32 = {
            let mapped = count_slice
                .get_mapped_range()
                .map_err(|_| GpuError::MapFailed)?;
            let c = u32::from_le_bytes([mapped[0], mapped[1], mapped[2], mapped[3]]);
            drop(mapped);
            c
        };
        staging_count.unmap();

        let result: Vec<u32> = {
            let mapped = idx_slice
                .get_mapped_range()
                .map_err(|_| GpuError::MapFailed)?;
            let rows = if survivor_count == 0 {
                Vec::new()
            } else {
                let byte_end = survivor_count as usize * size_of::<u32>();
                let indices: &[u32] = cast_slice(&mapped[..byte_end]);
                indices.to_vec()
            };
            drop(mapped);
            rows
        };
        staging_indices.unmap();

        Ok(result)
    }
}

fn uniform_buf(device: &wgpu::Device, contents: &[u8]) -> wgpu::Buffer {
    device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: None,
        contents,
        usage: wgpu::BufferUsages::UNIFORM,
    })
}

fn storage_buf(device: &wgpu::Device, contents: &[u8]) -> wgpu::Buffer {
    device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: None,
        contents,
        usage: wgpu::BufferUsages::STORAGE,
    })
}

fn uniform_entry(binding: u32) -> wgpu::BindGroupLayoutEntry {
    wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Uniform,
            has_dynamic_offset: false,
            min_binding_size: None,
        },
        count: None,
    }
}

fn storage_entry(binding: u32, read_only: bool) -> wgpu::BindGroupLayoutEntry {
    wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Storage { read_only },
            has_dynamic_offset: false,
            min_binding_size: None,
        },
        count: None,
    }
}

fn bind_entry(binding: u32, buffer: &'_ wgpu::Buffer) -> wgpu::BindGroupEntry<'_> {
    wgpu::BindGroupEntry {
        binding,
        resource: buffer.as_entire_binding(),
    }
}
