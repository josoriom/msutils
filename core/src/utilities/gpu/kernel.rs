use std::sync::mpsc;

use bytemuck::{Pod, Zeroable, bytes_of, cast_slice};
use wgpu::util::DeviceExt;

use crate::utilities::gpu::{
    context::GpuContext,
    processor::{FlattenedScans, GpuError},
};

const SHADER: &str = include_str!("shader.wgsl");
const WORKGROUP: u32 = 64;

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct Params {
    num_targets: u32,
    num_scans: u32,
    ppm_tol: f32,
    mz_tol: f32,
}

pub(crate) struct EicKernel {
    pipeline: wgpu::ComputePipeline,
    bind_group_layout: wgpu::BindGroupLayout,
}

impl EicKernel {
    pub(crate) fn new(ctx: &GpuContext) -> Self {
        let shader = ctx
            .device
            .create_shader_module(wgpu::ShaderModuleDescriptor {
                label: None,
                source: wgpu::ShaderSource::Wgsl(SHADER.into()),
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

        let pipeline_layout = ctx
            .device
            .create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                label: None,
                bind_group_layouts: &[Some(&bind_group_layout)],
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

        Self {
            pipeline,
            bind_group_layout,
        }
    }

    pub(crate) fn dispatch(
        &self,
        ctx: &GpuContext,
        targets: &[f32],
        scans: &FlattenedScans,
        scan_count: usize,
        ppm_tol: f32,
        mz_tol: f32,
    ) -> Result<Vec<f32>, GpuError> {
        let num_targets = targets.len() as u32;
        let num_scans = scan_count as u32;
        let output_bytes = (targets.len() * scan_count * size_of::<f32>()) as u64;

        let max_buffer = ctx.device.limits().max_buffer_size;
        if output_bytes > max_buffer {
            return Err(GpuError::BufferOverflow {
                requested: output_bytes,
                limit: max_buffer,
            });
        }

        let params = Params {
            num_targets,
            num_scans,
            ppm_tol,
            mz_tol,
        };

        let params_buf = uniform_buf(&ctx.device, bytes_of(&params));
        let mz_buf = storage_buf(&ctx.device, cast_slice(&scans.mz));
        let it_buf = storage_buf(&ctx.device, cast_slice(&scans.intensity));
        let off_buf = storage_buf(&ctx.device, cast_slice(&scans.offsets));
        let len_buf = storage_buf(&ctx.device, cast_slice(&scans.lengths));
        let tgt_buf = storage_buf(&ctx.device, cast_slice(targets));

        let out_buf = ctx.device.create_buffer(&wgpu::BufferDescriptor {
            label: None,
            size: output_bytes,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });

        let staging_buf = ctx.device.create_buffer(&wgpu::BufferDescriptor {
            label: None,
            size: output_bytes,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let bind_group = ctx.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: None,
            layout: &self.bind_group_layout,
            entries: &[
                bind_entry(0, &params_buf),
                bind_entry(1, &mz_buf),
                bind_entry(2, &it_buf),
                bind_entry(3, &off_buf),
                bind_entry(4, &len_buf),
                bind_entry(5, &tgt_buf),
                bind_entry(6, &out_buf),
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
            pass.set_bind_group(0, &bind_group, &[]);
            pass.dispatch_workgroups(num_targets.div_ceil(WORKGROUP), 1, 1);
        }

        encoder.copy_buffer_to_buffer(&out_buf, 0, &staging_buf, 0, output_bytes);
        let _index = ctx.queue.submit(std::iter::once(encoder.finish()));

        let slice = staging_buf.slice(..);
        let (tx, rx) = mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |r| {
            let _ = tx.send(r);
        });
        ctx.device
            .poll(wgpu::PollType::wait_indefinitely())
            .unwrap();
        rx.recv().unwrap().map_err(|_| GpuError::MapFailed)?;

        let mapped = slice.get_mapped_range();
        let result: Vec<f32> = cast_slice(&mapped).to_vec();
        drop(mapped);
        staging_buf.unmap();

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
