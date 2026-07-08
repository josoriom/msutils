use std::panic::catch_unwind;
use std::sync::OnceLock;
use wgpu;

#[derive(Clone, Debug)]
pub struct GpuContext {
    pub(crate) device: wgpu::Device,
    pub(crate) queue: wgpu::Queue,
    pub(crate) vram_bytes: u64,
}

impl GpuContext {
    pub fn try_init() -> Option<Self> {
        static GPU_CONTEXT: OnceLock<Option<GpuContext>> = OnceLock::new();

        GPU_CONTEXT
            .get_or_init(|| {
                catch_unwind(|| pollster::block_on(Self::init_async())).unwrap_or_else(|_| {
                    eprintln!(
                        "[gpu][WARN] GPU init panicked (no compatible async executor), falling back to CPU"
                    );
                    None
                })
            })
            .clone()
    }

    async fn init_async() -> Option<Self> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends: wgpu::Backends::all(),
            flags: wgpu::InstanceFlags::default(),
            memory_budget_thresholds: wgpu::MemoryBudgetThresholds::default(),
            backend_options: wgpu::BackendOptions::default(),
            display: None,
        });

        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                compatible_surface: None,
                force_fallback_adapter: false,
                apply_limit_buckets: false,
            })
            .await
            .ok()?;

        let info = adapter.get_info();
        eprintln!("[gpu] adapter: {} ({:?})", info.name, info.device_type);

        let (device, queue) = adapter
            .request_device(&wgpu::DeviceDescriptor {
                label: None,
                required_features: wgpu::Features::empty(),
                required_limits: adapter.limits(),
                memory_hints: wgpu::MemoryHints::Performance,
                ..Default::default()
            })
            .await
            .ok()?;
        let vram_bytes = vram_estimate(&adapter);
        eprintln!(
            "[gpu] initialized, estimated VRAM budget: {} MB",
            vram_bytes / (1024 * 1024)
        );

        Some(Self {
            vram_bytes,
            device,
            queue,
        })
    }
}

fn vram_estimate(adapter: &wgpu::Adapter) -> u64 {
    const MB: u64 = 1024 * 1024;
    match adapter.get_info().device_type {
        wgpu::DeviceType::DiscreteGpu => 4096 * MB,
        wgpu::DeviceType::IntegratedGpu => 512 * MB,
        _ => 256 * MB,
    }
}
