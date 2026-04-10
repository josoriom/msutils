#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
pub(crate) mod context;
#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
pub(crate) mod kernel;
#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
pub(crate) mod processor;

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
pub use context::GpuContext;
#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
pub use processor::GpuGridProcessor;
