#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use rayon::ThreadPoolBuilder;

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
pub(crate) fn run_with_cores<F, T>(cores: usize, work: F) -> T
where
    F: FnOnce() -> T + Send,
    T: Send,
{
    ThreadPoolBuilder::new()
        .num_threads(cores.max(1))
        .thread_name(|index| format!("msutils-{}", index))
        .build()
        .expect("failed to build rayon pool")
        .install(work)
}

#[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
pub(crate) fn run_with_cores<F, T>(_cores: usize, work: F) -> T
where
    F: FnOnce() -> T,
{
    work()
}
