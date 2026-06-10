#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use std::collections::HashMap;
#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use std::sync::{Arc, Mutex, OnceLock};

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use rayon::{ThreadPool, ThreadPoolBuilder};

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
fn shared_pools() -> &'static Mutex<HashMap<usize, Arc<ThreadPool>>> {
    static POOLS: OnceLock<Mutex<HashMap<usize, Arc<ThreadPool>>>> = OnceLock::new();
    POOLS.get_or_init(|| Mutex::new(HashMap::new()))
}

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
fn get_pool(cores: usize) -> Arc<ThreadPool> {
    let core_count = cores.max(1);
    let mut pools = shared_pools().lock().unwrap();
    pools
        .entry(core_count)
        .or_insert_with(|| {
            Arc::new(
                ThreadPoolBuilder::new()
                    .num_threads(core_count)
                    .thread_name(|index| format!("msutils-{}", index))
                    .build()
                    .expect("failed to build rayon pool"),
            )
        })
        .clone()
}

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
pub(crate) fn run_with_cores<F, T>(cores: usize, work: F) -> T
where
    F: FnOnce() -> T + Send,
    T: Send,
{
    get_pool(cores).install(work)
}

#[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
pub(crate) fn run_with_cores<F, T>(_cores: usize, work: F) -> T
where
    F: FnOnce() -> T,
{
    work()
}
