pub mod calculate_eic;
pub mod ion_image;
pub use calculate_eic::{
    EicOptions, EicReader, FastError, ScanTime, calculate_eic, get_scan_times, plan_eic_ranges,
    read_mz_window,
};

pub mod parallel;

pub mod find_features;
pub use find_features::find_features;

pub mod find_masses;

pub mod bridge;

pub mod find_noise_level;
pub use find_noise_level::{find_noise_level, find_noise_level_san_plot};

pub mod find_peaks;

pub mod shape_filter;

pub mod fit_peak;

pub mod functions;

pub mod calculate_baseline;
pub use calculate_baseline::calculate_baseline;

pub mod get_boundaries;

pub mod get_peak;

pub mod get_peaks_from_chrom;

pub mod get_peaks_from_eic;

pub mod scan_for_peaks;

pub mod structs;

pub mod cheminfo;

pub mod math;
pub use math::closest_index;

pub mod find_feature;
pub use find_feature::find_feature;

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
pub mod get_features;
pub mod mz_estimator;
#[cfg(test)]
mod tests;
