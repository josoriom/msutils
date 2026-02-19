pub mod calculate_eic;
pub use calculate_eic::{Eic, EicOptions, calculate_eic};

pub mod find_features;
pub use find_features::find_features;

pub mod find_noise_level;
pub use find_noise_level::find_noise_level;

pub mod find_peaks;

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

pub mod utilities;
pub use utilities::closest_index;

pub mod find_feature;
pub use find_feature::find_feature;
