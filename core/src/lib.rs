mod ffi;
pub mod utilities;

pub use utilities::structs::{DataXY, FromTo, Peak, Roi};

pub mod eic {
    pub use crate::utilities::calculate_eic::{
        EicOptions, EicReader, FastError, ScanQuery, TimeUnit, calculate_eic, get_scans,
        plan_eic_ranges,
    };
}

pub mod peaks {
    pub use crate::utilities::find_peaks::{
        ArtifactFilter, FindPeaksOptions, PeakFilter, find_peaks,
    };
    pub use crate::utilities::fit_peak::{
        PeakParameters, PeakSeed, PeakShape, draw_peak, fit_peak,
    };
    pub use crate::utilities::get_peak::get_peak;
    pub use crate::utilities::get_peaks_from_chrom::{ChromPeakRow, get_peaks_from_chrom};
    pub use crate::utilities::get_peaks_from_eic::{get_peaks_from_eic, plan_peaks_ranges};
}

pub mod boundaries {
    pub use crate::utilities::get_boundaries::{
        Boundaries, BoundariesOptions, Boundary, BoundaryMethod, get_boundaries,
    };
}

pub mod noise {
    pub use crate::utilities::find_noise_level::{find_noise_level, find_noise_level_san_plot};
}

pub mod baseline {
    pub use crate::utilities::calculate_baseline::{BaselineOptions, calculate_baseline};
}

pub mod features {
    pub use crate::utilities::find_feature::{FindFeatureOptions, find_feature};
    pub use crate::utilities::find_features::{FindFeaturesOptions, find_features};
    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    pub use crate::utilities::get_features::{AlignmentOptions, get_features};
    pub use crate::utilities::mz_estimator::MzEstimatorKind;
}
