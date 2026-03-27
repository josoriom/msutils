use crate::utilities::calculate_eic::{
    CentroidScan, EicOptions, collect_scans, compute_eic_for_mz,
};
use crate::utilities::find_features::max_intensity_centroid;
use crate::utilities::find_peaks::FindPeaksOptions;
use crate::utilities::get_peak::get_peak;
use crate::utilities::structs::{DataXY, EicRoi, FromTo, Peak, Roi};
use octo::MzML;

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use rayon::{ThreadPoolBuilder, prelude::*};

#[derive(Clone, Debug, Default)]
pub struct Feature {
    pub id: String,
    pub mz: f64,
    pub rt: f64,
    pub peak: Peak,
}

#[derive(Clone, Debug)]
pub struct FindFeatureOptions {
    pub scan_eic_options: Option<EicOptions>,
    pub eic_options: Option<EicOptions>,
    pub find_peaks: Option<FindPeaksOptions>,
}

impl Default for FindFeatureOptions {
    fn default() -> Self {
        Self {
            scan_eic_options: Some(EicOptions {
                ppm_tolerance: 10.0,
                mz_tolerance: 0.003,
                ..Default::default()
            }),
            eic_options: Some(EicOptions {
                ppm_tolerance: 20.0,
                mz_tolerance: 0.005,
                ..Default::default()
            }),
            find_peaks: Some(FindPeaksOptions::default()),
        }
    }
}

pub fn find_feature(
    mzml: &MzML,
    rois: &[&EicRoi],
    cores: usize,
    options: Option<FindFeatureOptions>,
) -> Vec<Option<Feature>> {
    if rois.is_empty() {
        return Vec::new();
    }

    #[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
    {
        return rois
            .iter()
            .map(|roi| find_one_feature(mzml, roi, options.clone()))
            .collect();
    }

    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    {
        if cores <= 1 {
            return rois
                .iter()
                .map(|roi| find_one_feature(mzml, roi, options.clone()))
                .collect();
        }
        let pool = ThreadPoolBuilder::new()
            .num_threads(cores)
            .build()
            .expect("thread pool");
        pool.install(|| {
            rois.par_iter()
                .map(|roi| find_one_feature(mzml, roi, options.clone()))
                .collect::<Vec<_>>()
        })
    }
}

pub fn find_one_feature(
    mzml: &MzML,
    roi: &EicRoi,
    options: Option<FindFeatureOptions>,
) -> Option<Feature> {
    let opts = options.unwrap_or_default();
    let scan_opts = opts.scan_eic_options.unwrap_or(EicOptions {
        ppm_tolerance: 10.0,
        mz_tolerance: 0.003,
        ..Default::default()
    });
    let eic_opts = opts.eic_options.unwrap_or(EicOptions {
        ppm_tolerance: 20.0,
        mz_tolerance: 0.005,
        ..Default::default()
    });
    let fp_opts = opts.find_peaks.unwrap_or_default();

    let time_window = FromTo {
        from: roi.rt - roi.window,
        to: roi.rt + roi.window,
    };
    let (rts, scans): (Vec<f64>, Vec<CentroidScan>) =
        collect_scans(mzml, time_window, eic_opts.time_unit, 1, false);
    if scans.is_empty() || rts.is_empty() {
        return None;
    }

    let refined_mz =
        max_intensity_centroid(&scans, &rts, time_window.from, time_window.to, scan_opts).unwrap();

    let y = compute_eic_for_mz(&scans, rts.len(), refined_mz, eic_opts);
    let data = DataXY { x: rts, y };

    let picked = get_peak(
        &data,
        &Roi {
            rt: roi.rt,
            window: roi.window,
        },
        Some(fp_opts),
    )?;
    Some(Feature {
        id: roi.id.clone(),
        mz: refined_mz,
        rt: picked.rt,
        peak: picked,
    })
}
