use rayon::ThreadPoolBuilder;
use rayon::prelude::*;

use crate::utilities::calculate_eic::{
    CentroidScan, EicOptions, collect_ms1_scans, compute_eic_for_mz,
};
use crate::utilities::find_features::{refine_mz_for_peak, weighted_centroid};
use crate::utilities::find_peaks::FindPeaksOptions;
use crate::utilities::get_peak::get_peak;
use crate::utilities::structs::{DataXY, EicRoi, FromTo, Peak, Roi};
use octo::MzML;

#[derive(Clone, Debug)]
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

    if cores <= 1 {
        let mut out = Vec::with_capacity(rois.len());
        for roi in rois {
            out.push(find_one_feature(mzml, roi, options.clone()));
        }
        return out;
    }

    let pool = ThreadPoolBuilder::new()
        .num_threads(cores)
        .build()
        .expect("thread pool");
    let opts = options.clone();
    pool.install(|| {
        rois.par_iter()
            .map(|roi| find_one_feature(mzml, roi, opts.clone()))
            .collect::<Vec<_>>()
    })
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
        collect_ms1_scans(mzml, time_window, eic_opts.time_unit);
    if scans.is_empty() || rts.is_empty() {
        return None;
    }

    let mut bins_scratch: Vec<f64> = Vec::new();
    let refined_approx = refine_mz_for_peak(
        &scans,
        &rts,
        roi.mz,
        time_window.from,
        time_window.to,
        scan_opts,
        &mut bins_scratch,
    )
    .unwrap_or(roi.mz);

    let refined_mz = weighted_centroid(
        &scans,
        &rts,
        time_window.from,
        time_window.to,
        refined_approx,
        scan_opts,
    )
    .unwrap_or(refined_approx);

    let y = compute_eic_for_mz(&scans, rts.len(), refined_mz, eic_opts);
    let data = DataXY { x: rts, y };

    let picked = get_peak(
        &data,
        Roi {
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
