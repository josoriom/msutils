use crate::utilities::{
    calculate_eic::{EicOptions, ScanQuery, get_eic_for_mz, get_scans, lower_bound, upper_bound},
    find_features::max_intensity_centroid,
    find_peaks::FindPeaksOptions,
    get_peak::get_peak,
    structs::{DataXY, EicRoi, FromTo, Peak, Roi},
};

use ionic::ScanSource;
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
    source: &mut impl ScanSource,
    rois: &[&EicRoi],
    cores: usize,
    options: Option<FindFeatureOptions>,
) -> Vec<Option<Feature>> {
    if rois.is_empty() {
        return Vec::new();
    }

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

    let mut rt_min = f64::MAX;
    let mut rt_max = f64::MIN;
    for roi in rois {
        if roi.rt.is_finite() && roi.window.is_finite() && roi.window > 0.0 {
            rt_min = rt_min.min(roi.rt - roi.window);
            rt_max = rt_max.max(roi.rt + roi.window);
        }
    }

    if !rt_min.is_finite() || !rt_max.is_finite() || rt_max <= rt_min {
        return vec![None; rois.len()];
    }

    let (rts, scans) = get_scans(
        source,
        ScanQuery::RtRange(FromTo {
            from: rt_min,
            to: rt_max,
        }),
        eic_opts.time_unit,
        1,
    );

    if scans.is_empty() || rts.is_empty() {
        return vec![None; rois.len()];
    }

    let process = |roi: &&EicRoi| -> Option<Feature> {
        if !roi.rt.is_finite() || !roi.mz.is_finite() || roi.window <= 0.0 {
            return None;
        }

        let local_from = roi.rt - roi.window;
        let local_to = roi.rt + roi.window;
        let start = lower_bound(&rts, local_from);
        let end = upper_bound(&rts, local_to).min(rts.len());
        if end <= start {
            return None;
        }

        let local_rts = &rts[start..end];
        let local_scans = &scans[start..end];

        let refined_mz = max_intensity_centroid(local_scans, local_rts, roi.rt, roi.mz, scan_opts)?;
        let y = get_eic_for_mz(local_scans, local_rts.len(), refined_mz, eic_opts);
        let data = DataXY {
            x: local_rts.to_vec(),
            y,
        };

        let picked = get_peak(
            &data,
            &Roi {
                rt: roi.rt,
                window: roi.window,
            },
            Some(fp_opts.clone()),
        )?;

        Some(Feature {
            id: roi.id.clone(),
            mz: refined_mz,
            rt: picked.rt,
            peak: picked,
        })
    };

    #[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
    {
        return rois.iter().map(process).collect();
    }

    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    {
        if cores <= 1 {
            return rois.iter().map(process).collect();
        }
        let pool = ThreadPoolBuilder::new()
            .num_threads(cores)
            .build()
            .expect("thread pool");
        pool.install(|| rois.par_iter().map(process).collect::<Vec<_>>())
    }
}
