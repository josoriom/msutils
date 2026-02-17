use octo::MzML;
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;

use crate::utilities::calculate_eic::TimeUnit;
use crate::utilities::calculate_eic::{EicOptions, collect_ms1_scans, compute_eic_for_mz};
use crate::utilities::find_peaks::FindPeaksOptions;
use crate::utilities::get_peak::get_peak;
use crate::utilities::structs::{DataXY, EicRoi, FromTo, Peak, Roi};

pub fn get_peaks_from_eic(
    mzml: &MzML,
    from_to: FromTo,
    rois: &[EicRoi],
    options: Option<FindPeaksOptions>,
    cores: usize,
) -> Option<Vec<(String, f64, f64, Peak)>> {
    if cores <= 1 || rois.len() < 2 {
        let mut out: Vec<(String, f64, f64, Peak)> = Vec::with_capacity(rois.len());
        for roi in rois {
            out.push(compute_one(&mzml, from_to, roi, &options));
        }
        return Some(out);
    }
    let pool = ThreadPoolBuilder::new().num_threads(cores).build().ok()?;
    Some(pool.install(|| {
        rois.par_iter()
            .map(|roi| compute_one(&mzml, from_to, roi, &options))
            .collect()
    }))
}

#[inline]
fn compute_one(
    mzml: &MzML,
    from_to: FromTo,
    roi: &EicRoi,
    options: &Option<FindPeaksOptions>,
) -> (String, f64, f64, Peak) {
    const HALF_WIDTH: f64 = 0.5;
    let local_from = (roi.rt - HALF_WIDTH).max(from_to.from);
    let local_to = (roi.rt + HALF_WIDTH).min(from_to.to);

    if !(local_to > local_from) {
        return (roi.id.clone(), roi.rt, roi.mz, Peak::default());
    }

    let local_from_to = FromTo {
        from: local_from,
        to: local_to,
    };

    let (rts, scans) = collect_ms1_scans(mzml, local_from_to, TimeUnit::Minutes);
    if rts.len() < 3 || scans.is_empty() {
        return (roi.id.clone(), roi.rt, roi.mz, Peak::default());
    }

    let eic_opts = EicOptions {
        ppm_tolerance: 20.0,
        mz_tolerance: 0.005,
        ..Default::default()
    };
    let y = compute_eic_for_mz(&scans, rts.len(), roi.mz, eic_opts);
    if y.len() < 3 || y.len() != rts.len() {
        return (roi.id.clone(), roi.rt, roi.mz, Peak::default());
    }

    let pk = match get_peak(
        &DataXY { x: rts, y },
        Roi {
            rt: roi.rt,
            window: roi.window,
        },
        options.clone(),
    ) {
        Some(p) => p,
        None => Peak::default(),
    };
    (roi.id.clone(), roi.rt, roi.mz, pk)
}
