use octo::MzML;
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;

use crate::utilities::calculate_eic::TimeUnit;
use crate::utilities::calculate_eic::{EicOptions, collect_ms1_scans, compute_eic_for_mz};
use crate::utilities::find_noise_level;
use crate::utilities::find_peaks::FilterPeaksOptions;
use crate::utilities::find_peaks::FindPeaksOptions;
use crate::utilities::get_peak::get_peak;
use crate::utilities::structs::{DataXY, EicRoi, FromTo, Peak, Roi};

pub fn get_peaks_from_eic<'a>(
    mzml: &MzML,
    from_to: FromTo,
    rois: &'a [EicRoi],
    options: Option<FindPeaksOptions>,
    cores: usize,
) -> Option<Vec<(&'a str, f64, f64, Peak)>> {
    if cores <= 1 || rois.len() < 2 {
        let mut out: Vec<(&'a str, f64, f64, Peak)> = Vec::with_capacity(rois.len());
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
fn compute_one<'a>(
    mzml: &MzML,
    from_to: FromTo,
    roi: &'a EicRoi,
    options: &Option<FindPeaksOptions>,
) -> (&'a str, f64, f64, Peak) {
    const HALF_WIDTH: f64 = 0.5;
    let from = (roi.rt - HALF_WIDTH).max(from_to.from);
    let to = (roi.rt + HALF_WIDTH).min(from_to.to);
    if !(to > from) {
        return (&roi.id, roi.rt, roi.mz, Peak::default());
    }
    let eic_opts = EicOptions::default();

    let (rts, scans) = collect_ms1_scans(mzml, FromTo { from, to }, TimeUnit::Minutes);
    let (rts_full, scans_full) = collect_ms1_scans(mzml, from_to, TimeUnit::Minutes);
    if rts.len() < 3 || scans.is_empty() {
        return (&roi.id, roi.rt, roi.mz, Peak::default());
    }

    let y = compute_eic_for_mz(&scans, rts.len(), roi.mz, eic_opts);
    let y_full = compute_eic_for_mz(&scans_full, rts_full.len(), roi.mz, eic_opts);
    if y.len() < 3 || y.len() != rts.len() {
        return (&roi.id, roi.rt, roi.mz, Peak::default());
    }

    let noise = find_noise_level(&y_full);
    let max_y = y.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let snr = if noise.intensity > 0.0 {
        max_y / noise.intensity
    } else {
        0.0
    };

    let mut local_options = options.clone().unwrap_or_default();
    let filter_settings = local_options
        .filter_peaks_options
        .get_or_insert_with(Default::default);

    let mut width_threshold = filter_settings.width_threshold.unwrap_or_default();
    let mut intensity_threshold = filter_settings.intensity_threshold.unwrap_or_default();

    if snr <= 5.0 {
        width_threshold /= 2;
        intensity_threshold /= 2.0;
    }

    local_options.filter_peaks_options = Some(FilterPeaksOptions {
        width_threshold: Some(width_threshold),
        intensity_threshold: Some(intensity_threshold),
        ..local_options.filter_peaks_options.unwrap_or_default()
    });

    let pk = match get_peak(
        &DataXY { x: rts, y },
        Roi {
            rt: roi.rt,
            window: roi.window,
        },
        Some(local_options),
    ) {
        Some(p) => p,
        None => Peak::default(),
    };
    (&roi.id, roi.rt, roi.mz, pk)
}
