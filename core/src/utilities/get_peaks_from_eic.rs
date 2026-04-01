use ionic::SpectrumSource;

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use rayon::{ThreadPoolBuilder, prelude::*};

use crate::utilities::calculate_eic::{
    CentroidScan, EicOptions, TimeUnit, collect_scans, compute_eic_for_mz, lower_bound, upper_bound,
};
use crate::utilities::find_noise_level;
use crate::utilities::find_peaks::{FilterPeaksOptions, FindPeaksOptions};
use crate::utilities::get_peak::get_peak;
use crate::utilities::structs::{DataXY, EicRoi, FromTo, Peak, Roi};

pub fn get_peaks_from_eic<'a>(
    source: &mut impl SpectrumSource,
    from_to: FromTo,
    rois: &'a [EicRoi],
    options: Option<FindPeaksOptions>,
    cores: usize,
) -> Option<Vec<(&'a str, f64, f64, Peak)>> {
    let (rts_full, scans_full) = collect_scans(source, from_to, TimeUnit::Minutes, 1);

    if rts_full.len() < 3 || scans_full.is_empty() {
        return Some(
            rois.iter()
                .map(|r| (r.id.as_str(), r.rt, r.mz, Peak::default()))
                .collect(),
        );
    }

    let f = |roi: &'a EicRoi| -> (&'a str, f64, f64, Peak) {
        compute_one(roi, &options, &rts_full, &scans_full, from_to)
    };

    #[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
    {
        return Some(rois.iter().map(f).collect());
    }

    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    {
        if cores <= 1 || rois.len() < 2 {
            return Some(rois.iter().map(f).collect());
        }
        let pool = ThreadPoolBuilder::new().num_threads(cores).build().ok()?;
        Some(pool.install(|| rois.par_iter().map(f).collect()))
    }
}

#[inline]
fn compute_one<'a>(
    roi: &'a EicRoi,
    options: &Option<FindPeaksOptions>,
    rts_full: &[f64],
    scans_full: &[CentroidScan],
    from_to: FromTo,
) -> (&'a str, f64, f64, Peak) {
    const HALF_WIDTH: f64 = 0.5;
    let local_from = (roi.rt - HALF_WIDTH).max(from_to.from);
    let local_to = (roi.rt + HALF_WIDTH).min(from_to.to);
    if !(local_to > local_from) {
        return (&roi.id, roi.rt, roi.mz, Peak::default());
    }

    let eic_opts = EicOptions::default();

    let start = lower_bound(rts_full, local_from);
    let end = upper_bound(rts_full, local_to).min(rts_full.len());
    if end.saturating_sub(start) < 3 {
        return (&roi.id, roi.rt, roi.mz, Peak::default());
    }

    let local_rts = &rts_full[start..end];
    let local_scans = &scans_full[start..end];

    let y_local = compute_eic_for_mz(local_scans, local_rts.len(), roi.mz, eic_opts);
    let y_full = compute_eic_for_mz(scans_full, rts_full.len(), roi.mz, eic_opts);

    if y_local.len() < 3 {
        return (&roi.id, roi.rt, roi.mz, Peak::default());
    }

    let noise = find_noise_level(&y_full);
    let max_y = y_local.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let snr = if noise.intensity > 0.0 {
        max_y / noise.intensity
    } else {
        0.0
    };

    let mut local_options = options.clone().unwrap_or_default();
    let filter = local_options
        .filter_peaks_options
        .get_or_insert_with(Default::default);
    let mut width_threshold = filter.width_threshold.unwrap_or_default();
    let mut intensity_threshold = filter.intensity_threshold.unwrap_or_default();
    if snr <= 5.0 {
        width_threshold /= 2;
        intensity_threshold /= 2.0;
    }
    local_options.filter_peaks_options = Some(FilterPeaksOptions {
        width_threshold: Some(width_threshold),
        intensity_threshold: Some(intensity_threshold),
        ..local_options.filter_peaks_options.unwrap_or_default()
    });

    let roi_hint = Roi {
        rt: roi.rt,
        window: roi.window,
    };

    let pk = match get_peak(
        &DataXY {
            x: rts_full.to_vec(),
            y: y_full,
        },
        &roi_hint,
        Some(local_options.clone()),
    ) {
        Some(p) => p,
        None => match get_peak(
            &DataXY {
                x: local_rts.to_vec(),
                y: y_local,
            },
            &roi_hint,
            Some(local_options),
        ) {
            Some(p) => p,
            None => Peak::default(),
        },
    };

    (&roi.id, roi.rt, roi.mz, pk)
}
