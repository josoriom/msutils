use ionic::ion::{ByteRange, IonReader};

use crate::utilities::{
    calculate_eic::{
        EicOptions, EicReader, FastError, MS1_LEVEL, get_scan_times, lower_bound, mz_tolerance_for,
        plan_window_ranges, read_mz_window, upper_bound,
    },
    find_peaks::FindPeaksOptions,
    get_peak::get_peak,
    structs::{DataXY, FromTo, Peak, Roi},
};

pub fn plan_peaks_ranges(
    ion: &mut IonReader,
    rois: &[Roi],
    from: f64,
    to: f64,
) -> Result<Vec<ByteRange>, FastError> {
    if rois.is_empty() {
        return Ok(Vec::new());
    }

    let options = EicOptions::default();
    let mut mz_from = f64::INFINITY;
    let mut mz_to = f64::NEG_INFINITY;
    for roi in rois {
        let Roi::Eic { mz, .. } = roi else {
            continue;
        };
        let tolerance = mz_tolerance_for(*mz, options);
        mz_from = mz_from.min(mz - tolerance);
        mz_to = mz_to.max(mz + tolerance);
    }

    plan_window_ranges(ion, from, to, mz_from, mz_to)
}

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use rayon::prelude::*;

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use crate::utilities::parallel::run_with_cores;

struct PeakJob {
    roi_idx: usize,
    rt: f64,
    range: f64,
    mz_from: f64,
    mz_to: f64,
}

pub fn get_peaks_from_eic<'a>(
    reader: &mut EicReader,
    from_to: FromTo,
    rois: &'a [Roi],
    options: Option<FindPeaksOptions>,
    cores: usize,
) -> Result<Vec<(&'a str, f64, f64, Peak)>, FastError> {
    if rois.is_empty() {
        return Ok(Vec::new());
    }

    let scan_times = get_scan_times(reader, from_to.from, from_to.to, MS1_LEVEL);

    if scan_times.len() < 3 {
        return Ok(rois
            .iter()
            .map(|roi| make_row(roi, Peak::default()))
            .collect());
    }

    let rts_full: Vec<f64> = scan_times.iter().map(|s| s.rt).collect();
    let scan_count = scan_times.len();

    let jobs = build_peak_jobs(rois);
    if jobs.is_empty() {
        return Ok(rois
            .iter()
            .map(|roi| make_row(roi, Peak::default()))
            .collect());
    }

    let eic_values = read_eics_one_pass(reader, &scan_times, &jobs)?;

    let peaks = compute_peaks_for_jobs(
        &jobs,
        &eic_values,
        scan_count,
        &rts_full,
        &options,
        from_to,
        cores,
    );

    let mut results = vec![Peak::default(); rois.len()];
    for (roi_idx, peak) in peaks {
        results[roi_idx] = peak;
    }

    Ok(rois
        .iter()
        .enumerate()
        .map(|(idx, roi)| make_row(roi, results[idx]))
        .collect())
}

fn make_row(roi: &Roi, peak: Peak) -> (&str, f64, f64, Peak) {
    match roi {
        Roi::Eic { id, rt, mz, .. } => (id.as_str(), *rt, *mz, peak),
        _ => ("", 0.0, 0.0, Peak::default()),
    }
}

fn build_peak_jobs(rois: &[Roi]) -> Vec<PeakJob> {
    let eic_opts = EicOptions::default();
    rois.iter()
        .enumerate()
        .filter_map(|(idx, roi)| {
            let Roi::Eic { rt, mz, range, .. } = roi else {
                return None;
            };
            let tolerance = mz_tolerance_for(*mz, eic_opts);
            Some(PeakJob {
                roi_idx: idx,
                rt: *rt,
                range: *range,
                mz_from: mz - tolerance,
                mz_to: mz + tolerance,
            })
        })
        .collect()
}

fn compute_peak_for_job(
    job: &PeakJob,
    y_full: &[f64],
    rts_full: &[f64],
    options: &Option<FindPeaksOptions>,
    from_to: FromTo,
) -> Peak {
    const HALF_WIDTH: f64 = 0.5;
    let local_from = (job.rt - HALF_WIDTH).max(from_to.from);
    let local_to = (job.rt + HALF_WIDTH).min(from_to.to);
    if local_to <= local_from {
        return Peak::default();
    }

    let start = lower_bound(rts_full, local_from);
    let end = upper_bound(rts_full, local_to).min(rts_full.len());
    if end.saturating_sub(start) < 3 {
        return Peak::default();
    }

    get_peak(
        &DataXY {
            x: rts_full.to_vec(),
            y: y_full.to_vec(),
        },
        &Roi::peak(job.rt, job.range),
        options.clone(),
    )
}

fn read_eics_one_pass(
    reader: &mut EicReader,
    scan_times: &[crate::utilities::calculate_eic::ScanTime],
    jobs: &[PeakJob],
) -> Result<Vec<f64>, FastError> {
    let scan_count = scan_times.len();
    let mut eic_values = vec![0.0; jobs.len() * scan_count];

    let mz_from = jobs
        .iter()
        .map(|job| job.mz_from)
        .fold(f64::INFINITY, f64::min);
    let mz_to = jobs
        .iter()
        .map(|job| job.mz_to)
        .fold(f64::NEG_INFINITY, f64::max);

    let mut mz_values = Vec::new();
    let mut intensity_values = Vec::new();

    for (scan_position, scan_time) in scan_times.iter().enumerate() {
        read_mz_window(
            reader,
            scan_time.index,
            mz_from,
            mz_to,
            &mut mz_values,
            &mut intensity_values,
        )?;

        for (job_index, job) in jobs.iter().enumerate() {
            let start = lower_bound(&mz_values, job.mz_from);
            let end = upper_bound(&mz_values, job.mz_to).min(mz_values.len());
            eic_values[job_index * scan_count + scan_position] =
                intensity_values[start..end].iter().sum();
        }
    }

    Ok(eic_values)
}

fn compute_peaks_for_jobs(
    jobs: &[PeakJob],
    eic_values: &[f64],
    scan_count: usize,
    rts_full: &[f64],
    options: &Option<FindPeaksOptions>,
    from_to: FromTo,
    cores: usize,
) -> Vec<(usize, Peak)> {
    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    if cores > 1 && jobs.len() > 1 {
        return run_with_cores(cores, || {
            compute_peaks_for_jobs_parallel(
                jobs, eic_values, scan_count, rts_full, options, from_to,
            )
        });
    }

    #[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
    let _ = cores;

    compute_peaks_for_jobs_serial(jobs, eic_values, scan_count, rts_full, options, from_to)
}

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
fn compute_peaks_for_jobs_parallel(
    jobs: &[PeakJob],
    eic_values: &[f64],
    scan_count: usize,
    rts_full: &[f64],
    options: &Option<FindPeaksOptions>,
    from_to: FromTo,
) -> Vec<(usize, Peak)> {
    (0..jobs.len())
        .into_par_iter()
        .map(|job_index| {
            compute_peak_for_job_index(
                job_index, jobs, eic_values, scan_count, rts_full, options, from_to,
            )
        })
        .collect()
}

fn compute_peaks_for_jobs_serial(
    jobs: &[PeakJob],
    eic_values: &[f64],
    scan_count: usize,
    rts_full: &[f64],
    options: &Option<FindPeaksOptions>,
    from_to: FromTo,
) -> Vec<(usize, Peak)> {
    (0..jobs.len())
        .map(|job_index| {
            compute_peak_for_job_index(
                job_index, jobs, eic_values, scan_count, rts_full, options, from_to,
            )
        })
        .collect()
}

fn compute_peak_for_job_index(
    job_index: usize,
    jobs: &[PeakJob],
    eic_values: &[f64],
    scan_count: usize,
    rts_full: &[f64],
    options: &Option<FindPeaksOptions>,
    from_to: FromTo,
) -> (usize, Peak) {
    let job = &jobs[job_index];
    let row_start = job_index * scan_count;
    let row_end = row_start + scan_count;
    let peak = compute_peak_for_job(
        job,
        &eic_values[row_start..row_end],
        rts_full,
        options,
        from_to,
    );
    (job.roi_idx, peak)
}
