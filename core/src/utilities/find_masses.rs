use crate::utilities::{
    calculate_eic::{SpectrumKind, lower_bound, upper_bound},
    find_peaks::{FindPeaksOptions, PeakFilter, find_peaks},
    mz_estimator::SampleMz,
    structs::DataXY,
};

pub fn find_masses(
    mz: &[f64],
    intensity: &[f64],
    mz_lower: f64,
    mz_upper: f64,
    kind: SpectrumKind,
) -> Vec<SampleMz> {
    if mz.len() != intensity.len() {
        return Vec::new();
    }
    let start = lower_bound(mz, mz_lower);
    let end = upper_bound(mz, mz_upper).min(mz.len());
    if end <= start {
        return Vec::new();
    }
    let mz = &mz[start..end];
    let intensity = &intensity[start..end];
    match kind {
        SpectrumKind::Centroid => centroid_masses(mz, intensity),
        SpectrumKind::Profile => profile_masses(mz, intensity),
    }
}

fn centroid_masses(mz: &[f64], intensity: &[f64]) -> Vec<SampleMz> {
    let mut masses = Vec::new();
    for (m, i) in mz.iter().zip(intensity.iter()) {
        if m.is_finite() && i.is_finite() && *i > 0.0 {
            masses.push(SampleMz {
                mz: *m,
                intensity: *i,
            });
        }
    }
    masses
}

fn profile_masses(mz: &[f64], intensity: &[f64]) -> Vec<SampleMz> {
    let profile = DataXY {
        x: mz.to_vec(),
        y: intensity.to_vec(),
    };
    let peaks = find_peaks(&profile, Some(mz_peak_options()));

    let mut masses = Vec::with_capacity(peaks.len());
    for peak in &peaks {
        if let Some(mass) = apex_in_range(mz, intensity, peak.from, peak.to) {
            masses.push(mass);
        }
    }
    if masses.is_empty()
        && let Some(mass) = apex_in_range(mz, intensity, f64::NEG_INFINITY, f64::INFINITY)
    {
        masses.push(mass);
    }
    masses
}

fn apex_in_range(mz: &[f64], intensity: &[f64], from: f64, to: f64) -> Option<SampleMz> {
    let mut apex: Option<SampleMz> = None;
    for (m, i) in mz.iter().zip(intensity.iter()) {
        if !m.is_finite() || !i.is_finite() || *i <= 0.0 || *m < from || *m > to {
            continue;
        }
        if apex.is_none_or(|top| *i > top.intensity) {
            apex = Some(SampleMz {
                mz: *m,
                intensity: *i,
            });
        }
    }
    apex
}

fn mz_peak_options() -> FindPeaksOptions {
    FindPeaksOptions {
        boundaries: Some(Default::default()),
        filter: Some(PeakFilter {
            min_peak_width_points: Some(3),
            auto_noise: Some(false),
            auto_baseline: Some(false),
            min_snr: Some(1.0),
            ..Default::default()
        }),
        baseline: None,
        artifact_filter: None,
    }
}
