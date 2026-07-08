use msutils::utilities::calculate_eic::{EicOptions, SpectrumKind};
use msutils::utilities::find_features::{
    Feature, FeatureError, FindFeaturesOptions, Scan, detect_features,
};
use msutils::utilities::find_peaks::{ArtifactFilter, FindPeaksOptions, PeakFilter};
use msutils::utilities::shape_filter::PeakShape;

const TARGET_MZ: f64 = 500.0;
const APEX_RT: f64 = 1.0;
const MIN_INTENSITY: f64 = 1000.0;

fn profile_scan(center: f64, points: usize, spacing: f64, apex: f64) -> Scan {
    let half = (points as f64 - 1.0) / 2.0;
    let mut mz = Vec::with_capacity(points);
    let mut intensity = Vec::with_capacity(points);
    for step in 0..points {
        let offset = step as f64 - half;
        mz.push(center + offset * spacing);
        let shape = (1.0 - (offset / (half + 0.5)).powi(2)).max(0.0);
        intensity.push(apex * shape);
    }
    Scan { mz, intensity }
}

fn profile_run(chrom_apex: f64, profile_points: usize) -> (Vec<f64>, Vec<Scan>) {
    let count = 41;
    let mut time = Vec::with_capacity(count);
    let mut scans = Vec::with_capacity(count);
    for step in 0..count {
        let rt = step as f64 * 0.05;
        let height = chrom_apex * (-0.5 * ((rt - APEX_RT) / 0.12).powi(2)).exp();
        time.push(rt);
        scans.push(profile_scan(TARGET_MZ, profile_points, 0.001, height));
    }
    (time, scans)
}

fn detection_options() -> FindFeaturesOptions {
    let mut options = FindFeaturesOptions::default();
    options.seed_eic_options = EicOptions {
        ppm_tolerance: 20.0,
        mz_tolerance: 0.005,
        ..Default::default()
    };
    options.final_eic_options = EicOptions {
        ppm_tolerance: 20.0,
        mz_tolerance: 0.005,
        ..Default::default()
    };
    options.peak_options = FindPeaksOptions {
        filter: Some(PeakFilter {
            auto_noise: Some(true),
            auto_baseline: Some(true),
            min_peak_width_points: Some(5),
            min_intensity: Some(MIN_INTENSITY),
            min_snr: Some(3.0),
            ..Default::default()
        }),
        artifact_filter: Some(ArtifactFilter {
            min_r2: 0.0,
            shape: PeakShape::EMG,
        }),
        ..Default::default()
    };
    options
}

fn feature_at_apex(features: &[Feature]) -> Option<&Feature> {
    features
        .iter()
        .find(|feature| (feature.rt - APEX_RT).abs() < 0.2)
}

fn detect(chrom_apex: f64, profile_points: usize, kind: SpectrumKind) -> Vec<Feature> {
    let (time, scans) = profile_run(chrom_apex, profile_points);
    let grid = vec![TARGET_MZ];
    detect_features(&scans, &grid, &time, kind, &detection_options(), 1)
        .expect("detection should run")
}

#[test]
fn detection_keeps_profile_peak_that_summed_intensity_clears() {
    let features = detect(700.0, 13, SpectrumKind::Profile);
    assert!(
        feature_at_apex(&features).is_some(),
        "a real profile peak whose summed window clears min_intensity must be detected, \
         even though its tallest single profile point stays below the threshold"
    );
}

#[test]
fn detection_reports_intensity_on_the_summed_scale() {
    let features = detect(700.0, 13, SpectrumKind::Profile);
    let feature = feature_at_apex(&features).expect("peak must be detected");
    assert!(
        feature.intensity > 3000.0,
        "detected intensity {} must be on the summed scale (whole profile peak), \
         not the single tallest point (~700)",
        feature.intensity
    );
}

#[test]
fn detection_keeps_strong_profile_peak() {
    let features = detect(3000.0, 13, SpectrumKind::Profile);
    assert!(
        feature_at_apex(&features).is_some(),
        "a strong profile peak whose tallest point already clears min_intensity stays detected"
    );
}

#[test]
fn detection_keeps_centroid_peak() {
    let features = detect(2000.0, 1, SpectrumKind::Centroid);
    assert!(
        feature_at_apex(&features).is_some(),
        "a centroid peak above min_intensity is detected the same way, since one point sums to itself"
    );
}

fn add_bell(center: f64, apex: f64, points: usize, spacing: f64, mz: &mut Vec<f64>, intensity: &mut Vec<f64>) {
    let half = (points as f64 - 1.0) / 2.0;
    for step in 0..points {
        let offset = step as f64 - half;
        mz.push(center + offset * spacing);
        let shape = (1.0 - (offset / (half + 0.5)).powi(2)).max(0.0);
        intensity.push(apex * shape);
    }
}

fn multi_ion_run(ions: &[(f64, f64)], profile_points: usize) -> (Vec<f64>, Vec<Scan>) {
    let count = 41;
    let mut time = Vec::with_capacity(count);
    let mut scans = Vec::with_capacity(count);
    for step in 0..count {
        let rt = step as f64 * 0.05;
        let factor = (-0.5 * ((rt - APEX_RT) / 0.12).powi(2)).exp();
        let mut mz = Vec::new();
        let mut intensity = Vec::new();
        for &(center, apex) in ions {
            add_bell(center, apex * factor, profile_points, 0.001, &mut mz, &mut intensity);
        }
        time.push(rt);
        scans.push(Scan { mz, intensity });
    }
    (time, scans)
}

fn grid_over(low: f64, high: f64, step: f64) -> Vec<f64> {
    let mut grid = Vec::new();
    let mut value = low;
    while value <= high + 1e-9 {
        grid.push(value);
        value += step;
    }
    grid
}

fn detect_over_grid(ions: &[(f64, f64)], profile_points: usize, grid: &[f64]) -> Vec<Feature> {
    let (time, scans) = multi_ion_run(ions, profile_points);
    detect_features(&scans, grid, &time, SpectrumKind::Profile, &detection_options(), 1)
        .expect("detection should run")
}

#[test]
fn detection_keeps_one_profile_bell_as_one_mass() {
    let grid = grid_over(499.99, 500.01, 0.005);
    let features = detect_over_grid(&[(TARGET_MZ, 3000.0)], 13, &grid);
    assert!(feature_at_apex(&features).is_some(), "the single profile bell is detected");
    let widest = features.iter().map(|feature| feature.mz).fold(f64::MIN, f64::max);
    let narrowest = features.iter().map(|feature| feature.mz).fold(f64::MAX, f64::min);
    assert!(
        widest - narrowest < 0.03,
        "one ion must stay one mass; scanning it across {} grid points spread the reported m/z by {}",
        grid.len(),
        widest - narrowest
    );
    assert!(
        features.len() <= 2,
        "one ion must not shatter into many features; got {}",
        features.len()
    );
}

#[test]
fn detection_resolves_two_separated_ions() {
    let grid = grid_over(499.98, 500.08, 0.005);
    let features = detect_over_grid(&[(500.0, 3000.0), (500.06, 3000.0)], 13, &grid);
    let near = |mz: f64| {
        features
            .iter()
            .any(|feature| (feature.mz - mz).abs() < 0.01 && (feature.rt - APEX_RT).abs() < 0.2)
    };
    assert!(near(500.0), "the first ion is detected as its own mass");
    assert!(near(500.06), "the second ion is detected as its own mass");
}

#[test]
fn detection_on_empty_scans_reports_no_candidates() {
    let time: Vec<f64> = (0..41).map(|step| step as f64 * 0.05).collect();
    let scans: Vec<Scan> = (0..41)
        .map(|_| Scan {
            mz: Vec::new(),
            intensity: Vec::new(),
        })
        .collect();
    let grid = vec![TARGET_MZ];
    let result = detect_features(&scans, &grid, &time, SpectrumKind::Profile, &detection_options(), 1);
    assert!(
        matches!(result, Err(FeatureError::NoCandidateMasses)),
        "empty scans must report no candidates, not panic"
    );
}
