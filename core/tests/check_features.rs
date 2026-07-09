mod helpers;

use helpers::{Eic, load_chromatograms_from, param_f64};
use msutils::utilities::find_peaks::{ArtifactFilter, FindPeaksOptions, PeakFilter};
use msutils::utilities::get_peak::get_peak;
use msutils::utilities::shape_filter::PeakShape;
use msutils::utilities::structs::{DataXY, Roi};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

const RT_WINDOW: f64 = 0.3;
const EXPECTED_FEATURES: usize = 80;

fn ion_path() -> PathBuf {
    if let Ok(path) = std::env::var("MISSING_ION") {
        return PathBuf::from(path);
    }
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join("check_features.ion")
}

fn peak_options() -> FindPeaksOptions {
    FindPeaksOptions {
        filter: Some(PeakFilter {
            auto_noise: Some(true),
            auto_baseline: Some(true),
            min_peak_width_points: Some(3),
            min_intensity: Some(0.0),
            min_snr: Some(2.0),
            ..Default::default()
        }),
        artifact_filter: Some(ArtifactFilter {
            min_r2: 0.0,
            shape: PeakShape::EMG,
        }),
        ..Default::default()
    }
}

fn detects_peak(eic: &Eic) -> bool {
    if eic.time.len() < 3 {
        return false;
    }
    let data = DataXY {
        x: eic.time.clone(),
        y: eic.intensity.clone(),
    };
    let rt = param_f64(eic, "rt");
    let span = eic.time[eic.time.len() - 1] - eic.time[0];
    match get_peak(&data, &Roi::new(rt, span), Some(peak_options())) {
        Some(peak) => (peak.rt - rt).abs() <= RT_WINDOW,
        None => false,
    }
}

struct Recovery {
    category: String,
    samples: usize,
    detected: usize,
    baseline: usize,
}

fn recovery_by_feature(cases: &BTreeMap<String, Eic>) -> BTreeMap<String, Recovery> {
    let mut by_feature: BTreeMap<String, Recovery> = BTreeMap::new();
    for eic in cases.values() {
        if eic.params.get("kind").map(String::as_str) != Some("missing_feature") {
            continue;
        }
        let compound = eic.params.get("compound").cloned().unwrap_or_default();
        let entry = by_feature.entry(compound).or_insert_with(|| Recovery {
            category: eic.params.get("category").cloned().unwrap_or_default(),
            samples: 0,
            detected: 0,
            baseline: param_f64(eic, "expected_min_detected") as usize,
        });
        entry.samples += 1;
        if detects_peak(eic) {
            entry.detected += 1;
        }
    }
    by_feature
}

#[test]
fn missing_feature_peaks_stay_detectable() {
    let path = ion_path();
    let cases = load_chromatograms_from(&path);
    assert!(!cases.is_empty(), "no cases in {:?}", path);

    let features = recovery_by_feature(&cases);
    assert_eq!(
        features.len(),
        EXPECTED_FEATURES,
        "expected {EXPECTED_FEATURES} missing features in the fixture"
    );

    let mut regressions = Vec::new();
    for (id, r) in &features {
        if r.detected < r.baseline {
            regressions.push(format!(
                "[{id}] {}/{} {} samples now detect the peak, below baseline {}",
                r.detected, r.samples, r.category, r.baseline
            ));
        }
    }

    assert!(
        regressions.is_empty(),
        "{} missing feature(s) lost recoverable peaks:\n{}",
        regressions.len(),
        regressions.join("\n")
    );
}
