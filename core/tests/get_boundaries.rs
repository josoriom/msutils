use std::{
    collections::BTreeMap,
    path::{Path, PathBuf},
};

use msutils::utilities::{
    find_peaks::{FindPeaksOptions, PeakFilter, find_peaks},
    get_boundaries::{BoundariesOptions, get_boundaries},
    structs::DataXY,
};

mod helpers;
use helpers::load_chromatograms_from;

struct Feature {
    id: String,
    fc_benchmark: f64,
    samples: Vec<Sample>,
}

struct Sample {
    sample: String,
    peak_x: f64,
    noise: f64,
    expected_from: usize,
    expected_to: usize,
    x: Vec<f64>,
    y: Vec<f64>,
}

const ACCURACY_TOLERANCE: f64 = 0.20;
const ALLOWED_GAP: usize = 4;
const MIN_DETECTION_RATE: f64 = 0.95;

fn ion_path() -> PathBuf {
    if let Ok(path) = std::env::var("TEST_ION") {
        return PathBuf::from(path);
    }
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join("test.ion")
}

fn load_cases() -> Vec<Feature> {
    let cases = load_chromatograms_from(&ion_path());
    let mut features: BTreeMap<String, Feature> = BTreeMap::new();
    for (_, eic) in cases {
        if eic.params.get("kind").map(String::as_str) != Some("boundaries") {
            continue;
        }
        let get = |key: &str| {
            eic.params
                .get(key)
                .unwrap_or_else(|| panic!("get_boundaries.ion case missing param {key}"))
        };
        let feature_id = get("feature_id").clone();
        let fc_benchmark: f64 = get("fc_benchmark").parse().expect("fc_benchmark");
        let sample = Sample {
            sample: get("sample").clone(),
            peak_x: get("peak_x").parse().expect("peak_x"),
            noise: get("noise").parse().expect("noise"),
            expected_from: get("expected_from").parse().expect("expected_from"),
            expected_to: get("expected_to").parse().expect("expected_to"),
            x: eic.time.clone(),
            y: eic.intensity.clone(),
        };
        features
            .entry(feature_id.clone())
            .or_insert_with(|| Feature {
                id: feature_id,
                fc_benchmark,
                samples: Vec::new(),
            })
            .samples
            .push(sample);
    }
    features.into_values().collect()
}

fn area(x: &[f64], y: &[f64], from: usize, to: usize) -> f64 {
    if to <= from || to >= x.len() {
        return 0.0;
    }
    let mut sum = 0.0;
    for i in from..to {
        sum += (x[i + 1] - x[i]) * (y[i] + y[i + 1]) / 2.0;
    }
    sum
}

fn measured_window(sample: &Sample) -> (usize, usize) {
    let data = DataXY {
        x: sample.x.clone(),
        y: sample.y.clone(),
    };
    let edges = get_boundaries(
        &data,
        sample.peak_x,
        Some(BoundariesOptions {
            noise: sample.noise,
            ..Default::default()
        }),
    );
    (edges.from.index.expect("from"), edges.to.index.expect("to"))
}

fn fold_change(feature: &Feature, window: impl Fn(&Sample) -> (usize, usize)) -> f64 {
    let mut group_a = Vec::new();
    let mut group_b = Vec::new();
    for sample in &feature.samples {
        let (from, to) = window(sample);
        let value = area(&sample.x, &sample.y, from, to);
        if value > 0.0 {
            if sample.sample.starts_with("SA") {
                group_a.push(value.log2());
            } else {
                group_b.push(value.log2());
            }
        }
    }
    if group_a.is_empty() || group_b.is_empty() {
        return f64::NAN;
    }
    let mean = |values: &[f64]| values.iter().sum::<f64>() / values.len() as f64;
    (mean(&group_b) - mean(&group_a)).exp2()
}

fn is_accurate(measured: f64, benchmark: f64) -> bool {
    measured.is_finite()
        && benchmark.abs() > 0.0
        && (measured - benchmark).abs() / benchmark.abs() < ACCURACY_TOLERANCE
}

#[test]
fn get_boundaries_recovers_quantification() {
    let features = load_cases();

    let mut recovered = 0usize;
    let mut ceiling = 0usize;
    let mut missed: Vec<String> = Vec::new();

    for feature in &features {
        let mine = fold_change(feature, measured_window);
        let ideal = fold_change(feature, |sample| (sample.expected_from, sample.expected_to));

        let mine_ok = is_accurate(mine, feature.fc_benchmark);
        let ideal_ok = is_accurate(ideal, feature.fc_benchmark);
        if mine_ok {
            recovered += 1;
        }
        if ideal_ok {
            ceiling += 1;
        }
        if ideal_ok && !mine_ok {
            missed.push(format!(
                "{}: benchmark {:.3}, ideal labels {:.3} (ok), get_boundaries {:.3} (off)",
                feature.id, feature.fc_benchmark, ideal, mine
            ));
        }
    }

    println!(
        "\nfold-change accuracy on {} features:\n  get_boundaries : {recovered}\n  ideal labels   : {ceiling}",
        features.len()
    );
    for line in &missed {
        println!("  {line}");
    }

    assert!(
        recovered + ALLOWED_GAP >= ceiling,
        "get_boundaries recovered {recovered} features vs the ideal-label ceiling {ceiling} (allowed gap {ALLOWED_GAP})"
    );
}

fn detection_options() -> FindPeaksOptions {
    FindPeaksOptions {
        filter: Some(PeakFilter {
            auto_noise: Some(true),
            auto_baseline: Some(true),
            min_peak_width_points: Some(5),
            min_intensity: Some(1000.0),
            min_snr: Some(2.0),
            ..Default::default()
        }),
        ..Default::default()
    }
}

#[test]
fn find_peaks_locates_apex_in_window() {
    let features = load_cases();
    let mut located = 0usize;
    let mut total = 0usize;

    for feature in &features {
        for sample in &feature.samples {
            total += 1;
            let data = DataXY {
                x: sample.x.clone(),
                y: sample.y.clone(),
            };
            let low = sample.x[sample.expected_from];
            let high = sample.x[sample.expected_to];
            let inside = find_peaks(&data, Some(detection_options()))
                .iter()
                .any(|peak| peak.rt >= low && peak.rt <= high);
            if inside {
                located += 1;
            }
        }
    }

    println!("\nfind_peaks located a peak inside the window: {located}/{total}");
    assert!(
        located as f64 >= total as f64 * MIN_DETECTION_RATE,
        "find_peaks located a peak inside the window in only {located} of {total} cases"
    );
}
