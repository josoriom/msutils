use msutils::utilities::find_peaks::{FindPeaksOptions, PeakFilter};
use msutils::utilities::get_peak::get_peak;
use msutils::utilities::structs::{DataXY, Roi};
use serde::Deserialize;
use std::path::PathBuf;

const HALF_WIDTH: f64 = 0.5;
const TOLERANCE: f64 = 0.01;
const MIN_MATCH_RATE: f64 = 0.99;

#[derive(Deserialize)]
struct Feature {
    id: String,
    rt: f64,
    eics: Vec<SampleEic>,
}

#[derive(Deserialize)]
struct SampleEic {
    sample: String,
    x: Vec<f64>,
    y: Vec<f64>,
    integral: Option<f64>,
}

fn options() -> FindPeaksOptions {
    FindPeaksOptions {
        filter: Some(PeakFilter {
            auto_noise: Some(true),
            auto_baseline: Some(true),
            min_peak_width_points: Some(5),
            min_intensity: Some(1000.0),
            min_snr: Some(3.0),
            ..Default::default()
        }),
        ..Default::default()
    }
}

fn path_from(env_key: &str, default: &str) -> PathBuf {
    std::env::var(env_key).map(PathBuf::from).unwrap_or_else(|_| PathBuf::from(default))
}

fn is_match(measured: f64, expected: f64) -> bool {
    let scale = expected.abs().max(1.0);
    (measured - expected).abs() / scale < TOLERANCE
}

fn check_stored_integrals(path: PathBuf, label: &str) {
    let text = std::fs::read_to_string(&path).expect("read eics json");
    let features: Vec<Feature> = serde_json::from_str(&text).expect("parse eics json");

    let mut matched = 0usize;
    let mut total = 0usize;
    let mut misses: Vec<String> = Vec::new();

    for feature in &features {
        for eic in &feature.eics {
            let expected = eic.integral.unwrap_or(0.0);
            total += 1;
            let data = DataXY {
                x: eic.x.clone(),
                y: eic.y.clone(),
            };
            let measured = get_peak(&data, &Roi::new(feature.rt, HALF_WIDTH), Some(options()))
                .map(|peak| peak.integral)
                .unwrap_or(0.0);
            if is_match(measured, expected) {
                matched += 1;
            } else if misses.len() < 20 {
                misses.push(format!(
                    "{} {}: get_peak {:.6e} vs stored {:.6e}",
                    feature.id, eic.sample, measured, expected
                ));
            }
        }
    }

    println!(
        "\n{label} integral parity: {matched}/{total} within {:.0}% of the stored integral",
        TOLERANCE * 100.0
    );
    for line in &misses {
        println!("  {line}");
    }

    assert!(
        matched as f64 >= total as f64 * MIN_MATCH_RATE,
        "{label}: get_peak reproduced only {matched}/{total} stored integrals within {:.0}%",
        TOLERANCE * 100.0
    );
}

#[test]
#[ignore]
fn qehf_get_peak_reproduces_stored_integrals() {
    check_stored_integrals(
        path_from(
            "MTBLS733_QEHF_JSON",
            "/Users/josoriom/Documents/pqc/mtbls733/quant_fail_eics.json",
        ),
        "qehf",
    );
}

#[test]
#[ignore]
fn ttof_get_peak_reproduces_stored_integrals() {
    check_stored_integrals(
        path_from(
            "MTBLS733_TTOF_JSON",
            "/Users/josoriom/Documents/pqc/mtbls733/quant_fail_ttof.json",
        ),
        "ttof",
    );
}
