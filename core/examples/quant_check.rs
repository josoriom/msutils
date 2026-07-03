use std::cmp::Ordering;
use std::env;
use std::fs;
use std::io::Write;

use msutils::utilities::find_peaks::{FindPeaksOptions, PeakFilter, find_peaks};
use msutils::utilities::structs::{DataXY, Peak};
use serde_json::Value;

const SEARCH_HALF_WIDTH: f64 = 0.5;

fn options() -> FindPeaksOptions {
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

fn numbers(value: &Value) -> Vec<f64> {
    value
        .as_array()
        .map(|items| items.iter().filter_map(Value::as_f64).collect())
        .unwrap_or_default()
}

fn pick_peak(peaks: &[Peak], rt: f64) -> Option<&Peak> {
    peaks
        .iter()
        .filter(|peak| (peak.rt - rt).abs() <= SEARCH_HALF_WIDTH)
        .max_by(|a, b| a.intensity.partial_cmp(&b.intensity).unwrap_or(Ordering::Equal))
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let text = fs::read_to_string(&args[1]).expect("read input json");
    let features: Value = serde_json::from_str(&text).expect("parse input json");

    let mut lines = vec![String::from("id\tsample\tfrom_rt\tto_rt")];
    for feature in features.as_array().expect("features array") {
        let id = feature["id"].as_str().unwrap_or_default();
        let rt = feature["rt"].as_f64().unwrap_or(f64::NAN);
        for eic in feature["eics"].as_array().expect("eics array") {
            let sample = eic["sample"].as_str().unwrap_or_default();
            let data = DataXY {
                x: numbers(&eic["x"]),
                y: numbers(&eic["y"]),
            };
            match pick_peak(&find_peaks(&data, Some(options())), rt) {
                Some(peak) => lines.push(format!("{id}\t{sample}\t{}\t{}", peak.from, peak.to)),
                None => lines.push(format!("{id}\t{sample}\tNA\tNA")),
            }
        }
    }

    let mut file = fs::File::create(&args[2]).expect("create output");
    file.write_all(lines.join("\n").as_bytes()).expect("write output");
    eprintln!("wrote {} rows", lines.len() - 1);
}
