use msutils::utilities::emg_filter::EmgFilter;
use msutils::utilities::functions::{emg_fn, gaussian_fn, sigma_from_fwhm};

use serde::Deserialize;
use std::{fs, path::Path};

fn linspace(from: f64, to: f64, points: usize) -> Vec<f64> {
    (0..points)
        .map(|i| from + (to - from) * i as f64 / (points - 1) as f64)
        .collect()
}

#[derive(Deserialize)]
struct Peak {
    name: String,
    x: Vec<f64>,
    y: Vec<f64>,
}

fn load_peaks_from(file: &str) -> Vec<Peak> {
    let path = Path::new(file!()).parent().unwrap().join("fixtures").join(file);
    let text = fs::read_to_string(&path).unwrap_or_else(|e| panic!("cannot read {:?}: {}", path, e));
    serde_json::from_str(&text).unwrap_or_else(|e| panic!("invalid JSON in {}: {}", file, e))
}

fn load_peaks() -> Vec<Peak> {
    load_peaks_from("emg_filter.json")
}

#[test]
fn real_peaks_fit_emg_with_high_r2() {
    let filter = EmgFilter::new(0.7);
    for peak in load_peaks() {
        let r2 = filter.score(&peak.x, &peak.y).expect("peak should be fittable");
        assert!(r2 >= 0.8, "{} fit r2 = {:.3}", peak.name, r2);
        assert!(filter.keeps_peak(&peak.x, &peak.y), "{} should pass at 0.7", peak.name);
    }
}

#[test]
fn tailing_void_peaks_fit_emg_with_high_r2() {
    let filter = EmgFilter::new(0.7);
    for peak in load_peaks_from("emg_tailing.json") {
        let r2 = filter.score(&peak.x, &peak.y).expect("peak should be fittable");
        assert!(r2 >= 0.75, "{} fit r2 = {:.3}", peak.name, r2);
        assert!(filter.keeps_peak(&peak.x, &peak.y), "{} should pass at 0.7", peak.name);
    }
}

#[test]
fn flat_noise_is_rejected() {
    let y = vec![
        0.40, 0.55, 0.45, 0.60, 0.50, 0.58, 0.48, 0.61, 0.52, 0.57, 0.49, 0.60, 0.50, 0.56, 0.47,
        0.59, 0.51, 0.55, 0.46, 0.58, 0.50, 0.54,
    ];
    let x: Vec<f64> = (0..y.len()).map(|i| i as f64 * 0.01).collect();
    let filter = EmgFilter::new(0.7);
    let r2 = filter.score(&x, &y).expect("scored");
    assert!(r2 < 0.7, "flat noise r2 = {:.3} should be < 0.7", r2);
    assert!(!filter.keeps_peak(&x, &y), "flat noise should be rejected at 0.7");
}

#[test]
fn empty_or_short_input_is_unfittable() {
    let filter = EmgFilter::new(0.7);
    assert_eq!(filter.score(&[], &[]), None);
    assert_eq!(filter.score(&[0.0, 1.0, 2.0], &[1.0, 2.0, 1.0]), None);
    assert!(filter.keeps_peak(&[], &[]));
}

#[test]
fn all_zero_region_is_unfittable_and_kept() {
    let filter = EmgFilter::new(0.7);
    let x = linspace(0.0, 1.0, 20);
    let y = vec![0.0; 20];
    assert_eq!(filter.score(&x, &y), None);
    assert!(filter.keeps_peak(&x, &y));
}

#[test]
fn non_finite_values_do_not_panic() {
    let filter = EmgFilter::new(0.7);
    let x = linspace(4.0, 6.0, 40);
    let mut y: Vec<f64> = x.iter().map(|&v| gaussian_fn(v, 1000.0, 5.0, sigma_from_fwhm(0.2))).collect();
    y[8] = f64::NAN;
    y[20] = f64::INFINITY;
    let _ = filter.score(&x, &y);
    let _ = filter.keeps_peak(&x, &y);
}

#[test]
fn clean_gaussian_passes() {
    let filter = EmgFilter::new(0.7);
    let x = linspace(4.0, 6.0, 120);
    let y: Vec<f64> = x.iter().map(|&v| gaussian_fn(v, 1000.0, 5.0, sigma_from_fwhm(0.2))).collect();
    let r2 = filter.score(&x, &y).expect("fittable");
    assert!(r2 >= 0.95, "clean gaussian r2 = {:.3}", r2);
    assert!(filter.keeps_peak(&x, &y));
}

#[test]
fn clean_emg_passes() {
    let filter = EmgFilter::new(0.7);
    let x = linspace(4.0, 7.0, 150);
    let y: Vec<f64> = x.iter().map(|&v| emg_fn(v, 1000.0, 4.9, sigma_from_fwhm(0.15), 0.1)).collect();
    let r2 = filter.score(&x, &y).expect("fittable");
    assert!(r2 >= 0.95, "clean emg r2 = {:.3}", r2);
    assert!(filter.keeps_peak(&x, &y));
}

#[test]
fn off_for_zero_and_negative_threshold() {
    assert!(EmgFilter::new(0.0).is_off());
    assert!(EmgFilter::new(-0.5).is_off());
    assert!(!EmgFilter::new(0.5).is_off());
}

#[test]
fn zero_threshold_turns_filter_off() {
    let filter = EmgFilter::new(0.0);
    assert!(filter.is_off());
    assert!(filter.keeps_peak(&[0.0, 1.0, 2.0, 3.0], &[5.0, 5.0, 5.0, 5.0]));
    assert!(filter.keeps_peak(&[], &[]));
}
