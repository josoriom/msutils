use msutils::utilities::fit_peak::{PeakSeed, PeakShape, draw_peak, fit_peak};
use msutils::utilities::functions::{emg_fn, gaussian_fn, sigma_from_fwhm};
use msutils::utilities::structs::DataXY;

use serde::Deserialize;
use std::{fs, path::Path};

fn linspace(from: f64, to: f64, points: usize) -> Vec<f64> {
    (0..points)
        .map(|i| from + (to - from) * i as f64 / (points - 1) as f64)
        .collect()
}

fn apex(data: &DataXY) -> PeakSeed {
    let mut rt = data.x[0];
    let mut intensity = data.y[0];
    for (x, y) in data.x.iter().zip(data.y.iter()) {
        if *y > intensity {
            intensity = *y;
            rt = *x;
        }
    }
    PeakSeed { rt, intensity }
}

#[test]
fn gaussian_round_trip() {
    let rt = 5.0;
    let intensity = 1000.0;
    let fwhm = 0.2;
    let sigma = sigma_from_fwhm(fwhm);
    let x = linspace(4.0, 6.0, 200);
    let y = x.iter().map(|&v| gaussian_fn(v, intensity, rt, sigma)).collect();
    let data = DataXY { x, y };

    let params = fit_peak(&data, &PeakSeed { rt, intensity }, PeakShape::Gaussian).expect("fit");
    assert_eq!(params.shape, PeakShape::Gaussian);
    assert!(params.r2 > 0.99, "r2 = {:.4}", params.r2);
    assert!((params.fwhm - fwhm).abs() < 0.02, "fwhm = {:.4}", params.fwhm);

    let drawn = draw_peak(&data, &params);
    assert_eq!(drawn.x, data.x);
    let peak = drawn.y.iter().cloned().fold(f64::MIN, f64::max);
    assert!((peak - intensity).abs() < intensity * 0.02, "peak = {:.1}", peak);
}

#[test]
fn emg_round_trip() {
    let x = linspace(4.0, 7.0, 300);
    let y = x
        .iter()
        .map(|&v| emg_fn(v, 1000.0, 4.9, sigma_from_fwhm(0.15), 0.1))
        .collect();
    let data = DataXY { x, y };
    let seed = apex(&data);

    let params = fit_peak(&data, &seed, PeakShape::EMG).expect("fit");
    assert_eq!(params.shape, PeakShape::EMG);
    assert!(params.r2 > 0.98, "r2 = {:.4}", params.r2);

    let drawn = draw_peak(&data, &params);
    assert_eq!(drawn.x.len(), data.x.len());
    assert_eq!(drawn.x, data.x);
}

#[derive(Deserialize)]
struct Peak {
    name: String,
    x: Vec<f64>,
    y: Vec<f64>,
}

#[test]
fn real_peaks_fit_emg() {
    let path = Path::new(file!()).parent().unwrap().join("fixtures").join("emg_filter.json");
    let text = fs::read_to_string(&path).unwrap_or_else(|e| panic!("cannot read {:?}: {}", path, e));
    let peaks: Vec<Peak> = serde_json::from_str(&text).expect("invalid JSON");

    for peak in peaks {
        let data = DataXY {
            x: peak.x.clone(),
            y: peak.y.clone(),
        };
        let seed = apex(&data);
        let params = fit_peak(&data, &seed, PeakShape::EMG).expect(&peak.name);
        assert!(params.r2 >= 0.9, "{} r2 = {:.3}", peak.name, params.r2);

        let drawn = draw_peak(&data, &params);
        assert_eq!(drawn.x, data.x);
    }
}
