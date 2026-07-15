use msutils::utilities::{
    fit_peak::{PeakSeed, PeakShape, draw_peak, fit_peak},
    functions::{emg_fn, gaussian_fn, sigma_from_fwhm},
    structs::DataXY,
};

mod helpers;
use helpers::load_chromatograms;

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
    let y = x
        .iter()
        .map(|&v| gaussian_fn(v, intensity, rt, sigma))
        .collect();
    let data = DataXY { x, y };

    let params = fit_peak(&data, &PeakSeed { rt, intensity }, PeakShape::Gaussian).expect("fit");
    assert_eq!(params.shape, PeakShape::Gaussian);
    assert!(params.r2 > 0.99, "r2 = {:.4}", params.r2);
    assert!(
        (params.fwhm - fwhm).abs() < 0.02,
        "fwhm = {:.4}",
        params.fwhm
    );

    let drawn = draw_peak(&data, &params);
    assert_eq!(drawn.x, data.x);
    let peak = drawn.y.iter().cloned().fold(f64::MIN, f64::max);
    assert!(
        (peak - intensity).abs() < intensity * 0.02,
        "peak = {:.1}",
        peak
    );
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

#[test]
fn real_peaks_fit_emg() {
    let chromatograms = load_chromatograms("chromatogram.ion");
    let peaks: Vec<(String, helpers::Eic)> = chromatograms
        .into_iter()
        .filter(|(id, _)| id.starts_with("emg_filter/"))
        .collect();
    assert!(!peaks.is_empty(), "no emg_filter chromatograms");

    for (id, eic) in peaks {
        let data = DataXY {
            x: eic.time.clone(),
            y: eic.intensity.clone(),
        };
        let seed = apex(&data);
        let params = fit_peak(&data, &seed, PeakShape::EMG).expect(&id);
        assert!(params.r2 >= 0.9, "{} r2 = {:.3}", id, params.r2);

        let drawn = draw_peak(&data, &params);
        assert_eq!(drawn.x, data.x);
    }
}
