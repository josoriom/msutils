use quantion::utilities::{
    calculate_baseline::BaselineOptions,
    find_peaks::{ArtifactFilter, FindPeaksOptions, PeakFilter, find_peaks},
    functions::{emg_fn, gaussian_fn, sigma_from_fwhm},
    shape_filter::{PeakShape, ShapeFilter},
    structs::DataXY,
};

mod helpers;
use helpers::{Eic, load_chromatograms, param_bool, param_f64};

const RT_TOLERANCE: f64 = 0.35;

fn linspace(from: f64, to: f64, points: usize) -> Vec<f64> {
    (0..points)
        .map(|i| from + (to - from) * i as f64 / (points - 1) as f64)
        .collect()
}

fn emg(min_r2: f64) -> ShapeFilter {
    ShapeFilter::new(min_r2, PeakShape::EMG)
}

fn group(prefix: &str) -> Vec<(String, Eic)> {
    let entries: Vec<(String, Eic)> = load_chromatograms("chromatogram.ion")
        .into_iter()
        .filter(|(id, _)| id.starts_with(prefix))
        .collect();
    assert!(!entries.is_empty(), "no chromatograms under {prefix}");
    entries
}

#[test]
fn real_peaks_fit_with_high_r2() {
    let filter = emg(0.7);
    for (id, eic) in group("emg_filter/") {
        let r2 = filter
            .score(&eic.time, &eic.intensity)
            .expect("peak should be fittable");
        assert!(r2 >= 0.8, "{id} fit r2 = {r2:.3}");
        assert!(
            filter.keeps_peak(&eic.time, &eic.intensity),
            "{id} should pass at 0.7"
        );
    }
}

#[test]
fn tailing_void_peaks_fit_with_high_r2() {
    let filter = emg(0.7);
    for (id, eic) in group("emg_tailing/") {
        let r2 = filter
            .score(&eic.time, &eic.intensity)
            .expect("peak should be fittable");
        assert!(r2 >= 0.7, "{id} fit r2 = {r2:.3}");
        assert!(
            filter.keeps_peak(&eic.time, &eic.intensity),
            "{id} should pass at 0.7"
        );
    }
}

#[test]
fn synthetic_cases_match_expectations() {
    for (id, eic) in group("synthetic/") {
        assert_eq!(
            emg(0.7).keeps_peak(&eic.time, &eic.intensity),
            param_bool(&eic, "expect_pass_07"),
            "{id} at threshold 0.7"
        );
        assert_eq!(
            emg(0.9).keeps_peak(&eic.time, &eic.intensity),
            param_bool(&eic, "expect_pass_09"),
            "{id} at threshold 0.9"
        );
    }
}

#[test]
fn flat_noise_is_rejected() {
    let y = vec![
        0.40, 0.55, 0.45, 0.60, 0.50, 0.58, 0.48, 0.61, 0.52, 0.57, 0.49, 0.60, 0.50, 0.56, 0.47,
        0.59, 0.51, 0.55, 0.46, 0.58, 0.50, 0.54,
    ];
    let x: Vec<f64> = (0..y.len()).map(|i| i as f64 * 0.01).collect();
    let filter = emg(0.7);
    let r2 = filter.score(&x, &y).expect("scored");
    assert!(r2 < 0.7, "flat noise r2 = {:.3} should be < 0.7", r2);
    assert!(
        !filter.keeps_peak(&x, &y),
        "flat noise should be rejected at 0.7"
    );
}

#[test]
fn empty_or_short_input_is_unfittable() {
    let filter = emg(0.7);
    assert_eq!(filter.score(&[], &[]), None);
    assert_eq!(filter.score(&[0.0, 1.0, 2.0], &[1.0, 2.0, 1.0]), None);
    assert!(filter.keeps_peak(&[], &[]));
}

#[test]
fn all_zero_region_is_unfittable_and_kept() {
    let filter = emg(0.7);
    let x = linspace(0.0, 1.0, 20);
    let y = vec![0.0; 20];
    assert_eq!(filter.score(&x, &y), None);
    assert!(filter.keeps_peak(&x, &y));
}

#[test]
fn non_finite_values_do_not_panic() {
    let filter = emg(0.7);
    let x = linspace(4.0, 6.0, 40);
    let mut y: Vec<f64> = x
        .iter()
        .map(|&v| gaussian_fn(v, 1000.0, 5.0, sigma_from_fwhm(0.2)))
        .collect();
    y[8] = f64::NAN;
    y[20] = f64::INFINITY;
    let _ = filter.score(&x, &y);
    let _ = filter.keeps_peak(&x, &y);
}

#[test]
fn clean_gaussian_passes() {
    let filter = emg(0.7);
    let x = linspace(4.0, 6.0, 120);
    let y: Vec<f64> = x
        .iter()
        .map(|&v| gaussian_fn(v, 1000.0, 5.0, sigma_from_fwhm(0.2)))
        .collect();
    let r2 = filter.score(&x, &y).expect("fittable");
    assert!(r2 >= 0.95, "clean gaussian r2 = {:.3}", r2);
    assert!(filter.keeps_peak(&x, &y));
}

#[test]
fn clean_gaussian_passes_with_gaussian_shape() {
    let filter = ShapeFilter::new(0.7, PeakShape::Gaussian);
    let x = linspace(4.0, 6.0, 120);
    let y: Vec<f64> = x
        .iter()
        .map(|&v| gaussian_fn(v, 1000.0, 5.0, sigma_from_fwhm(0.2)))
        .collect();
    let r2 = filter.score(&x, &y).expect("fittable");
    assert!(r2 >= 0.95, "clean gaussian (gaussian shape) r2 = {:.3}", r2);
    assert!(filter.keeps_peak(&x, &y));
}

#[test]
fn clean_emg_passes() {
    let filter = emg(0.7);
    let x = linspace(4.0, 7.0, 150);
    let y: Vec<f64> = x
        .iter()
        .map(|&v| emg_fn(v, 1000.0, 4.9, sigma_from_fwhm(0.15), 0.1))
        .collect();
    let r2 = filter.score(&x, &y).expect("fittable");
    assert!(r2 >= 0.95, "clean emg r2 = {:.3}", r2);
    assert!(filter.keeps_peak(&x, &y));
}

#[test]
fn off_for_zero_and_negative_threshold() {
    assert!(emg(0.0).is_off());
    assert!(emg(-0.5).is_off());
    assert!(!emg(0.5).is_off());
}

#[test]
fn zero_threshold_turns_filter_off() {
    let filter = emg(0.0);
    assert!(filter.is_off());
    assert!(filter.keeps_peak(&[0.0, 1.0, 2.0, 3.0], &[5.0, 5.0, 5.0, 5.0]));
    assert!(filter.keeps_peak(&[], &[]));
}

fn find_peaks_options(min_r2: f64) -> FindPeaksOptions {
    FindPeaksOptions {
        boundaries: Some(Default::default()),
        filter: Some(PeakFilter {
            min_integral: None,
            min_intensity: Some(1000.0),
            min_peak_width_points: Some(5),
            noise: None,
            auto_noise: Some(true),
            auto_baseline: Some(true),
            allow_overlap: Some(false),
            min_snr: Some(2.0),
            noise_method: None,
            kernel_size: None,
        }),
        baseline: Some(BaselineOptions {
            lambda: None,
            max_iterations: None,
            edge_slope_level: Some(1),
        }),
        artifact_filter: Some(ArtifactFilter {
            min_r2,
            shape: PeakShape::EMG,
        }),
    }
}

fn has_peak_near(eic: &Eic, rt: f64, min_r2: f64) -> bool {
    let data = DataXY {
        x: eic.time.clone(),
        y: eic.intensity.clone(),
    };
    find_peaks(&data, Some(find_peaks_options(min_r2)))
        .iter()
        .any(|peak| (peak.rt - rt).abs() <= RT_TOLERANCE)
}

#[test]
fn real_targets_match_expected_verdicts() {
    let cases = group("real/");
    assert_eq!(cases.len(), 8, "expected 8 real targets");
    for (id, eic) in cases {
        let rt = param_f64(&eic, "rt");
        let threshold = param_f64(&eic, "threshold");
        let expect_keep = param_bool(&eic, "expect_keep");
        assert!(
            has_peak_near(&eic, rt, 0.0),
            "{id} must be detected without filter"
        );
        assert_eq!(
            has_peak_near(&eic, rt, threshold),
            expect_keep,
            "{id} at {threshold}"
        );
    }
}
