mod helpers;

use std::path::{Path, PathBuf};

use helpers::{Eic, load_chromatograms_from, param_bool, param_f64};
use quantion::utilities::{
    find_noise_level::find_noise_level,
    find_peaks::{ArtifactFilter, FindPeaksOptions, PeakFilter},
    get_peak::get_peak,
    shape_filter::PeakShape,
    structs::{DataXY, Roi},
};

fn ion_path() -> PathBuf {
    if let Ok(path) = std::env::var("TEST_ION") {
        return PathBuf::from(path);
    }
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join("test.ion")
}

fn shape_from(eic: &Eic) -> PeakShape {
    match eic.params.get("opt.shape").map(String::as_str) {
        Some("Gaussian") => PeakShape::Gaussian,
        _ => PeakShape::EMG,
    }
}

fn options_from(eic: &Eic) -> FindPeaksOptions {
    FindPeaksOptions {
        filter: Some(PeakFilter {
            auto_noise: Some(param_bool(eic, "opt.auto_noise")),
            auto_baseline: Some(param_bool(eic, "opt.auto_baseline")),
            min_peak_width_points: Some(param_f64(eic, "opt.min_peak_width_points") as usize),
            min_intensity: Some(param_f64(eic, "opt.min_intensity")),
            min_snr: Some(param_f64(eic, "opt.min_snr")),
            ..Default::default()
        }),
        artifact_filter: Some(ArtifactFilter {
            min_r2: param_f64(eic, "opt.min_r2"),
            shape: shape_from(eic),
        }),
        ..Default::default()
    }
}

fn assert_get_peak(id: &str, eic: &Eic) {
    let description = eic.params.get("description").cloned().unwrap_or_default();
    let data = DataXY {
        x: eic.time.clone(),
        y: eic.intensity.clone(),
    };
    let target = param_f64(eic, "opt.target");
    let roi = Roi::peak(target, param_f64(eic, "opt.roi_half_width"));
    let peak = get_peak(&data, &roi, Some(options_from(eic)));
    let found = peak.intensity > 0.0;

    if param_bool(eic, "exp.must_find") {
        assert!(found, "[{id}] expected a peak but got none. {description}");
        let rt_within = param_f64(eic, "exp.rt_within");
        assert!(
            (peak.rt - target).abs() <= rt_within,
            "[{id}] found rt {} but expected within {} of {}",
            peak.rt,
            rt_within,
            target
        );
        let min_intensity = param_f64(eic, "exp.min_intensity");
        assert!(
            peak.intensity >= min_intensity,
            "[{id}] peak intensity {} below expected minimum {}",
            peak.intensity,
            min_intensity
        );
    } else {
        assert!(!found, "[{id}] expected no peak but one was found");
    }

    if let Some(fraction) = eic.params.get("exp.noise_below_apex_fraction") {
        let fraction: f64 = fraction.parse().expect("numeric noise fraction");
        let apex = eic.intensity.iter().cloned().fold(f64::MIN, f64::max);
        let noise = find_noise_level(&eic.intensity).intensity;
        assert!(
            noise < fraction * apex,
            "[{id}] find_noise_level reported noise {noise}, not below {fraction} x apex {apex}"
        );
    }
}

#[test]
fn bug_cases_hold_their_contract() {
    let path = ion_path();
    let cases = load_chromatograms_from(&path);
    assert!(!cases.is_empty(), "no cases found in {:?}", path);

    let mut tested = 0usize;
    let mut failures = Vec::new();
    for (id, eic) in &cases {
        if eic.params.get("kind").map(String::as_str) != Some("get_peak") {
            continue;
        }
        tested += 1;
        let result =
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| assert_get_peak(id, eic)));
        if result.is_err() {
            failures.push(id.clone());
        }
    }

    assert!(tested > 0, "no get_peak cases in {:?}", path);
    assert!(
        failures.is_empty(),
        "{}/{} get_peak cases violated their contract: {:?}",
        failures.len(),
        tested,
        failures
    );
}
