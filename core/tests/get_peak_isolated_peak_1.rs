// Bug #1 guard: a lone clean peak on a flat baseline must be found. The bug was that
// find_noise_level returned the peak's own apex as the noise floor, so get_peak rejected it.
// Fully synthetic (no fixture) so it runs identically here and in quantion/core/tests/.

use quantion::utilities::{
    find_noise_level::find_noise_level,
    find_peaks::{ArtifactFilter, FindPeaksOptions, PeakFilter},
    get_peak::get_peak,
    shape_filter::PeakShape,
    structs::{DataXY, Roi},
};

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
        artifact_filter: Some(ArtifactFilter {
            min_r2: 0.0,
            shape: PeakShape::EMG,
        }),
        ..Default::default()
    }
}

fn lone_peak_eic(points: usize, center_index: usize, apex: f64, width_points: f64) -> DataXY {
    let x: Vec<f64> = (0..points).map(|i| 28.0 + 0.01 * i as f64).collect();
    let center = x[center_index];
    let y: Vec<f64> = x
        .iter()
        .map(|&t| {
            let distance = (t - center) / (0.01 * width_points);
            (apex * (-0.5 * distance * distance).exp()).max(0.0)
        })
        .collect();
    DataXY { x, y }
}

#[test]
fn finds_isolated_lone_peak_1() {
    let eic = lone_peak_eic(200, 100, 350_000.0, 4.0);
    let center = eic.x[100];
    let peak = get_peak(&eic, &Roi::peak(center, 1.0), Some(options()));
    assert!(
        peak.intensity > 0.0,
        "get_peak must find a lone 350k peak on a flat baseline"
    );
    assert!(
        (peak.rt - center).abs() < 0.1,
        "found rt {} far from {center}",
        peak.rt
    );
    assert!(
        peak.intensity > 100_000.0,
        "intensity {} too low",
        peak.intensity
    );
}

#[test]
fn find_noise_level_must_not_report_lone_peak_as_noise_1() {
    let apex = 300_000.0;
    let mut signal = vec![0.0f64; 200];
    for (index, value) in signal.iter_mut().enumerate() {
        let distance = (index as f64 - 100.0) / 4.0;
        *value = apex * (-0.5 * distance * distance).exp();
    }
    let noise = find_noise_level(&signal).intensity;
    assert!(
        noise < 0.1 * apex,
        "find_noise_level reported noise {noise} for a lone {apex}-tall peak on a flat baseline; \
         the floor should be far below the peak, not equal to it"
    );
}
