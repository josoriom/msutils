use quantion::utilities::{
    find_peaks::{ArtifactFilter, FindPeaksOptions, PeakFilter},
    get_peak::get_peak,
    shape_filter::PeakShape,
    structs::{DataXY, Roi},
};

const STEP_MINUTES: f64 = 0.006;

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

fn gaussian(time: f64, center: f64, height: f64, width: f64) -> f64 {
    let distance = (time - center) / width;
    height * (-0.5 * distance * distance).exp()
}

fn ripple(index: usize, amplitude: f64) -> f64 {
    let value = index.wrapping_mul(2654435761) % 1009;
    amplitude * (value as f64 / 1009.0)
}

fn build_eic(from: f64, to: f64, signal: impl Fn(usize, f64) -> f64) -> DataXY {
    let mut x = Vec::new();
    let mut y = Vec::new();
    let mut index = 0usize;
    let mut time = from;
    while time <= to {
        x.push(time);
        y.push(signal(index, time).max(0.0));
        time += STEP_MINUTES;
        index += 1;
    }
    DataXY { x, y }
}

fn find_at(eic: &DataXY, target: f64) -> f64 {
    let peak = get_peak(eic, &Roi::peak(target, 0.5), Some(options()));
    assert!(peak.intensity > 0.0, "expected a peak near rt {target}");
    assert!(
        (peak.rt - target).abs() < 0.1,
        "found rt {} but expected near {target}",
        peak.rt
    );
    peak.rt
}

#[test]
fn finds_peak_with_elevated_left_edge() {
    let eic = build_eic(26.0, 30.0, |index, time| {
        let peak = gaussian(time, 28.63, 300_000.0, 0.03);
        let left_tail = gaussian(time, 26.2, 150_000.0, 0.6);
        peak + left_tail + ripple(index, 20_000.0)
    });
    find_at(&eic, 28.63);
}

#[test]
fn finds_weak_peak_despite_distant_noise() {
    let eic = build_eic(0.4, 35.0, |index, time| {
        let peak = gaussian(time, 10.0, 250_000.0, 0.03);
        let distant = if (2.0..4.0).contains(&time) || (28.0..32.0).contains(&time) {
            ripple(index, 300_000.0)
        } else {
            0.0
        };
        peak + distant + ripple(index, 10_000.0)
    });
    find_at(&eic, 10.0);
}

#[test]
fn finds_peak_at_start_of_run() {
    let eic = build_eic(0.4, 35.0, |index, time| {
        gaussian(time, 0.5, 1_000_000.0, 0.03) + ripple(index, 20_000.0)
    });
    find_at(&eic, 0.5);
}

#[test]
fn finds_peak_at_end_of_run() {
    let eic = build_eic(0.4, 35.0, |index, time| {
        gaussian(time, 34.8, 1_000_000.0, 0.03) + ripple(index, 20_000.0)
    });
    find_at(&eic, 34.8);
}

#[test]
fn picks_peak_closest_to_target_over_taller_neighbor() {
    let eic = build_eic(26.0, 31.0, |index, time| {
        let taller_neighbor = gaussian(time, 28.4, 1_000_000.0, 0.04);
        let target = gaussian(time, 28.63, 300_000.0, 0.03);
        taller_neighbor + target + ripple(index, 20_000.0)
    });
    let peak = get_peak(&eic, &Roi::peak(28.63, 0.2), Some(options()));
    assert!(peak.intensity > 0.0, "expected a peak near rt 28.63");
    assert!(
        (peak.rt - 28.4).abs() > 0.1,
        "selected the taller neighbor at 28.4 instead of the target at 28.63",
    );
    assert!(
        (peak.rt - 28.63).abs() < 0.1,
        "found rt {} but expected near 28.63",
        peak.rt
    );
}

#[test]
fn finds_spiky_peak() {
    let eic = build_eic(28.0, 29.3, |index, time| {
        let peak = gaussian(time, 28.63, 250_000.0, 0.04);
        let spikes = ripple(index, 120_000.0);
        peak + spikes
    });
    find_at(&eic, 28.63);
}

#[test]
fn finds_peak_on_rising_baseline() {
    let eic = build_eic(0.4, 35.0, |index, time| {
        let baseline = 50_000.0 + 5_000.0 * time;
        let peak = gaussian(time, 15.0, 400_000.0, 0.03);
        baseline + peak + ripple(index, 15_000.0)
    });
    find_at(&eic, 15.0);
}
