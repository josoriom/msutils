use crate::utilities::{
    find_peaks::{FindPeaksOptions, find_peaks},
    structs::{DataXY, Peak, Roi},
};

const LOCAL_WINDOW_MINUTES: f64 = 2.0;
const DEFAULT_MAX_DRIFT: f64 = 0.5;

pub fn get_peak(data: &DataXY, roi: &Roi, options: Option<FindPeaksOptions>) -> Option<Peak> {
    let target = roi.rt;

    let local;
    let section: &DataXY = if target.is_finite() {
        local = crop_around(data, target, LOCAL_WINDOW_MINUTES);
        &local
    } else {
        data
    };

    if section.x.is_empty() {
        return None;
    }

    let peaks = find_peaks(section, options);
    let max_drift = if roi.range.is_finite() && roi.range > 0.0 {
        roi.range
    } else {
        DEFAULT_MAX_DRIFT
    };

    pick_peak(&peaks, target, max_drift).copied()
}

fn position_weight(peak_rt: f64, target: f64, max_drift: f64) -> f64 {
    let drift = (peak_rt - target) / max_drift;
    1.0 - drift * drift
}

fn pick_peak(peaks: &[Peak], target: f64, max_drift: f64) -> Option<&Peak> {
    let mut best_peak: Option<&Peak> = None;
    let mut best_score = 0.0;
    for peak in peaks {
        if !peak.rt.is_finite() {
            continue;
        }
        let weight = position_weight(peak.rt, target, max_drift);
        if weight <= 0.0 {
            continue;
        }
        let score = peak.intensity * weight;
        if best_peak.is_none() || score > best_score {
            best_peak = Some(peak);
            best_score = score;
        }
    }
    best_peak
}

fn crop_around(data: &DataXY, center: f64, range: f64) -> DataXY {
    if data.x.is_empty() {
        return DataXY {
            x: Vec::new(),
            y: Vec::new(),
        };
    }
    let eic_min = data.x.iter().copied().fold(f64::INFINITY, f64::min);
    let eic_max = data.x.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let from = (center - range).max(eic_min);
    let to = (center + range).min(eic_max);
    let mut x = Vec::new();
    let mut y = Vec::new();
    for (rt, value) in data.x.iter().zip(data.y.iter()) {
        if *rt >= from && *rt <= to {
            x.push(*rt);
            y.push(*value);
        }
    }
    DataXY { x, y }
}

#[cfg(test)]
mod tests {
    use crate::utilities::{
        calculate_baseline::BaselineOptions,
        find_peaks::{FindPeaksOptions, PeakFilter},
        get_peak::get_peak,
        structs::{DataXY, Roi},
    };

    fn bell(rt: f64, center: f64, width: f64, height: f64) -> f64 {
        height * (-0.5 * ((rt - center) / width).powi(2)).exp()
    }

    fn default_options() -> FindPeaksOptions {
        FindPeaksOptions {
            boundaries: Some(Default::default()),
            filter: Some(PeakFilter {
                min_integral: None,
                min_intensity: Some(500.0),
                min_peak_width_points: Some(3),
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
            artifact_filter: Some(Default::default()),
        }
    }

    fn build_eic(
        from: f64,
        to: f64,
        points: usize,
        mut height_at: impl FnMut(f64) -> f64,
    ) -> DataXY {
        let mut x = Vec::with_capacity(points);
        let mut y = Vec::with_capacity(points);
        for index in 0..points {
            let rt = from + (to - from) * index as f64 / (points - 1) as f64;
            x.push(rt);
            y.push(height_at(rt).max(0.0));
        }
        DataXY { x, y }
    }

    #[test]
    fn finds_strong_peak_when_window_has_local_baseline() {
        let center = 28.6;
        let grass = [26.8, 27.2, 27.6, 29.6, 30.0, 30.4];
        let eic = build_eic(26.5, 30.7, 420, |rt| {
            let mut height = bell(rt, center, 0.05, 350_000.0);
            for grass_center in grass {
                height += bell(rt, grass_center, 0.02, 3_000.0);
            }
            height
        });
        let roi = Roi {
            rt: center,
            range: 1.0,
        };
        let peak =
            get_peak(&eic, &roi, Some(default_options())).expect("strong peak should be found");
        assert!((peak.rt - center).abs() < 0.05);
        assert!(peak.intensity > 100_000.0);
    }

    #[test]
    fn finds_isolated_peak_without_local_baseline_1() {
        let center = 28.6;
        let eic = build_eic(28.1, 29.1, 120, |rt| bell(rt, center, 0.05, 350_000.0));
        let roi = Roi {
            rt: center,
            range: 1.0,
        };
        let peak = get_peak(&eic, &roi, Some(default_options()))
            .expect("isolated peak on a flat baseline should be found");
        assert!((peak.rt - center).abs() < 0.05);
        assert!(peak.intensity > 100_000.0);
    }
}
