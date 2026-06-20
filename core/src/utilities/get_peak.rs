use crate::utilities::find_peaks::{FindPeaksOptions, find_peaks};
use crate::utilities::structs::{DataXY, Peak, Roi};

const LOCAL_WINDOW_MINUTES: f64 = 2.0;

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
    let half_width = if roi.half_width.is_finite() && roi.half_width > 0.0 {
        roi.half_width
    } else {
        0.0
    };

    find_closest_peak(&peaks, target, half_width).copied()
}

fn find_closest_peak(peaks: &[Peak], target: f64, half_width: f64) -> Option<&Peak> {
    let mut closest: Option<&Peak> = None;
    for peak in peaks {
        if half_width > 0.0 && (!peak.rt.is_finite() || (peak.rt - target).abs() > half_width) {
            continue;
        }
        match closest {
            None => closest = Some(peak),
            Some(best) => {
                let best_distance = (best.rt - target).abs();
                let peak_distance = (peak.rt - target).abs();
                let same_distance = (peak_distance - best_distance).abs() <= f64::EPSILON;
                if peak_distance < best_distance
                    || (same_distance && peak.intensity > best.intensity)
                {
                    closest = Some(peak);
                }
            }
        }
    }
    closest
}

fn crop_around(data: &DataXY, center: f64, half_width: f64) -> DataXY {
    if data.x.is_empty() {
        return DataXY {
            x: Vec::new(),
            y: Vec::new(),
        };
    }
    let eic_min = data.x.iter().copied().fold(f64::INFINITY, f64::min);
    let eic_max = data.x.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let from = (center - half_width).max(eic_min);
    let to = (center + half_width).min(eic_max);
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
