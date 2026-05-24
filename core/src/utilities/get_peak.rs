use crate::utilities::find_peaks::{FindPeaksOptions, find_peaks};
use crate::utilities::structs::{DataXY, Peak, Roi};

pub fn get_peak(data: &DataXY, roi: &Roi, options: Option<FindPeaksOptions>) -> Option<Peak> {
    let peaks = find_peaks(data, options);

    if peaks.is_empty() {
        return None;
    }

    let target = roi.rt;
    let w = if roi.half_width.is_finite() && roi.half_width > 0.0 {
        roi.half_width
    } else {
        0.0
    };

    if w > 0.0 {
        let lo = target - w;
        let hi = target + w;
        let mut best: Option<&Peak> = None;
        for p in &peaks {
            if p.rt.is_finite() && p.rt >= lo && p.rt <= hi {
                match best {
                    None => best = Some(p),
                    Some(b) => {
                        let db = (b.rt - target).abs();
                        let dp = (p.rt - target).abs();
                        if dp < db || ((dp - db).abs() <= f64::EPSILON && p.intensity > b.intensity)
                        {
                            best = Some(p);
                        }
                    }
                }
            }
        }

        return best.copied();
    }

    let mut best_any = &peaks[0];
    for p in &peaks {
        let db = (best_any.rt - target).abs();
        let dp = (p.rt - target).abs();
        if dp < db || ((dp - db).abs() <= f64::EPSILON && p.intensity > best_any.intensity) {
            best_any = p;
        }
    }
    Some(*best_any)
}
