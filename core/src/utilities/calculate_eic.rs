use std::{cmp::Ordering, sync::Arc};

use ionic::{ScanMeta, SpectrumSource};
use serde::Serialize;

use crate::utilities::structs::{FromTo, Peak};

const MS1_LEVEL: u8 = 1;

#[derive(Clone, Copy, Debug)]
pub struct EicOptions {
    pub ppm_tolerance: f64,
    pub mz_tolerance: f64,
    pub time_unit: TimeUnit,
}

impl Default for EicOptions {
    fn default() -> Self {
        Self {
            ppm_tolerance: 20.0,
            mz_tolerance: 0.005,
            time_unit: TimeUnit::Minutes,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub enum TimeUnit {
    Seconds,
    Minutes,
}

impl TimeUnit {
    #[inline]
    pub fn to_minutes(self, value: f64) -> f64 {
        match self {
            TimeUnit::Seconds => value / 60.0,
            TimeUnit::Minutes => value,
        }
    }
}

impl Default for TimeUnit {
    fn default() -> Self {
        TimeUnit::Minutes
    }
}

pub struct Eic {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
}

impl Eic {
    fn empty() -> Self {
        Self {
            x: Vec::new(),
            y: Vec::new(),
        }
    }
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct SpectrumSummary {
    pub rt_seconds: f64,
    pub base_peak_mz: f64,
    pub selected_ion_mz: f64,
    pub base_peak_int: f64,
    pub total_ion_current: f64,
    pub ms_level: u8,
    pub polarity: u8,
}

impl SpectrumSummary {
    #[inline]
    pub fn unknown() -> Self {
        Self {
            rt_seconds: f64::NAN,
            base_peak_mz: f64::NAN,
            selected_ion_mz: f64::NAN,
            base_peak_int: f64::NAN,
            total_ion_current: f64::NAN,
            ms_level: 0,
            polarity: 0,
        }
    }

    #[inline]
    pub fn from_scan_meta(rt_seconds: f64, meta: &ScanMeta) -> Self {
        Self {
            rt_seconds,
            base_peak_mz: meta.base_peak_mz,
            selected_ion_mz: meta.selected_ion_mz,
            base_peak_int: meta.base_peak_int,
            total_ion_current: meta.total_ion_current,
            ms_level: meta.ms_level,
            polarity: meta.polarity,
        }
    }
}

#[derive(Clone, Debug)]
pub struct CentroidScan {
    pub rt: f64,
    pub mz: Arc<[f64]>,
    pub intensity: Arc<[f64]>,
    pub metadata: SpectrumSummary,
}

pub fn calculate_eic(
    source: &mut impl SpectrumSource,
    target_mass: f64,
    time_range: FromTo,
    options: EicOptions,
) -> Eic {
    if !target_mass.is_finite() || target_mass <= 0.0 {
        return Eic::empty();
    }
    let tolerance = mz_tolerance_for(target_mass, options);
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Eic::empty();
    }
    let mz_lower = target_mass - tolerance;
    let mz_upper = target_mass + tolerance;
    let rt_min = options
        .time_unit
        .to_minutes(time_range.from.min(time_range.to));
    let rt_max = options
        .time_unit
        .to_minutes(time_range.from.max(time_range.to));
    let mut x = Vec::new();
    let mut y = Vec::new();
    source.for_each_scan_in_range(
        rt_min,
        rt_max,
        MS1_LEVEL,
        &mut |rt, _meta, mz, intensity| {
            x.push(rt);
            y.push(summed_intensity_in_window(
                mz, intensity, mz_lower, mz_upper,
            ));
        },
    );
    Eic { x, y }
}

pub fn compute_eic_for_mz(
    scans: &[CentroidScan],
    scan_count: usize,
    target_mz: f64,
    options: EicOptions,
) -> Vec<f64> {
    if !target_mz.is_finite() || target_mz <= 0.0 {
        return vec![0.0; scan_count];
    }
    let tolerance = mz_tolerance_for(target_mz, options);
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return vec![0.0; scan_count];
    }
    let mz_lower = target_mz - tolerance;
    let mz_upper = target_mz + tolerance;
    let mut intensities = vec![0.0f64; scan_count];
    for (index, scan) in scans.iter().take(scan_count).enumerate() {
        intensities[index] =
            summed_intensity_in_window(&scan.mz, &scan.intensity, mz_lower, mz_upper);
    }
    intensities
}

pub fn collect_scans(
    source: &mut impl SpectrumSource,
    time_range: FromTo,
    time_unit: TimeUnit,
    ms_level: u8,
) -> (Vec<f64>, Vec<CentroidScan>) {
    let is_single_point = (time_range.from - time_range.to).abs() < f64::EPSILON;
    let target_rt = time_unit.to_minutes(time_range.from);
    let (rt_min, rt_max) = if is_single_point {
        (0.0_f64, f64::MAX)
    } else {
        (
            time_unit.to_minutes(time_range.from.min(time_range.to)),
            time_unit.to_minutes(time_range.from.max(time_range.to)),
        )
    };
    let mut scans = Vec::new();

    source.for_each_scan_in_range(rt_min, rt_max, ms_level, &mut |rt, meta, mz, intensity| {
        scans.push(CentroidScan {
            rt,
            mz: Arc::from(mz),
            intensity: Arc::from(intensity),
            metadata: SpectrumSummary::from_scan_meta(rt * 60.0, meta),
        });
    });

    scans.sort_unstable_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap_or(Ordering::Equal));
    if is_single_point {
        return closest_scan_to(scans, target_rt);
    }
    let retention_times = scans.iter().map(|s| s.rt).collect();
    (retention_times, scans)
}

pub fn with_eic_apex_intensity(rt: &[f64], y: &[f64], mut p: Peak) -> Peak {
    let apex_intensity = max_in_range(rt, y, p.from, p.to);
    if apex_intensity.is_finite() && apex_intensity > 0.0 {
        p.intensity = apex_intensity;
    }
    p
}

#[inline]
pub fn lower_bound(values: &[f64], target: f64) -> usize {
    let mut low = 0usize;
    let mut high = values.len();
    while low < high {
        let mid = (low + high) / 2;
        if values[mid] < target {
            low = mid + 1;
        } else {
            high = mid;
        }
    }
    low
}

#[inline]
pub fn upper_bound(values: &[f64], target: f64) -> usize {
    let mut low = 0usize;
    let mut high = values.len();
    while low < high {
        let mid = (low + high) / 2;
        if values[mid] <= target {
            low = mid + 1;
        } else {
            high = mid;
        }
    }
    low
}

fn closest_scan_to(mut scans: Vec<CentroidScan>, target_rt: f64) -> (Vec<f64>, Vec<CentroidScan>) {
    let closest_index = scans
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            (a.rt - target_rt)
                .abs()
                .partial_cmp(&(b.rt - target_rt).abs())
                .unwrap_or(Ordering::Equal)
        })
        .map(|(i, _)| i);
    match closest_index {
        Some(i) => {
            let s = scans.swap_remove(i);
            (vec![s.rt], vec![s])
        }
        None => (Vec::new(), Vec::new()),
    }
}

#[inline]
fn mz_tolerance_for(target_mz: f64, options: EicOptions) -> f64 {
    let ppm_window = if options.ppm_tolerance > 0.0 {
        (options.ppm_tolerance * 1e-6) * target_mz.abs()
    } else {
        0.0
    };
    ppm_window.max(options.mz_tolerance.max(0.0))
}

#[inline]
fn summed_intensity_in_window(mz: &[f64], intensity: &[f64], mz_lower: f64, mz_upper: f64) -> f64 {
    if mz.is_empty() || intensity.is_empty() || mz.len() != intensity.len() {
        return 0.0;
    }
    let start = lower_bound(mz, mz_lower);
    let end = upper_bound(mz, mz_upper);
    if end <= start {
        return 0.0;
    }
    unsafe { intensity.get_unchecked(start..end) }.iter().sum()
}

fn max_in_range(retention_times: &[f64], intensities: &[f64], from_rt: f64, to_rt: f64) -> f64 {
    let start = lower_bound(retention_times, from_rt);
    let end = upper_bound(retention_times, to_rt).min(intensities.len());
    if start >= intensities.len() || end <= start {
        return 0.0;
    }
    let mut maximum = 0.0f64;
    for &value in &intensities[start..end] {
        if value > maximum {
            maximum = value;
        }
    }
    maximum
}
