use std::{cmp::Ordering, sync::Arc};

use ionic::ion::{Ion, IonError};
use ionic::mzml::structs::MzML;
use ionic::{ScanSource, ScanSummary};
use serde::Serialize;

use crate::utilities::structs::{FromTo, Peak, ser_finite_f64};

pub(crate) const MS1_LEVEL: u8 = 1;

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

#[derive(Clone, Copy, Debug, Default)]
pub enum TimeUnit {
    Seconds,
    #[default]
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

pub enum EicReader<'a> {
    Ion(&'a mut Ion),
    Mzml(&'a mut MzML),
}

#[derive(Clone, Copy, Debug)]
pub struct ScanTime {
    pub index: usize,
    pub rt: f64,
}

#[derive(Clone, Debug)]
pub enum FastError {
    InvalidRequest,
    UnsupportedBackend,
    MissingBounds,
    BadBoundsChecksum,
    MalformedBounds(String),
    ReadFailed(String),
}

impl std::fmt::Display for FastError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FastError::InvalidRequest => write!(f, "invalid EIC request"),
            FastError::UnsupportedBackend => write!(f, "unsupported file backend"),
            FastError::MissingBounds => write!(
                f,
                "spectrum has no segment bounds (A3); re-encode the file with the current Ionic"
            ),
            FastError::BadBoundsChecksum => write!(f, "spectrum bounds (A3) checksum failed"),
            FastError::MalformedBounds(reason) => {
                write!(f, "spectrum bounds (A3) malformed: {}", reason)
            }
            FastError::ReadFailed(reason) => write!(f, "read failed: {}", reason),
        }
    }
}

impl std::error::Error for FastError {}

impl From<IonError> for FastError {
    fn from(error: IonError) -> Self {
        match error {
            IonError::MissingSpectrumBounds => FastError::MissingBounds,
            IonError::BadSpectrumBoundsChecksum => FastError::BadBoundsChecksum,
            IonError::MalformedSpectrumBounds(reason) => FastError::MalformedBounds(reason),
            other => FastError::ReadFailed(other.to_string()),
        }
    }
}

#[derive(Clone, Debug)]
pub struct EicRequest {
    pub target_mz: f64,
    pub rt_from: f64,
    pub rt_to: f64,
    pub ppm_tolerance: f64,
    pub mz_tolerance: f64,
    pub time_unit: TimeUnit,
    pub ms_level: u8,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct SpectrumSummary {
    #[serde(serialize_with = "ser_finite_f64")]
    pub rt_seconds: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    pub base_peak_mz: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    pub selected_ion_mz: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    pub base_peak_int: f64,
    #[serde(serialize_with = "ser_finite_f64")]
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
    pub fn from_summary(s: &ScanSummary) -> Self {
        Self {
            rt_seconds: s.rt * 60.0,
            base_peak_mz: s.base_peak_mz,
            selected_ion_mz: s.selected_ion_mz,
            base_peak_int: s.base_peak_int,
            total_ion_current: s.total_ion_current,
            ms_level: s.ms_level,
            polarity: s.polarity,
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct CentroidScan {
    #[serde(serialize_with = "ser_finite_f64")]
    pub rt: f64,
    pub mz: Arc<[f64]>,
    pub intensity: Arc<[f64]>,
    pub metadata: SpectrumSummary,
}

pub fn get_scan_times(
    reader: &mut EicReader,
    rt_from: f64,
    rt_to: f64,
    ms_level: u8,
) -> Vec<ScanTime> {
    let mut result = Vec::new();
    let rt_min = rt_from.min(rt_to);
    let rt_max = rt_from.max(rt_to);

    match reader {
        EicReader::Ion(ion) => {
            ion.for_each_summary(&mut |idx, summary| {
                if summary.rt >= rt_min && summary.rt <= rt_max && summary.ms_level == ms_level {
                    result.push(ScanTime {
                        index: idx,
                        rt: summary.rt,
                    });
                }
            });
        }
        EicReader::Mzml(mzml) => {
            mzml.for_each_summary(&mut |idx, summary| {
                if summary.rt >= rt_min && summary.rt <= rt_max && summary.ms_level == ms_level {
                    result.push(ScanTime {
                        index: idx,
                        rt: summary.rt,
                    });
                }
            });
        }
    }
    result
}

pub fn read_mz_window(
    reader: &mut EicReader,
    scan_index: usize,
    mz_from: f64,
    mz_to: f64,
    mz_out: &mut Vec<f64>,
    intensity_out: &mut Vec<f64>,
) -> Result<(), FastError> {
    mz_out.clear();
    intensity_out.clear();

    match reader {
        EicReader::Ion(ion) => {
            let window = ion
                .read_spectrum_mz_window_strict(scan_index, mz_from, mz_to)
                .map_err(FastError::from)?;
            mz_out.extend_from_slice(&window.mz);
            intensity_out.extend_from_slice(&window.intensity);
            Ok(())
        }
        EicReader::Mzml(mzml) => {
            let mut temp_mz = Vec::new();
            let mut temp_intensity = Vec::new();

            if !mzml.load_scan(scan_index, &mut temp_mz, &mut temp_intensity) {
                return Err(FastError::ReadFailed("failed to load scan".to_string()));
            }

            for (mz, intensity) in temp_mz.iter().zip(temp_intensity.iter()) {
                if *mz >= mz_from && *mz <= mz_to {
                    mz_out.push(*mz);
                    intensity_out.push(*intensity);
                }
            }
            Ok(())
        }
    }
}

pub fn calculate_eic(
    reader: &mut EicReader,
    target_mz: f64,
    time_range: FromTo,
    options: EicOptions,
) -> Result<Eic, FastError> {
    if !target_mz.is_finite() || target_mz <= 0.0 {
        return Err(FastError::InvalidRequest);
    }

    let tolerance = mz_tolerance_for(target_mz, options);
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(FastError::InvalidRequest);
    }

    let mz_lower = target_mz - tolerance;
    let mz_upper = target_mz + tolerance;

    let rt_min = options
        .time_unit
        .to_minutes(time_range.from.min(time_range.to));
    let rt_max = options
        .time_unit
        .to_minutes(time_range.from.max(time_range.to));

    let scan_times = get_scan_times(reader, rt_min, rt_max, MS1_LEVEL);

    if scan_times.is_empty() {
        return Ok(Eic::empty());
    }

    let mut x = Vec::new();
    let mut y = Vec::new();
    let mut mz_buf = Vec::new();
    let mut intensity_buf = Vec::new();

    for scan_time in scan_times {
        x.push(scan_time.rt);
        read_mz_window(
            reader,
            scan_time.index,
            mz_lower,
            mz_upper,
            &mut mz_buf,
            &mut intensity_buf,
        )?;

        let intensity_sum: f64 = intensity_buf.iter().sum();
        y.push(intensity_sum);
    }

    Ok(Eic { x, y })
}

pub fn calculate_eic_from_scan_source(
    source: &mut impl ScanSource,
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
    source.for_each_in_range(rt_min, rt_max, MS1_LEVEL, |summary, mz, intensity| {
        x.push(summary.rt);
        y.push(summed_intensity_in_window(
            mz, intensity, mz_lower, mz_upper,
        ));
    });
    Eic { x, y }
}

pub fn get_eic_for_mz(
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

#[derive(Debug, Clone, Copy)]
pub enum ScanQuery {
    RtRange(FromTo),
    ClosestRt(f64),
    MzRange(FromTo),
    ClosestMz(f64),
}

pub fn get_scans(
    source: &mut impl ScanSource,
    query: ScanQuery,
    time_unit: TimeUnit,
    ms_level: u8,
) -> (Vec<f64>, Vec<CentroidScan>) {
    match query {
        ScanQuery::RtRange(range) => get_by_rt_range(source, range, time_unit, ms_level),
        ScanQuery::ClosestRt(rt) => get_closest_by_rt(source, time_unit.to_minutes(rt), ms_level),
        ScanQuery::MzRange(range) => get_by_mz_range(source, range, ms_level),
        ScanQuery::ClosestMz(mz) => get_closest_by_mz(source, mz, ms_level),
    }
}

fn get_by_rt_range(
    source: &mut impl ScanSource,
    range: FromTo,
    time_unit: TimeUnit,
    ms_level: u8,
) -> (Vec<f64>, Vec<CentroidScan>) {
    let rt_min = time_unit.to_minutes(range.from.min(range.to));
    let rt_max = time_unit.to_minutes(range.from.max(range.to));
    let mut candidates: Vec<(usize, ScanSummary)> = Vec::new();
    source.for_each_summary(&mut |index, summary| {
        if ms_level != 0 && summary.ms_level != ms_level {
            return;
        }
        if !summary.rt.is_finite() || summary.rt < rt_min || summary.rt > rt_max {
            return;
        }
        candidates.push((index, summary));
    });
    candidates.sort_unstable_by(|a, b| a.1.rt.partial_cmp(&b.1.rt).unwrap_or(Ordering::Equal));
    load_selected(source, candidates)
}

fn get_closest_by_rt(
    source: &mut impl ScanSource,
    target_rt: f64,
    ms_level: u8,
) -> (Vec<f64>, Vec<CentroidScan>) {
    let mut by_dist: Vec<(f64, usize, ScanSummary)> = Vec::new();
    source.for_each_summary(&mut |index, summary| {
        if ms_level != 0 && summary.ms_level != ms_level {
            return;
        }
        if !summary.rt.is_finite() {
            return;
        }
        by_dist.push(((summary.rt - target_rt).abs(), index, summary));
    });
    by_dist.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal));
    load_first(source, by_dist.into_iter().map(|(_, i, s)| (i, s)))
}

fn get_by_mz_range(
    source: &mut impl ScanSource,
    range: FromTo,
    ms_level: u8,
) -> (Vec<f64>, Vec<CentroidScan>) {
    let mz_min = range.from.min(range.to);
    let mz_max = range.from.max(range.to);
    let mut candidates: Vec<(usize, ScanSummary)> = Vec::new();
    source.for_each_summary(&mut |index, summary| {
        if ms_level != 0 && summary.ms_level != ms_level {
            return;
        }
        let mz = summary.selected_ion_mz;
        if !mz.is_finite() || mz < mz_min || mz > mz_max {
            return;
        }
        candidates.push((index, summary));
    });
    candidates.sort_unstable_by(|a, b| a.1.rt.partial_cmp(&b.1.rt).unwrap_or(Ordering::Equal));
    load_selected(source, candidates)
}

fn get_closest_by_mz(
    source: &mut impl ScanSource,
    target_mz: f64,
    ms_level: u8,
) -> (Vec<f64>, Vec<CentroidScan>) {
    let mut by_dist: Vec<(f64, usize, ScanSummary)> = Vec::new();
    source.for_each_summary(&mut |index, summary| {
        if ms_level != 0 && summary.ms_level != ms_level {
            return;
        }
        if !summary.selected_ion_mz.is_finite() {
            return;
        }
        by_dist.push(((summary.selected_ion_mz - target_mz).abs(), index, summary));
    });
    by_dist.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal));
    load_first(source, by_dist.into_iter().map(|(_, i, s)| (i, s)))
}

fn load_first(
    source: &mut impl ScanSource,
    candidates: impl Iterator<Item = (usize, ScanSummary)>,
) -> (Vec<f64>, Vec<CentroidScan>) {
    let mut mz_buf = Vec::new();
    let mut int_buf = Vec::new();
    for (index, summary) in candidates {
        if !source.load_scan(index, &mut mz_buf, &mut int_buf) {
            continue;
        }
        let len = mz_buf.len().min(int_buf.len());
        return (
            vec![summary.rt],
            vec![CentroidScan {
                rt: summary.rt,
                mz: Arc::from(&mz_buf[..len]),
                intensity: Arc::from(&int_buf[..len]),
                metadata: SpectrumSummary::from_summary(&summary),
            }],
        );
    }
    (Vec::new(), Vec::new())
}

fn load_selected(
    source: &mut impl ScanSource,
    candidates: Vec<(usize, ScanSummary)>,
) -> (Vec<f64>, Vec<CentroidScan>) {
    let mut mz_buf = Vec::new();
    let mut int_buf = Vec::new();
    let mut rts = Vec::with_capacity(candidates.len());
    let mut scans = Vec::with_capacity(candidates.len());
    for (index, summary) in candidates {
        if !source.load_scan(index, &mut mz_buf, &mut int_buf) {
            continue;
        }
        let len = mz_buf.len().min(int_buf.len());
        rts.push(summary.rt);
        scans.push(CentroidScan {
            rt: summary.rt,
            mz: Arc::from(&mz_buf[..len]),
            intensity: Arc::from(&int_buf[..len]),
            metadata: SpectrumSummary::from_summary(&summary),
        });
    }
    (rts, scans)
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

#[inline]
pub fn mz_tolerance_for(target_mz: f64, options: EicOptions) -> f64 {
    let ppm_window = if options.ppm_tolerance > 0.0 {
        (options.ppm_tolerance * 1e-6) * target_mz.abs()
    } else {
        0.0
    };
    ppm_window.max(options.mz_tolerance.max(0.0))
}

#[inline]
pub(crate) fn summed_intensity_in_window(
    mz: &[f64],
    intensity: &[f64],
    mz_lower: f64,
    mz_upper: f64,
) -> f64 {
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

#[cfg(test)]
mod tests {
    use super::*;
    use ionic::ScanSummary;

    struct MockSource {
        scans: Vec<(ScanSummary, Vec<f64>, Vec<f64>)>,
    }

    impl MockSource {
        fn new(scans: Vec<(ScanSummary, Vec<f64>, Vec<f64>)>) -> Self {
            Self { scans }
        }
    }

    impl ScanSource for MockSource {
        fn for_each_summary(&mut self, cb: &mut dyn FnMut(usize, ScanSummary)) {
            for (i, (summary, _, _)) in self.scans.iter().enumerate() {
                cb(i, *summary);
            }
        }

        fn load_scan(&mut self, index: usize, mz: &mut Vec<f64>, intensity: &mut Vec<f64>) -> bool {
            let Some((_, m, i)) = self.scans.get(index) else {
                return false;
            };
            *mz = m.clone();
            *intensity = i.clone();
            !m.is_empty()
        }
    }

    fn make_scan(
        rt: f64,
        ms_level: u8,
        selected_ion_mz: f64,
        mz: Vec<f64>,
        intensity: Vec<f64>,
    ) -> (ScanSummary, Vec<f64>, Vec<f64>) {
        (
            ScanSummary {
                rt,
                ms_level,
                polarity: 0,
                selected_ion_mz,
                base_peak_mz: f64::NAN,
                base_peak_int: f64::NAN,
                total_ion_current: f64::NAN,
            },
            mz,
            intensity,
        )
    }

    #[test]
    fn lower_bound_empty() {
        assert_eq!(lower_bound(&[], 5.0), 0);
    }

    #[test]
    fn lower_bound_before_all() {
        assert_eq!(lower_bound(&[1.0, 2.0, 3.0], 0.0), 0);
    }

    #[test]
    fn lower_bound_after_all() {
        assert_eq!(lower_bound(&[1.0, 2.0, 3.0], 4.0), 3);
    }

    #[test]
    fn lower_bound_exact_match() {
        assert_eq!(lower_bound(&[1.0, 2.0, 3.0], 2.0), 1);
    }

    #[test]
    fn upper_bound_empty() {
        assert_eq!(upper_bound(&[], 5.0), 0);
    }

    #[test]
    fn upper_bound_exact_match() {
        assert_eq!(upper_bound(&[1.0, 2.0, 3.0], 2.0), 2);
    }

    #[test]
    fn upper_bound_after_all() {
        assert_eq!(upper_bound(&[1.0, 2.0, 3.0], 5.0), 3);
    }

    #[test]
    fn summed_intensity_sums_window() {
        let mz = vec![100.0, 200.0, 300.0];
        let int = vec![10.0, 20.0, 30.0];
        assert_eq!(summed_intensity_in_window(&mz, &int, 150.0, 250.0), 20.0);
    }

    #[test]
    fn summed_intensity_empty_input() {
        assert_eq!(summed_intensity_in_window(&[], &[], 0.0, 100.0), 0.0);
    }

    #[test]
    fn summed_intensity_no_match_in_window() {
        let mz = vec![100.0, 200.0];
        let int = vec![10.0, 20.0];
        assert_eq!(summed_intensity_in_window(&mz, &int, 300.0, 400.0), 0.0);
    }

    #[test]
    fn summed_intensity_multiple_in_window() {
        let mz = vec![100.0, 150.0, 200.0, 250.0];
        let int = vec![10.0, 15.0, 20.0, 25.0];
        assert_eq!(summed_intensity_in_window(&mz, &int, 100.0, 200.0), 45.0);
    }

    #[test]
    fn mz_tolerance_ppm_dominates() {
        let opts = EicOptions {
            ppm_tolerance: 10.0,
            mz_tolerance: 0.001,
            ..Default::default()
        };
        let tol = mz_tolerance_for(500.0, opts);
        assert!((tol - 0.005).abs() < 1e-9);
    }

    #[test]
    fn mz_tolerance_abs_dominates() {
        let opts = EicOptions {
            ppm_tolerance: 1.0,
            mz_tolerance: 0.1,
            ..Default::default()
        };
        let tol = mz_tolerance_for(500.0, opts);
        assert!((tol - 0.1).abs() < 1e-9);
    }

    #[test]
    fn get_by_rt_range_returns_matching_scans() {
        let mut source = MockSource::new(vec![
            make_scan(1.0, 1, f64::NAN, vec![100.0], vec![10.0]),
            make_scan(2.0, 1, f64::NAN, vec![200.0], vec![20.0]),
            make_scan(3.0, 1, f64::NAN, vec![300.0], vec![30.0]),
        ]);
        let (rts, scans) = get_scans(
            &mut source,
            ScanQuery::RtRange(FromTo { from: 1.0, to: 2.5 }),
            TimeUnit::Minutes,
            1,
        );
        assert_eq!(rts, vec![1.0, 2.0]);
        assert_eq!(scans.len(), 2);
    }

    #[test]
    fn get_by_rt_range_filters_wrong_ms_level() {
        let mut source = MockSource::new(vec![
            make_scan(1.0, 1, f64::NAN, vec![100.0], vec![10.0]),
            make_scan(2.0, 2, f64::NAN, vec![200.0], vec![20.0]),
        ]);
        let (rts, scans) = get_scans(
            &mut source,
            ScanQuery::RtRange(FromTo { from: 0.0, to: 5.0 }),
            TimeUnit::Minutes,
            1,
        );
        assert_eq!(rts, vec![1.0]);
        assert_eq!(scans.len(), 1);
    }

    #[test]
    fn get_by_rt_range_empty_when_no_match() {
        let mut source =
            MockSource::new(vec![make_scan(5.0, 1, f64::NAN, vec![100.0], vec![10.0])]);
        let (rts, scans) = get_scans(
            &mut source,
            ScanQuery::RtRange(FromTo { from: 0.0, to: 2.0 }),
            TimeUnit::Minutes,
            1,
        );
        assert!(rts.is_empty());
        assert!(scans.is_empty());
    }

    #[test]
    fn get_closest_rt_picks_nearest_scan() {
        let mut source = MockSource::new(vec![
            make_scan(1.0, 1, f64::NAN, vec![100.0], vec![10.0]),
            make_scan(3.0, 1, f64::NAN, vec![300.0], vec![30.0]),
            make_scan(5.0, 1, f64::NAN, vec![500.0], vec![50.0]),
        ]);
        let (rts, scans) = get_scans(&mut source, ScanQuery::ClosestRt(2.8), TimeUnit::Minutes, 1);
        assert_eq!(rts, vec![3.0]);
        assert_eq!(scans.len(), 1);
    }

    #[test]
    fn get_closest_rt_falls_back_when_best_fails_to_load() {
        let mut source = MockSource::new(vec![
            make_scan(1.0, 1, f64::NAN, vec![], vec![]),
            make_scan(2.0, 1, f64::NAN, vec![200.0], vec![20.0]),
        ]);
        let (rts, scans) = get_scans(&mut source, ScanQuery::ClosestRt(1.0), TimeUnit::Minutes, 1);
        assert_eq!(rts, vec![2.0]);
        assert_eq!(scans.len(), 1);
    }

    #[test]
    fn get_closest_rt_empty_when_source_empty() {
        let mut source = MockSource::new(vec![]);
        let (rts, scans) = get_scans(&mut source, ScanQuery::ClosestRt(1.0), TimeUnit::Minutes, 1);
        assert!(rts.is_empty());
        assert!(scans.is_empty());
    }

    #[test]
    fn get_by_mz_range_returns_matching_scans() {
        let mut source = MockSource::new(vec![
            make_scan(1.0, 2, 500.0, vec![500.0], vec![100.0]),
            make_scan(2.0, 2, 600.0, vec![600.0], vec![200.0]),
            make_scan(3.0, 2, 700.0, vec![700.0], vec![300.0]),
        ]);
        let (rts, scans) = get_scans(
            &mut source,
            ScanQuery::MzRange(FromTo {
                from: 490.0,
                to: 650.0,
            }),
            TimeUnit::Minutes,
            2,
        );
        assert_eq!(rts.len(), 2);
        assert!(rts.contains(&1.0) && rts.contains(&2.0));
        assert_eq!(scans.len(), 2);
    }

    #[test]
    fn get_by_mz_range_skips_nan_selected_ion_mz() {
        let mut source = MockSource::new(vec![
            make_scan(1.0, 2, f64::NAN, vec![100.0], vec![10.0]),
            make_scan(2.0, 2, 600.0, vec![600.0], vec![200.0]),
        ]);
        let (rts, scans) = get_scans(
            &mut source,
            ScanQuery::MzRange(FromTo {
                from: 0.0,
                to: 1000.0,
            }),
            TimeUnit::Minutes,
            2,
        );
        assert_eq!(rts, vec![2.0]);
        assert_eq!(scans.len(), 1);
    }

    #[test]
    fn get_closest_mz_picks_nearest_selected_ion() {
        let mut source = MockSource::new(vec![
            make_scan(1.0, 2, 500.0, vec![500.0], vec![100.0]),
            make_scan(2.0, 2, 600.0, vec![600.0], vec![200.0]),
            make_scan(3.0, 2, 700.0, vec![700.0], vec![300.0]),
        ]);
        let (rts, scans) = get_scans(
            &mut source,
            ScanQuery::ClosestMz(610.0),
            TimeUnit::Minutes,
            2,
        );
        assert_eq!(rts, vec![2.0]);
        assert_eq!(scans[0].metadata.selected_ion_mz, 600.0);
    }

    #[test]
    fn get_closest_mz_falls_back_when_best_fails_to_load() {
        let mut source = MockSource::new(vec![
            make_scan(1.0, 2, 500.0, vec![], vec![]),
            make_scan(2.0, 2, 510.0, vec![510.0], vec![200.0]),
        ]);
        let (rts, scans) = get_scans(
            &mut source,
            ScanQuery::ClosestMz(500.0),
            TimeUnit::Minutes,
            2,
        );
        assert_eq!(rts, vec![2.0]);
        assert_eq!(scans.len(), 1);
    }

    #[test]
    fn get_closest_mz_empty_when_no_finite_selected_ion_mz() {
        let mut source =
            MockSource::new(vec![make_scan(1.0, 2, f64::NAN, vec![100.0], vec![10.0])]);
        let (rts, scans) = get_scans(
            &mut source,
            ScanQuery::ClosestMz(500.0),
            TimeUnit::Minutes,
            2,
        );
        assert!(rts.is_empty());
        assert!(scans.is_empty());
    }

    #[test]
    fn ms_level_zero_allows_all_levels() {
        let mut source = MockSource::new(vec![
            make_scan(1.0, 1, f64::NAN, vec![100.0], vec![10.0]),
            make_scan(2.0, 2, f64::NAN, vec![200.0], vec![20.0]),
        ]);
        let (rts, scans) = get_scans(
            &mut source,
            ScanQuery::RtRange(FromTo { from: 0.0, to: 5.0 }),
            TimeUnit::Minutes,
            0,
        );
        assert_eq!(rts.len(), 2);
        assert_eq!(scans.len(), 2);
    }
}
