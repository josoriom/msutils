use std::{cmp::Ordering, sync::Arc};

use ionic::{
    BinaryData, BinaryDataArray, BinaryDataArrayList, MzML, PrecursorList, ScanList, ScanMeta,
    Spectrum, SpectrumSource,
};
use serde::Serialize;

use crate::utilities::structs::{FromTo, Peak};

const MS1_LEVEL: u8 = 1;
const ACC_MZ_ARRAY: &str = "MS:1000514";
const ACC_INTENSITY_ARRAY: &str = "MS:1000515";
const ACC_SCAN_START_TIME: &str = "MS:1000016";
const ACC_MS_LEVEL: &str = "MS:1000511";
const ACC_BASE_PEAK_MZ: &str = "MS:1000504";
const ACC_BASE_PEAK_INT: &str = "MS:1000505";
const ACC_TOTAL_ION_CURRENT: &str = "MS:1000285";
const ACC_SELECTED_ION_MZ: &str = "MS:1000744";
const ACC_POSITIVE_SCAN: &str = "MS:1000130";
const ACC_NEGATIVE_SCAN: &str = "MS:1000129";
const UO_MIN: &str = "UO:0000031";
const UO_SEC: &str = "UO:0000010";
const UO_MS: &str = "UO:0000028";

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

pub fn collect_scans_with_metadata(
    mzml: &MzML,
    time_range: FromTo,
    time_unit: TimeUnit,
    ms_level: u8,
) -> (Vec<f64>, Vec<CentroidScan>) {
    let spectra = all_spectra(mzml);
    let is_single_point = (time_range.from - time_range.to).abs() < f64::EPSILON;
    let window = time_window_in_minutes(time_range, time_unit);
    let target_rt = time_unit.to_minutes(time_range.from);
    let mut scans: Vec<CentroidScan> = Vec::with_capacity(spectra.len() / 2);

    let has_filter = mzml.filter_record.len() == spectra.len();

    for (i, spectrum) in spectra.iter().enumerate() {
        if !is_ms_level(spectrum, ms_level) {
            continue;
        }
        let Some(rt) = spectrum_retention_time_minutes(spectrum) else {
            continue;
        };
        if !is_single_point && !window.contains(rt) {
            continue;
        }
        let Some((masses, intensities)) = extract_mz_and_intensity(spectrum) else {
            continue;
        };

        let summary = if has_filter {
            let r = &mzml.filter_record[i];
            SpectrumSummary {
                rt_seconds: r.rt_seconds,
                base_peak_mz: r.base_peak_mz,
                selected_ion_mz: r.selected_ion_mz,
                base_peak_int: r.base_peak_int,
                total_ion_current: r.total_ion_current,
                ms_level: r.ms_level,
                polarity: r.polarity,
            }
        } else {
            extract_spectrum_summary(rt * 60.0, spectrum)
        };

        scans.push(CentroidScan {
            rt,
            mz: Arc::from(masses.into_boxed_slice()),
            intensity: Arc::from(intensities.into_boxed_slice()),
            metadata: summary,
        });
    }

    scans.sort_unstable_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap_or(Ordering::Equal));
    if is_single_point {
        return closest_scan_to(scans, target_rt);
    }
    let retention_times = scans.iter().map(|s| s.rt).collect();
    (retention_times, scans)
}

fn extract_spectrum_summary(rt_seconds: f64, spectrum: &Spectrum) -> SpectrumSummary {
    let mut s = SpectrumSummary::unknown();
    s.rt_seconds = rt_seconds;
    s.ms_level = ms_level_of(spectrum).unwrap_or(0);

    for cv in &spectrum.cv_params {
        match cv.accession.as_deref() {
            Some(ACC_BASE_PEAK_MZ) => {
                s.base_peak_mz = cv
                    .value
                    .as_deref()
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(f64::NAN);
            }
            Some(ACC_BASE_PEAK_INT) => {
                s.base_peak_int = cv
                    .value
                    .as_deref()
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(f64::NAN);
            }
            Some(ACC_TOTAL_ION_CURRENT) => {
                s.total_ion_current = cv
                    .value
                    .as_deref()
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(f64::NAN);
            }
            Some(ACC_POSITIVE_SCAN) => s.polarity = 1,
            Some(ACC_NEGATIVE_SCAN) => s.polarity = 2,
            _ => {}
        }
    }

    if let Some(pl) = spectrum_precursor_list(spectrum) {
        if let Some(precursor) = pl.precursors.first() {
            if let Some(sil) = &precursor.selected_ion_list {
                if let Some(si) = sil.selected_ions.first() {
                    for cv in &si.cv_params {
                        if cv.accession.as_deref() == Some(ACC_SELECTED_ION_MZ) {
                            s.selected_ion_mz = cv
                                .value
                                .as_deref()
                                .and_then(|v| v.parse().ok())
                                .unwrap_or(f64::NAN);
                        }
                    }
                }
            }
        }
    }

    s
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

struct TimeWindow {
    start: f64,
    end: f64,
}

impl TimeWindow {
    fn contains(&self, rt: f64) -> bool {
        rt >= self.start && rt <= self.end
    }
}

fn time_window_in_minutes(time_range: FromTo, time_unit: TimeUnit) -> TimeWindow {
    TimeWindow {
        start: time_unit.to_minutes(time_range.from.min(time_range.to)),
        end: time_unit.to_minutes(time_range.from.max(time_range.to)),
    }
}

fn all_spectra(mzml: &MzML) -> &[Spectrum] {
    mzml.run
        .spectrum_list
        .as_ref()
        .map(|l| l.spectra.as_slice())
        .unwrap_or(&[])
}

fn is_ms_level(spectrum: &Spectrum, level: u8) -> bool {
    ms_level_of(spectrum) == Some(level)
}

fn ms_level_of(spectrum: &Spectrum) -> Option<u8> {
    if let Some(level) = spectrum.ms_level {
        return u8::try_from(level).ok();
    }
    spectrum
        .cv_params
        .iter()
        .find(|p| p.accession.as_deref() == Some(ACC_MS_LEVEL))
        .and_then(|p| p.value.as_deref()?.parse::<u8>().ok())
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
fn spectrum_retention_time_minutes(spectrum: &Spectrum) -> Option<f64> {
    spectrum_scan_list(spectrum).and_then(retention_time_from_scan_list)
}

fn retention_time_from_scan_list(scan_list: &ScanList) -> Option<f64> {
    for scan in &scan_list.scans {
        for param in &scan.cv_params {
            if param.accession.as_deref() == Some(ACC_SCAN_START_TIME) {
                let raw = param.value.as_deref()?.parse::<f64>().ok()?;
                return convert_rt_to_minutes(
                    raw,
                    param.unit_accession.as_deref(),
                    param.unit_name.as_deref(),
                );
            }
        }
    }
    for param in &scan_list.cv_params {
        if param.accession.as_deref() == Some(ACC_SCAN_START_TIME) {
            let raw = param.value.as_deref()?.parse::<f64>().ok()?;
            return convert_rt_to_minutes(
                raw,
                param.unit_accession.as_deref(),
                param.unit_name.as_deref(),
            );
        }
    }
    None
}

#[inline]
fn convert_rt_to_minutes(
    raw: f64,
    unit_accession: Option<&str>,
    unit_name: Option<&str>,
) -> Option<f64> {
    if !raw.is_finite() {
        return None;
    }
    match unit_accession {
        Some(UO_MIN) => Some(raw),
        Some(UO_SEC) => Some(raw / 60.0),
        Some(UO_MS) => Some(raw / 60_000.0),
        _ => match unit_name {
            Some("minute") | Some("minutes") => Some(raw),
            Some("second") | Some("seconds") => Some(raw / 60.0),
            Some("millisecond") | Some("milliseconds") => Some(raw / 60_000.0),
            _ => None,
        },
    }
}

#[inline]
fn spectrum_scan_list(spectrum: &Spectrum) -> Option<&ScanList> {
    spectrum.scan_list.as_ref().or_else(|| {
        spectrum
            .spectrum_description
            .as_ref()
            .and_then(|d| d.scan_list.as_ref())
    })
}

#[inline]
fn spectrum_precursor_list(spectrum: &Spectrum) -> Option<&PrecursorList> {
    spectrum.precursor_list.as_ref().or_else(|| {
        spectrum
            .spectrum_description
            .as_ref()
            .and_then(|d| d.precursor_list.as_ref())
    })
}

fn extract_mz_and_intensity(spectrum: &Spectrum) -> Option<(Vec<f64>, Vec<f64>)> {
    let (mz_array, intensity_array) = spectrum_xy(spectrum)?;
    let mz_binary = mz_array.binary.as_ref()?;
    let intensity_binary = intensity_array.binary.as_ref()?;
    let count = bin_len(mz_binary).min(bin_len(intensity_binary));
    if count == 0 {
        return None;
    }
    let mut masses: Vec<f64> = Vec::with_capacity(count);
    let mut intensities: Vec<f64> = Vec::with_capacity(count);
    let mut last_mass = f64::NEG_INFINITY;
    let mut is_sorted = true;
    match (mz_binary, intensity_binary) {
        (BinaryData::F64(mz), BinaryData::F64(int)) => {
            for (&mass, &intensity) in mz.iter().zip(int.iter()).take(count) {
                record_scan_point(
                    &mut masses,
                    &mut intensities,
                    mass,
                    intensity,
                    &mut last_mass,
                    &mut is_sorted,
                );
            }
        }
        (BinaryData::F64(mz), BinaryData::F32(int)) => {
            for (&mass, &intensity) in mz.iter().zip(int.iter()).take(count) {
                record_scan_point(
                    &mut masses,
                    &mut intensities,
                    mass,
                    intensity as f64,
                    &mut last_mass,
                    &mut is_sorted,
                );
            }
        }
        (BinaryData::F32(mz), BinaryData::F64(int)) => {
            for (&mass, &intensity) in mz.iter().zip(int.iter()).take(count) {
                record_scan_point(
                    &mut masses,
                    &mut intensities,
                    mass as f64,
                    intensity,
                    &mut last_mass,
                    &mut is_sorted,
                );
            }
        }
        (BinaryData::F32(mz), BinaryData::F32(int)) => {
            for (&mass, &intensity) in mz.iter().zip(int.iter()).take(count) {
                record_scan_point(
                    &mut masses,
                    &mut intensities,
                    mass as f64,
                    intensity as f64,
                    &mut last_mass,
                    &mut is_sorted,
                );
            }
        }
        _ => {
            for index in 0..count {
                let mass = bin_get_f64(mz_binary, index);
                let intensity = bin_get_f64(intensity_binary, index);
                record_scan_point(
                    &mut masses,
                    &mut intensities,
                    mass,
                    intensity,
                    &mut last_mass,
                    &mut is_sorted,
                );
            }
        }
    }
    if masses.is_empty() {
        return None;
    }
    if !is_sorted {
        let (sm, si) = sort_masses_and_intensities(&masses, &intensities);
        return Some((sm, si));
    }
    Some((masses, intensities))
}

#[inline]
fn record_scan_point(
    masses: &mut Vec<f64>,
    intensities: &mut Vec<f64>,
    mass: f64,
    intensity: f64,
    last_mass: &mut f64,
    is_sorted: &mut bool,
) {
    if !mass.is_finite() || !intensity.is_finite() {
        return;
    }
    if mass < *last_mass {
        *is_sorted = false;
    }
    *last_mass = mass;
    masses.push(mass);
    intensities.push(intensity);
}

fn sort_masses_and_intensities(masses: &[f64], intensities: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let mut indices: Vec<usize> = (0..masses.len()).collect();
    indices.sort_unstable_by(|&a, &b| masses[a].partial_cmp(&masses[b]).unwrap_or(Ordering::Equal));
    (
        indices.iter().map(|&i| masses[i]).collect(),
        indices.iter().map(|&i| intensities[i]).collect(),
    )
}

#[inline]
fn spectrum_xy(spectrum: &Spectrum) -> Option<(&BinaryDataArray, &BinaryDataArray)> {
    let array_list = spectrum.binary_data_array_list.as_ref()?;
    Some((
        find_array_by_accession(array_list, ACC_MZ_ARRAY)?,
        find_array_by_accession(array_list, ACC_INTENSITY_ARRAY)?,
    ))
}

#[inline]
fn find_array_by_accession<'a>(
    array_list: &'a BinaryDataArrayList,
    accession: &str,
) -> Option<&'a BinaryDataArray> {
    array_list.binary_data_arrays.iter().find(|a| {
        a.cv_params
            .iter()
            .any(|p| p.accession.as_deref() == Some(accession))
    })
}

#[inline]
fn bin_len(binary: &BinaryData) -> usize {
    match binary {
        BinaryData::F64(v) => v.len(),
        BinaryData::F32(v) => v.len(),
        BinaryData::F16(v) => v.len(),
        BinaryData::I64(v) => v.len(),
        BinaryData::I32(v) => v.len(),
        BinaryData::I16(v) => v.len(),
    }
}

#[inline]
fn bin_get_f64(binary: &BinaryData, index: usize) -> f64 {
    match binary {
        BinaryData::F64(v) => v[index],
        BinaryData::F32(v) => v[index] as f64,
        BinaryData::F16(v) => v[index] as f64,
        BinaryData::I64(v) => v[index] as f64,
        BinaryData::I32(v) => v[index] as f64,
        BinaryData::I16(v) => v[index] as f64,
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

#[cfg(test)]
mod tests {
    use crate::OwnedIonReader;
    use crate::utilities::calculate_eic::{EicOptions, calculate_eic as calculate_eic_rs};
    use crate::utilities::structs::FromTo;
    use ionic::parse_mzml;
    use std::sync::Arc;

    const ION_BYTES: &[u8] = include_bytes!("./data/QCid_70-1000mz_10eV.ion");
    const MZML_BYTES: &[u8] = include_bytes!("./data/QCid_70-1000mz_10eV.mzML");

    fn eic_options() -> EicOptions {
        EicOptions {
            ppm_tolerance: 20.0,
            mz_tolerance: 0.005,
            ..Default::default()
        }
    }

    fn time_range() -> FromTo {
        FromTo {
            from: 0.0,
            to: 60.0,
        }
    }

    #[test]
    fn roundtrip_eic_lazy_equals_full() {
        let target_mz = 524.2;
        let mut mzml = parse_mzml(MZML_BYTES).unwrap();
        let eic_full = calculate_eic_rs(&mut mzml, target_mz, time_range(), eic_options());
        let arc: Arc<[u8]> = Arc::from(ION_BYTES);
        let mut reader = OwnedIonReader::new(arc).unwrap();
        let eic_lazy = calculate_eic_rs(&mut reader.inner, target_mz, time_range(), eic_options());
        assert_eq!(eic_full.x.len(), eic_lazy.x.len(), "x length mismatch");
        assert_eq!(eic_full.y.len(), eic_lazy.y.len(), "y length mismatch");
        for (a, b) in eic_full.x.iter().zip(eic_lazy.x.iter()) {
            assert_eq!(a.to_bits(), b.to_bits(), "x value mismatch: {a} vs {b}");
        }
        for (a, b) in eic_full.y.iter().zip(eic_lazy.y.iter()) {
            assert_eq!(a.to_bits(), b.to_bits(), "y value mismatch: {a} vs {b}");
        }
    }

    #[test]
    fn memory_comparison() {
        fn rss_kb() -> u64 {
            let pid = std::process::id();
            let out = std::process::Command::new("ps")
                .args(["-o", "rss=", "-p", &pid.to_string()])
                .output()
                .unwrap();
            String::from_utf8_lossy(&out.stdout)
                .trim()
                .parse()
                .unwrap_or(0)
        }
        let before = rss_kb();
        let mut mzml = parse_mzml(MZML_BYTES).unwrap();
        let _eic = calculate_eic_rs(&mut mzml, 524.2, time_range(), eic_options());
        eprintln!("FULL — RSS: {} KB", rss_kb() - before);
        drop(mzml);
        let before = rss_kb();
        let arc: Arc<[u8]> = Arc::from(ION_BYTES);
        let mut reader = OwnedIonReader::new(arc).unwrap();
        let _eic = calculate_eic_rs(&mut reader.inner, 524.2, time_range(), eic_options());
        eprintln!("LAZY — RSS: {} KB", rss_kb() - before);
    }
}
