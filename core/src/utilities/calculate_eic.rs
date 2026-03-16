use std::{cmp::Ordering, sync::Arc};

use octo::{BinaryData, BinaryDataArray, BinaryDataArrayList, MzML, Spectrum};
use serde::Serialize;

use crate::utilities::structs::{FromTo, Peak};

const MS1_LEVEL: u8 = 1;

const ACC_MZ_ARRAY: &str = "MS:1000514";
const ACC_INTENSITY_ARRAY: &str = "MS:1000515";
const ACC_SCAN_START_TIME: &str = "MS:1000016";
const ACC_MS_LEVEL: &str = "MS:1000511";

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
pub struct ScanMetadataEntry {
    pub section: &'static str,
    pub accession: Option<String>,
    pub name: String,
    pub value: String,
    pub unit_accession: Option<String>,
    pub unit_name: Option<String>,
}

impl ScanMetadataEntry {
    #[inline]
    fn from_cv_param(section: &'static str, param: &octo::CvParam) -> Option<Self> {
        let value = param.value.as_deref()?.trim();
        if value.is_empty() {
            return None;
        }
        Some(Self {
            section,
            accession: param.accession.clone(),
            name: param.name.clone(),
            value: value.to_owned(),
            unit_accession: param.unit_accession.clone(),
            unit_name: param.unit_name.clone(),
        })
    }
}

#[derive(Clone, Debug)]
pub struct CentroidScan {
    pub rt: f64,
    pub mz: Arc<[f64]>,
    pub intensity: Arc<[f64]>,
    pub metadata: Vec<ScanMetadataEntry>,
}

pub fn calculate_eic(
    mzml: &MzML,
    target_mass: f64,
    time_range: FromTo,
    options: EicOptions,
) -> Eic {
    let (retention_times, scans) =
        collect_scans(mzml, time_range, options.time_unit, MS1_LEVEL, false);

    if scans.is_empty() {
        return Eic::empty();
    }

    let scan_count = retention_times.len();

    Eic {
        x: retention_times,
        y: compute_eic_for_mz(&scans, scan_count, target_mass, options),
    }
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
    mzml: &MzML,
    time_range: FromTo,
    time_unit: TimeUnit,
    ms_level: u8,
    collect_metadata: bool,
) -> (Vec<f64>, Vec<CentroidScan>) {
    let spectra = all_spectra(mzml);
    let is_single_point = (time_range.from - time_range.to).abs() < f64::EPSILON;
    let window = time_window_in_minutes(time_range, time_unit);
    let target_rt = time_unit.to_minutes(time_range.from);

    let mut scans: Vec<CentroidScan> = Vec::with_capacity(spectra.len() / 2);

    for spectrum in spectra {
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

        let metadata = if collect_metadata {
            collect_scan_metadata(spectrum)
        } else {
            Vec::new()
        };

        scans.push(CentroidScan {
            rt,
            mz: Arc::from(masses.into_boxed_slice()),
            intensity: Arc::from(intensities.into_boxed_slice()),
            metadata,
        });
    }

    scans.sort_unstable_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap_or(Ordering::Equal));

    if is_single_point {
        return closest_scan_to(scans, target_rt);
    }

    let retention_times = scans.iter().map(|scan| scan.rt).collect();
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
        .map(|list| list.spectra.as_slice())
        .unwrap_or(&[])
}

fn is_ms_level(spectrum: &Spectrum, level: u8) -> bool {
    ms_level_of(spectrum) == Some(level)
}

fn ms_level_of(spectrum: &Spectrum) -> Option<u8> {
    if let Some(level) = spectrum.ms_level {
        return u8::try_from(level).ok();
    }
    for param in &spectrum.cv_params {
        if param.accession.as_deref() == Some(ACC_MS_LEVEL) {
            return param.value.as_deref().and_then(|v| v.parse::<u8>().ok());
        }
    }
    None
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
        .map(|(index, _)| index);

    match closest_index {
        Some(index) => {
            let scan = scans.swap_remove(index);
            (vec![scan.rt], vec![scan])
        }
        None => (Vec::new(), Vec::new()),
    }
}

#[inline]
fn spectrum_retention_time_minutes(spectrum: &Spectrum) -> Option<f64> {
    spectrum_scan_list(spectrum).and_then(retention_time_from_scan_list)
}

fn retention_time_from_scan_list(scan_list: &octo::ScanList) -> Option<f64> {
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
fn spectrum_scan_list(spectrum: &Spectrum) -> Option<&octo::ScanList> {
    spectrum.scan_list.as_ref().or_else(|| {
        spectrum
            .spectrum_description
            .as_ref()
            .and_then(|description| description.scan_list.as_ref())
    })
}

#[inline]
fn spectrum_precursor_list(spectrum: &Spectrum) -> Option<&octo::PrecursorList> {
    spectrum.precursor_list.as_ref().or_else(|| {
        spectrum
            .spectrum_description
            .as_ref()
            .and_then(|description| description.precursor_list.as_ref())
    })
}

#[inline]
fn spectrum_product_list(spectrum: &Spectrum) -> Option<&octo::ProductList> {
    spectrum.product_list.as_ref().or_else(|| {
        spectrum
            .spectrum_description
            .as_ref()
            .and_then(|description| description.product_list.as_ref())
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
        let (sorted_masses, sorted_intensities) =
            sort_masses_and_intensities(&masses, &intensities);
        return Some((sorted_masses, sorted_intensities));
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
    let sorted_masses = indices.iter().map(|&i| masses[i]).collect();
    let sorted_intensities = indices.iter().map(|&i| intensities[i]).collect();
    (sorted_masses, sorted_intensities)
}

#[inline]
fn spectrum_xy(spectrum: &Spectrum) -> Option<(&BinaryDataArray, &BinaryDataArray)> {
    let array_list = spectrum.binary_data_array_list.as_ref()?;
    let mz_array = find_array_by_accession(array_list, ACC_MZ_ARRAY)?;
    let intensity_array = find_array_by_accession(array_list, ACC_INTENSITY_ARRAY)?;
    Some((mz_array, intensity_array))
}

#[inline]
fn find_array_by_accession<'a>(
    array_list: &'a BinaryDataArrayList,
    accession: &str,
) -> Option<&'a BinaryDataArray> {
    array_list.binary_data_arrays.iter().find(|array| {
        array
            .cv_params
            .iter()
            .any(|param| param.accession.as_deref() == Some(accession))
    })
}

#[inline]
fn bin_len(binary: &BinaryData) -> usize {
    match binary {
        BinaryData::F64(values) => values.len(),
        BinaryData::F32(values) => values.len(),
        BinaryData::F16(values) => values.len(),
        BinaryData::I64(values) => values.len(),
        BinaryData::I32(values) => values.len(),
        BinaryData::I16(values) => values.len(),
    }
}

#[inline]
fn bin_get_f64(binary: &BinaryData, index: usize) -> f64 {
    match binary {
        BinaryData::F64(values) => values[index],
        BinaryData::F32(values) => values[index] as f64,
        BinaryData::F16(values) => values[index] as f64,
        BinaryData::I64(values) => values[index] as f64,
        BinaryData::I32(values) => values[index] as f64,
        BinaryData::I16(values) => values[index] as f64,
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

fn collect_scan_metadata(spectrum: &Spectrum) -> Vec<ScanMetadataEntry> {
    let mut metadata = Vec::with_capacity(18);

    push_metadata_value("spectrum", "id", Some(spectrum.id.as_str()), &mut metadata);
    push_valued_cv_params("spectrum", &spectrum.cv_params, &mut metadata);

    if let Some(scan_list) = spectrum_scan_list(spectrum) {
        push_valued_cv_params("scan_list", &scan_list.cv_params, &mut metadata);

        for scan in &scan_list.scans {
            push_valued_cv_params("scan", &scan.cv_params, &mut metadata);

            if let Some(scan_window_list) = scan.scan_window_list.as_ref() {
                for scan_window in &scan_window_list.scan_windows {
                    push_valued_cv_params("scan_window", &scan_window.cv_params, &mut metadata);
                }
            }
        }
    }

    if let Some(precursor_list) = spectrum_precursor_list(spectrum) {
        push_valued_cv_params("precursor_list", &precursor_list.cv_params, &mut metadata);

        for precursor in &precursor_list.precursors {
            if let Some(isolation_window) = precursor.isolation_window.as_ref() {
                push_valued_cv_params(
                    "isolation_window",
                    &isolation_window.cv_params,
                    &mut metadata,
                );
            }
            if let Some(selected_ion_list) = precursor.selected_ion_list.as_ref() {
                for selected_ion in &selected_ion_list.selected_ions {
                    push_valued_cv_params("selected_ion", &selected_ion.cv_params, &mut metadata);
                }
            }
            if let Some(activation) = precursor.activation.as_ref() {
                push_valued_cv_params("activation", &activation.cv_params, &mut metadata);
            }
        }
    }

    if let Some(product_list) = spectrum_product_list(spectrum) {
        push_valued_cv_params("product_list", &product_list.cv_params, &mut metadata);

        for product in &product_list.products {
            push_valued_cv_params("product", &product.cv_params, &mut metadata);

            if let Some(isolation_window) = product.isolation_window.as_ref() {
                push_valued_cv_params(
                    "product_isolation_window",
                    &isolation_window.cv_params,
                    &mut metadata,
                );
            }
        }
    }

    metadata
}

#[inline]
fn push_valued_cv_params(
    section: &'static str,
    params: &[octo::CvParam],
    out: &mut Vec<ScanMetadataEntry>,
) {
    for param in params {
        if let Some(entry) = ScanMetadataEntry::from_cv_param(section, param) {
            out.push(entry);
        }
    }
}

#[inline]
fn push_metadata_value(
    section: &'static str,
    name: &'static str,
    value: Option<&str>,
    out: &mut Vec<ScanMetadataEntry>,
) {
    let Some(value) = value.map(str::trim).filter(|v| !v.is_empty()) else {
        return;
    };
    out.push(ScanMetadataEntry {
        section,
        accession: None,
        name: name.to_owned(),
        value: value.to_owned(),
        unit_accession: None,
        unit_name: None,
    });
}
