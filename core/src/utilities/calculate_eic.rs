use std::{cmp::Ordering, sync::Arc};

use octo::{BinaryData, BinaryDataArray, BinaryDataArrayList, MzML, Spectrum, decode};

use crate::utilities::structs::{FromTo, Peak};

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

pub fn calculate_eic_from_bin1(
    bin1: &[u8],
    target_mass: &f64,
    from_to: FromTo,
    options: EicOptions,
) -> Result<Eic, &'static str> {
    let mzml = decode(bin1).map_err(|_| "decode BIN1 failed")?;
    calculate_eic_from_mzml(&mzml, target_mass, from_to, options)
}

pub fn calculate_eic_from_mzml(
    mzml: &MzML,
    target_mass: &f64,
    from_to: FromTo,
    options: EicOptions,
) -> Result<Eic, &'static str> {
    let (times, scans) = collect_ms1_scans(mzml, from_to, options.time_unit);

    if scans.is_empty() || times.is_empty() {
        return Ok(Eic {
            x: Vec::new(),
            y: Vec::new(),
        });
    }

    let y = compute_eic_for_mz(&scans, times.len(), *target_mass, options);
    Ok(Eic { x: times, y })
}

#[derive(Clone)]
pub struct CentroidScan {
    pub rt: f64,
    pub mz: Arc<[f64]>,
    pub intensity: Arc<[f64]>,
}

const ACC_MZ_ARRAY: &str = "MS:1000514";
const ACC_INTENSITY_ARRAY: &str = "MS:1000515";
const ACC_SCAN_START_TIME: &str = "MS:1000016";
const ACC_MS_LEVEL: &str = "MS:1000511";

const UO_MIN: &str = "UO:0000031";
const UO_SEC: &str = "UO:0000010";
const UO_MS: &str = "UO:0000028";

#[inline]
fn ms_level_u8(s: &Spectrum) -> Option<u8> {
    if let Some(v) = s.ms_level {
        return u8::try_from(v).ok();
    }
    for cv in &s.cv_params {
        if cv.accession.as_deref() == Some(ACC_MS_LEVEL) {
            return cv.value.as_deref().and_then(|v| v.parse::<u8>().ok());
        }
    }
    None
}

#[inline]
fn convert_rt_to_minutes(raw: f64, unit_acc: Option<&str>, unit_name: Option<&str>) -> Option<f64> {
    if !raw.is_finite() {
        return None;
    }

    let seconds = match unit_acc {
        Some(UO_SEC) => Some(raw),
        Some(UO_MIN) => Some(raw * 60.0),
        Some(UO_MS) => Some(raw / 1000.0),
        _ => match unit_name {
            Some("second") | Some("seconds") => Some(raw),
            Some("minute") | Some("minutes") => Some(raw * 60.0),
            Some("millisecond") | Some("milliseconds") => Some(raw / 1000.0),
            _ => None,
        },
    }?;

    Some(seconds / 60.0)
}

#[inline]
fn scan_start_time_minutes_from_scan_list(sl: &octo::ScanList) -> Option<f64> {
    for sc in &sl.scans {
        for cv in &sc.cv_params {
            if cv.accession.as_deref() == Some(ACC_SCAN_START_TIME) {
                let raw = cv.value.as_deref()?.parse::<f64>().ok()?;
                return convert_rt_to_minutes(
                    raw,
                    cv.unit_accession.as_deref(),
                    cv.unit_name.as_deref(),
                );
            }
        }
    }

    for cv in &sl.cv_params {
        if cv.accession.as_deref() == Some(ACC_SCAN_START_TIME) {
            let raw = cv.value.as_deref()?.parse::<f64>().ok()?;
            return convert_rt_to_minutes(
                raw,
                cv.unit_accession.as_deref(),
                cv.unit_name.as_deref(),
            );
        }
    }

    None
}

#[inline]
fn spectrum_rt_minutes(s: &Spectrum) -> Option<f64> {
    if let Some(sl) = s.scan_list.as_ref() {
        if let Some(rt) = scan_start_time_minutes_from_scan_list(sl) {
            return Some(rt);
        }
    }

    if let Some(sd) = s.spectrum_description.as_ref() {
        if let Some(sl) = sd.scan_list.as_ref() {
            if let Some(rt) = scan_start_time_minutes_from_scan_list(sl) {
                return Some(rt);
            }
        }
    }

    None
}

pub fn compute_eic_for_mz(
    scans: &[CentroidScan],
    rt_len: usize,
    center: f64,
    opts: EicOptions,
) -> Vec<f64> {
    if !center.is_finite() || center <= 0.0 {
        return vec![0.0; rt_len];
    }

    let tol_ppm = if opts.ppm_tolerance > 0.0 {
        (opts.ppm_tolerance * 1e-6) * center.abs()
    } else {
        0.0
    };
    let tol_abs = opts.mz_tolerance.max(0.0);
    let tol = tol_ppm.max(tol_abs);

    if !tol.is_finite() || tol <= 0.0 {
        return vec![0.0; rt_len];
    }

    let lo = center - tol;
    let hi = center + tol;

    let broad = (tol * 5.0).max(0.01);
    let blo = center - broad;
    let bhi = center + broad;

    let mut y = vec![0.0f64; rt_len];
    let n = rt_len.min(scans.len());

    for (i, s) in scans.iter().take(n).enumerate() {
        let mzs = s.mz.as_ref();
        let ints = s.intensity.as_ref();

        if mzs.is_empty() || ints.is_empty() || mzs.len() != ints.len() {
            continue;
        }

        let j0 = lower_bound(mzs, lo);
        let j1 = upper_bound(mzs, hi);

        if j1 > j0 {
            let mut acc = 0.0f64;
            for v in &ints[j0..j1] {
                acc += *v;
            }
            y[i] = acc;
            continue;
        }

        let k0 = lower_bound(mzs, blo);
        let k1 = upper_bound(mzs, bhi);

        if k1 > k0 {
            let mut best = 0.0f64;
            for v in &ints[k0..k1] {
                if *v > best {
                    best = *v;
                }
            }
            y[i] = best;
        }
    }

    y
}

pub fn collect_ms1_scans(
    mzml: &MzML,
    time_window: FromTo,
    time_unit: TimeUnit,
) -> (Vec<f64>, Vec<CentroidScan>) {
    let spectra: &[Spectrum] = mzml
        .run
        .spectrum_list
        .as_ref()
        .map(|sl| sl.spectra.as_slice())
        .unwrap_or(&[]);

    let win_from_raw = time_window.from.min(time_window.to);
    let win_to_raw = time_window.from.max(time_window.to);

    let win_from = time_unit.to_minutes(win_from_raw);
    let win_to = time_unit.to_minutes(win_to_raw);

    let mut scans: Vec<CentroidScan> = Vec::with_capacity(spectra.len() / 2);

    for s in spectra {
        if !matches!(ms_level_u8(s), Some(1)) {
            continue;
        }

        let rt = match spectrum_rt_minutes(s) {
            Some(v) => v,
            None => continue,
        };

        if rt < win_from || rt > win_to {
            continue;
        }

        let (mz_bda, in_bda) = match spectrum_xy(s) {
            Some(v) => v,
            None => continue,
        };

        let mz_bin = match mz_bda.binary.as_ref() {
            Some(v) => v,
            None => continue,
        };
        let in_bin = match in_bda.binary.as_ref() {
            Some(v) => v,
            None => continue,
        };

        let n = bin_len(mz_bin).min(bin_len(in_bin));
        if n == 0 {
            continue;
        }

        let mut mzs: Vec<f64> = Vec::with_capacity(n);
        let mut ints: Vec<f64> = Vec::with_capacity(n);

        let mut last_m = f64::NEG_INFINITY;
        let mut sorted = true;

        match (mz_bin, in_bin) {
            (BinaryData::F64(mzv), BinaryData::F64(itv)) => {
                for i in 0..n {
                    let m = mzv[i];
                    let it = itv[i];
                    if m.is_finite() && it.is_finite() {
                        if m < last_m {
                            sorted = false;
                        }
                        last_m = m;
                        mzs.push(m);
                        ints.push(it);
                    }
                }
            }
            (BinaryData::F64(mzv), BinaryData::F32(itv)) => {
                for i in 0..n {
                    let m = mzv[i];
                    let it = itv[i] as f64;
                    if m.is_finite() && it.is_finite() {
                        if m < last_m {
                            sorted = false;
                        }
                        last_m = m;
                        mzs.push(m);
                        ints.push(it);
                    }
                }
            }
            (BinaryData::F32(mzv), BinaryData::F64(itv)) => {
                for i in 0..n {
                    let m = mzv[i] as f64;
                    let it = itv[i];
                    if m.is_finite() && it.is_finite() {
                        if m < last_m {
                            sorted = false;
                        }
                        last_m = m;
                        mzs.push(m);
                        ints.push(it);
                    }
                }
            }
            (BinaryData::F32(mzv), BinaryData::F32(itv)) => {
                for i in 0..n {
                    let m = mzv[i] as f64;
                    let it = itv[i] as f64;
                    if m.is_finite() && it.is_finite() {
                        if m < last_m {
                            sorted = false;
                        }
                        last_m = m;
                        mzs.push(m);
                        ints.push(it);
                    }
                }
            }
            _ => {
                for i in 0..n {
                    let m = bin_get_f64(mz_bin, i);
                    let it = bin_get_f64(in_bin, i);
                    if m.is_finite() && it.is_finite() {
                        if m < last_m {
                            sorted = false;
                        }
                        last_m = m;
                        mzs.push(m);
                        ints.push(it);
                    }
                }
            }
        }

        if mzs.is_empty() {
            continue;
        }

        if !sorted {
            let mut pairs: Vec<(f64, f64)> = mzs.into_iter().zip(ints.into_iter()).collect();
            pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal));

            let mut mz_sorted: Vec<f64> = Vec::with_capacity(pairs.len());
            let mut in_sorted: Vec<f64> = Vec::with_capacity(pairs.len());

            for (m, it) in pairs {
                mz_sorted.push(m);
                in_sorted.push(it);
            }

            scans.push(CentroidScan {
                rt,
                mz: Arc::from(mz_sorted.into_boxed_slice()),
                intensity: Arc::from(in_sorted.into_boxed_slice()),
            });
        } else {
            scans.push(CentroidScan {
                rt,
                mz: Arc::from(mzs.into_boxed_slice()),
                intensity: Arc::from(ints.into_boxed_slice()),
            });
        }
    }

    scans.sort_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap_or(Ordering::Equal));

    let mut rt = Vec::with_capacity(scans.len());
    for s in &scans {
        rt.push(s.rt);
    }

    (rt, scans)
}

#[inline]
fn spectrum_xy(s: &Spectrum) -> Option<(&BinaryDataArray, &BinaryDataArray)> {
    let bal = s.binary_data_array_list.as_ref()?;
    let mz_bda = find_bda_by_accession(bal, ACC_MZ_ARRAY)?;
    let in_bda = find_bda_by_accession(bal, ACC_INTENSITY_ARRAY)?;
    Some((mz_bda, in_bda))
}

#[inline]
fn find_bda_by_accession<'a>(
    bal: &'a BinaryDataArrayList,
    accession: &str,
) -> Option<&'a BinaryDataArray> {
    bal.binary_data_arrays.iter().find(|bda| {
        bda.cv_params
            .iter()
            .any(|cv| cv.accession.as_deref() == Some(accession))
    })
}

#[inline]
fn bin_len(bin: &BinaryData) -> usize {
    match bin {
        BinaryData::F64(v) => v.len(),
        BinaryData::F32(v) => v.len(),
        BinaryData::F16(v) => v.len(),
        BinaryData::I64(v) => v.len(),
        BinaryData::I32(v) => v.len(),
        BinaryData::I16(v) => v.len(),
    }
}

#[inline]
fn bin_get_f64(bin: &BinaryData, i: usize) -> f64 {
    match bin {
        BinaryData::F64(v) => v[i],
        BinaryData::F32(v) => v[i] as f64,
        BinaryData::F16(v) => v[i] as f64,
        BinaryData::I64(v) => v[i] as f64,
        BinaryData::I32(v) => v[i] as f64,
        BinaryData::I16(v) => v[i] as f64,
    }
}

fn max_in_range(rt: &[f64], y: &[f64], from_rt: f64, to_rt: f64) -> f64 {
    let i0 = lower_bound(rt, from_rt);
    let mut i1 = upper_bound(rt, to_rt);
    if i0 >= y.len() {
        return 0.0;
    }
    if i1 > y.len() {
        i1 = y.len();
    }
    if i1 <= i0 {
        return 0.0;
    }
    let mut m = y[i0];
    let mut i = i0 + 1;
    while i < i1 {
        let v = y[i];
        if v > m {
            m = v;
        }
        i += 1;
    }
    m
}

pub fn with_eic_apex_intensity(rt: &[f64], y: &[f64], mut p: Peak) -> Peak {
    let a = max_in_range(rt, y, p.from, p.to);
    if a.is_finite() && a > 0.0 {
        p.intensity = a;
    }
    p
}

#[inline]
pub fn lower_bound(a: &[f64], x: f64) -> usize {
    let mut lo = 0usize;
    let mut hi = a.len();
    while lo < hi {
        let mid = (lo + hi) / 2;
        if a[mid] < x {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    lo
}

#[inline]
pub fn upper_bound(a: &[f64], x: f64) -> usize {
    let mut lo = 0usize;
    let mut hi = a.len();
    while lo < hi {
        let mid = (lo + hi) / 2;
        if a[mid] <= x {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    lo
}
