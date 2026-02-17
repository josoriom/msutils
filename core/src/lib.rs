use core::ffi::c_int;

use serde_json::json;
use std::{
    panic::{AssertUnwindSafe, catch_unwind},
    ptr, slice,
};

use octo::{bin_to_mzml as bin_to_mzml_rs, decode, encode, parse_mzml as parse_mzml_rs};

pub mod utilities;
use utilities::{
    calculate_baseline::{BaselineOptions, calculate_baseline as calculate_baseline_rs},
    calculate_eic::TimeUnit,
    calculate_eic::{
        EicOptions, calculate_eic_from_bin1, collect_ms1_scans as collect_ms1_scans_rs,
    },
    find_feature::{FindFeatureOptions, find_feature as find_feature_rs},
    find_features::{FindFeaturesOptions, MzScanGrid, find_features as find_features_rs},
    find_noise_level::find_noise_level as find_noise_level_rs,
    find_peaks::{FilterPeaksOptions, FindPeaksOptions, find_peaks as find_peaks_rs},
    get_boundaries::BoundariesOptions,
    get_peak::get_peak as get_peak_rs,
    get_peaks_from_chrom::get_peaks_from_chrom as get_peaks_from_chrom_rs,
    get_peaks_from_eic::get_peaks_from_eic as get_peaks_from_eic_rs,
    structs::{ChromRoi, EicRoi},
    structs::{DataXY, FromTo, Roi},
};

const OK: c_int = 0;
const ERR_INVALID_ARGS: c_int = 1;
const ERR_PANIC: c_int = 2;
const ERR_PARSE: c_int = 4;
const ERR_ENCODE: c_int = 5;

#[repr(C)]
pub struct Buf {
    pub ptr: *mut u8,
    pub len: usize,
}

#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct CPeakPOptions {
    pub integral_threshold: f64,
    pub intensity_threshold: f64,
    pub width_threshold: c_int,
    pub noise: f64,
    pub auto_noise: c_int,
    pub auto_baseline: c_int,
    pub baseline_window: c_int,
    pub baseline_window_factor: c_int,
    pub allow_overlap: c_int,
    pub window_size: c_int,
    pub sn_ratio: f64,
}

#[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
#[link(wasm_import_module = "env")]
unsafe extern "C" {
    fn js_log(ptr: *const u8, len: usize);
}

#[inline]
pub fn log_json<T: serde::Serialize>(v: &T) {
    if let Ok(s) = serde_json::to_string_pretty(v) {
        #[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
        unsafe {
            js_log(s.as_ptr(), s.len());
        }

        #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
        eprintln!("{s}");
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn alloc(size: usize) -> *mut u8 {
    if size == 0 {
        return core::ptr::null_mut();
    }
    let mut v = Vec::<u8>::with_capacity(size);
    let p = v.as_mut_ptr();
    core::mem::forget(v);
    p
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn free_(ptr_raw: *mut u8, len: usize) {
    if !ptr_raw.is_null() {
        let slice = unsafe { core::slice::from_raw_parts_mut(ptr_raw, len) };
        drop(unsafe { Box::<[u8]>::from_raw(slice) });
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn parse_mzml(
    data_ptr: *const u8,
    data_len: usize,
    out_data: *mut Buf,
) -> c_int {
    if data_ptr.is_null() || out_data.is_null() {
        return ERR_INVALID_ARGS;
    }
    let res = catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let data = unsafe { slice::from_raw_parts(data_ptr, data_len) };
        let parsed = parse_mzml_rs(data, false).map_err(|_| ERR_PARSE)?;
        let bin: Vec<u8> = encode(&parsed, 0, false);
        write_buf(out_data, bin.into_boxed_slice());
        Ok(())
    }));
    match res {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => ERR_PANIC,
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn parse_bin(
    bin_ptr: *const u8,
    bin_len: usize,
    out_blob: *mut Buf,
) -> c_int {
    if bin_ptr.is_null() || out_blob.is_null() {
        return ERR_INVALID_ARGS;
    }
    let res = catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let bin = unsafe { slice::from_raw_parts(bin_ptr, bin_len) };
        let mzml = decode(bin).map_err(|_| ERR_PARSE)?;
        let out_bytes = encode(&mzml, 0, false);
        write_buf(out_blob, out_bytes.into_boxed_slice());
        Ok(())
    }));
    match res {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => ERR_PANIC,
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn get_peak(
    x_ptr: *const f64,
    y_ptr: *const f64,
    len: usize,
    rt: f64,
    range: f64,
    options: *const CPeakPOptions,
    out_json: *mut Buf,
) -> c_int {
    if x_ptr.is_null() || y_ptr.is_null() || out_json.is_null() || len < 3 {
        return ERR_INVALID_ARGS;
    }
    let res = catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let xs = unsafe { slice::from_raw_parts(x_ptr, len) };
        let ys = unsafe { slice::from_raw_parts(y_ptr, len) };
        let data = DataXY {
            x: xs.to_vec(),
            y: ys.to_vec(),
        };

        let fp_opts = build_find_peaks_options(options);
        let roi = Roi { rt, window: range };
        let peak = get_peak_rs(&data, roi, Some(fp_opts));

        let s = match peak {
            Some(p) => serde_json::json!({
                "from": p.from,
                "to": p.to,
                "rt": p.rt,
                "integral": p.integral,
                "intensity": p.intensity,
                "ratio": p.ratio,
                "np": p.np
            })
            .to_string(),
            None => r#"{"from":0,"to":0,"rt":0,"integral":0,"intensity":0,"ratio":0,"np":0}"#
                .to_string(),
        };
        write_buf(out_json, s.into_bytes().into_boxed_slice());
        Ok(())
    }));
    match res {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => ERR_PANIC,
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn get_peaks_from_eic(
    bin_ptr: *const u8,
    bin_len: usize,
    rts_ptr: *const f64,
    mzs_ptr: *const f64,
    ranges_ptr: *const f64,
    ids_off_ptr: *const u32,
    ids_len_ptr: *const u32,
    ids_buf_ptr: *const u8,
    ids_buf_len: usize,
    n_items: usize,
    from_left: f64,
    to_right: f64,
    options: *const CPeakPOptions,
    cores: usize,
    out_json: *mut Buf,
) -> i32 {
    if bin_ptr.is_null()
        || rts_ptr.is_null()
        || mzs_ptr.is_null()
        || ranges_ptr.is_null()
        || out_json.is_null()
        || n_items == 0
    {
        return ERR_INVALID_ARGS;
    }
    let run = || -> Result<(), i32> {
        let bytes = unsafe { std::slice::from_raw_parts(bin_ptr, bin_len) };
        let rts = unsafe { std::slice::from_raw_parts(rts_ptr, n_items) };
        let mzs = unsafe { std::slice::from_raw_parts(mzs_ptr, n_items) };
        let ranges = unsafe { std::slice::from_raw_parts(ranges_ptr, n_items) };

        let has_ids = !(ids_off_ptr.is_null()
            || ids_len_ptr.is_null()
            || ids_buf_ptr.is_null()
            || ids_buf_len == 0);
        let (offs, lens, ibuf) = if has_ids {
            (
                unsafe { std::slice::from_raw_parts(ids_off_ptr, n_items) },
                unsafe { std::slice::from_raw_parts(ids_len_ptr, n_items) },
                Some(unsafe { std::slice::from_raw_parts(ids_buf_ptr, ids_buf_len) }),
            )
        } else {
            (&[][..], &[][..], None)
        };

        let mut items: Vec<EicRoi> = Vec::with_capacity(n_items);
        for i in 0..n_items {
            let rt = rts[i];
            let mz = mzs[i];
            let win = ranges[i];
            let ok = rt.is_finite() && mz.is_finite() && win.is_finite() && win > 0.0;

            let id = if let Some(buf) = ibuf {
                if has_ids {
                    let o = offs[i] as usize;
                    let l = lens[i] as usize;
                    if o.checked_add(l).map_or(true, |e| e > buf.len()) {
                        String::new()
                    } else {
                        std::str::from_utf8(&buf[o..o + l])
                            .unwrap_or("")
                            .to_string()
                    }
                } else {
                    String::new()
                }
            } else {
                String::new()
            };

            if ok {
                items.push(EicRoi {
                    id,
                    rt,
                    mz,
                    window: win,
                });
            } else {
                items.push(EicRoi {
                    id: String::new(),
                    rt: 0.0,
                    mz: 0.0,
                    window: 0.0,
                });
            }
        }

        let window = FromTo {
            from: from_left,
            to: to_right,
        };
        let fp = build_find_peaks_options(options);
        let peaks = get_peaks_from_eic_rs(bytes, window, items.as_slice(), Some(fp), cores)
            .ok_or(ERR_PARSE)?;

        let mut arr = Vec::with_capacity(peaks.len());
        for (id, ort, mz, p) in peaks {
            arr.push(serde_json::json!({
                "id": id,
                "mz": mz,
                "ort": ort,
                "rt": p.rt,
                "from": p.from,
                "to": p.to,
                "intensity": p.intensity,
                "integral": p.integral,
                "noise": p.noise
            }));
        }
        let s = serde_json::to_string(&arr).map_err(|_| ERR_PARSE)?;
        write_buf(out_json, s.into_bytes().into_boxed_slice());
        Ok(())
    };
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(run)) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => ERR_PANIC,
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn get_peaks_from_chrom(
    bin_ptr: *const u8,
    bin_len: usize,
    idxs_ptr: *const u32,
    rts_ptr: *const f64,
    ranges_ptr: *const f64,
    n_items: usize,
    options: *const CPeakPOptions,
    cores: usize,
    out_json: *mut Buf,
) -> i32 {
    if bin_ptr.is_null()
        || idxs_ptr.is_null()
        || rts_ptr.is_null()
        || ranges_ptr.is_null()
        || out_json.is_null()
        || n_items == 0
    {
        return ERR_INVALID_ARGS;
    }
    let run = || -> Result<(), i32> {
        let bin = unsafe { std::slice::from_raw_parts(bin_ptr, bin_len) };
        let idxs = unsafe { std::slice::from_raw_parts(idxs_ptr, n_items) };
        let rts = unsafe { std::slice::from_raw_parts(rts_ptr, n_items) };
        let wins = unsafe { std::slice::from_raw_parts(ranges_ptr, n_items) };
        let mzml = decode(bin).map_err(|_| ERR_PARSE)?;
        let chroms = &mzml
            .run
            .chromatogram_list
            .as_ref()
            .ok_or(ERR_PARSE)?
            .chromatograms;

        let mut items = Vec::with_capacity(n_items);
        for i in 0..n_items {
            let iu = idxs[i];
            if iu == u32::MAX {
                items.push(ChromRoi {
                    id: String::new(),
                    idx: usize::MAX,
                    rt: 0.0,
                    window: 0.0,
                });
                continue;
            }
            let idx = iu as usize;
            if idx >= chroms.len() {
                items.push(ChromRoi {
                    id: String::new(),
                    idx,
                    rt: 0.0,
                    window: 0.0,
                });
                continue;
            }
            let id = chroms[idx].id.clone();
            items.push(ChromRoi {
                id,
                idx,
                rt: rts[i],
                window: wins[i],
            });
        }

        let fp = build_find_peaks_options(options);
        let list =
            get_peaks_from_chrom_rs(&mzml, items.as_slice(), Some(fp), cores).ok_or(ERR_PARSE)?;

        let mut out = Vec::with_capacity(list.len());
        for (index, id, ort, rt, from_, to_, intensity, integral, total_area, timestamp) in list {
            out.push(serde_json::json!({
                "index": index,
                "id": id,
                "ort": ort,
                "rt": rt,
                "from": from_,
                "to": to_,
                "intensity": intensity,
                "integral":  integral,
                "total_area": total_area,
                "timestamp": timestamp,
            }));
        }
        let s = serde_json::to_string(&out).map_err(|_| ERR_PARSE)?;
        write_buf(out_json, s.into_bytes().into_boxed_slice());
        Ok(())
    };
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(run)) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => ERR_PANIC,
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn find_peaks(
    x_ptr: *const f64,
    y_ptr: *const f64,
    len: usize,
    options: *const CPeakPOptions,
    out_json: *mut Buf,
) -> c_int {
    if x_ptr.is_null() || y_ptr.is_null() || out_json.is_null() {
        return ERR_INVALID_ARGS;
    }
    let run = || -> Result<(), c_int> {
        let xs = unsafe { slice::from_raw_parts(x_ptr, len) };
        let ys = unsafe { slice::from_raw_parts(y_ptr, len) };
        if xs.len() != ys.len() || xs.len() < 3 {
            return Err(ERR_INVALID_ARGS);
        }
        let data = DataXY {
            x: xs.to_vec(),
            y: ys.to_vec(),
        };
        let opts = build_find_peaks_options(options);
        let peaks = find_peaks_rs(&data, Some(opts));
        let list: Vec<_> = peaks
            .iter()
            .map(|p| {
                json!({
                    "from": p.from,
                    "to": p.to,
                    "rt": p.rt,
                    "integral": p.integral,
                    "intensity": p.intensity,
                    "ratio": p.ratio,
                    "np": p.np,
                    "noise": p.noise
                })
            })
            .collect();
        let s = serde_json::to_string(&list).map_err(|_| ERR_PARSE)?;
        write_buf(out_json, s.into_bytes().into_boxed_slice());
        Ok(())
    };
    match catch_unwind(AssertUnwindSafe(run)) {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => ERR_PANIC,
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn find_noise_level(y_ptr: *const f32, len: usize) -> f32 {
    if y_ptr.is_null() || len == 0 {
        return f32::INFINITY;
    }
    let compute = || {
        let ys = unsafe { slice::from_raw_parts(y_ptr, len) };
        find_noise_level_rs(ys)
    };
    match catch_unwind(AssertUnwindSafe(compute)) {
        Ok(noise) => noise as f32,
        Err(_) => f32::INFINITY,
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn bin_to_json(
    bin_ptr: *const u8,
    bin_len: usize,
    out_json: *mut Buf,
) -> c_int {
    if bin_ptr.is_null() || out_json.is_null() {
        return ERR_INVALID_ARGS;
    }
    let res = catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let bin: &[u8] = unsafe { slice::from_raw_parts(bin_ptr, bin_len) };
        let mzml = decode(bin).map_err(|_| ERR_PARSE)?;
        let s = serde_json::to_string(&mzml).map_err(|_| ERR_PARSE)?;
        write_buf(out_json, s.into_bytes().into_boxed_slice());
        Ok(())
    }));

    match res {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => ERR_PANIC,
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn bin_to_mzml(
    bin_ptr: *const u8,
    bin_len: usize,
    out_mzml: *mut Buf,
) -> c_int {
    if bin_ptr.is_null() || out_mzml.is_null() {
        return ERR_INVALID_ARGS;
    }
    let res = catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let bin = unsafe { slice::from_raw_parts(bin_ptr, bin_len) };
        let mzml = decode(bin).map_err(|_| ERR_PARSE)?;
        let xml_bytes = bin_to_mzml_rs(&mzml).map_err(|_| ERR_ENCODE)?;
        write_buf(out_mzml, xml_bytes.into_bytes().into_boxed_slice());

        Ok(())
    }));
    match res {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => ERR_PANIC,
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn calculate_eic(
    bin_ptr: *const u8,
    bin_len: usize,
    target: f64,
    from_time: f64,
    to_time: f64,
    ppm_tolerance: f64,
    mz_tolerance: f64,
    out_x: *mut Buf,
    out_y: *mut Buf,
) -> c_int {
    if bin_ptr.is_null() || out_x.is_null() || out_y.is_null() {
        return ERR_INVALID_ARGS;
    }
    let res = catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let bin = unsafe { slice::from_raw_parts(bin_ptr, bin_len) };

        let eic = calculate_eic_from_bin1(
            bin,
            &target,
            FromTo {
                from: from_time,
                to: to_time,
            },
            EicOptions {
                ppm_tolerance,
                mz_tolerance,
                ..Default::default()
            },
        )
        .map_err(|_| ERR_PARSE)?;

        let x_bytes = f64_slice_to_u8_box(&eic.x);
        let y_bytes = f64_slice_to_u8_box(&eic.y);
        write_buf(out_x, x_bytes);
        write_buf(out_y, y_bytes);
        Ok(())
    }));
    match res {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => ERR_PANIC,
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn find_features(
    data_ptr: *const u8,
    data_len: usize,
    from_time: f64,
    to_time: f64,
    eic_ppm_tolerance: f64,
    eic_mz_tolerance: f64,
    grid_start: f64,
    grid_end: f64,
    grid_step: f64,
    peak_opts: *const CPeakPOptions,
    cores: c_int,
    out_json: *mut Buf,
) -> c_int {
    if data_ptr.is_null() || out_json.is_null() || !from_time.is_finite() || !to_time.is_finite() {
        return ERR_INVALID_ARGS;
    }
    if !(to_time > from_time) {
        return ERR_INVALID_ARGS;
    }

    let run = || -> Result<(), c_int> {
        let bytes = unsafe { slice::from_raw_parts(data_ptr, data_len) };

        let mzml = decode(bytes).map_err(|_| ERR_PARSE)?;

        let mut eic_opts = EicOptions::default();
        if eic_ppm_tolerance.is_finite() && eic_ppm_tolerance >= 0.0 {
            eic_opts.ppm_tolerance = eic_ppm_tolerance;
        }
        if eic_mz_tolerance.is_finite() && eic_mz_tolerance >= 0.0 {
            eic_opts.mz_tolerance = eic_mz_tolerance;
        }

        let mut mzr = MzScanGrid::default();
        if grid_start.is_finite() {
            mzr.mz_min = grid_start;
        }
        if grid_end.is_finite() {
            mzr.mz_max = grid_end;
        }
        if grid_step > 0.0 {
            mzr.step_size = grid_step as f64;
        }

        let fp_opts = build_find_peaks_options(peak_opts);

        let feats = find_features_rs(
            &mzml,
            FromTo {
                from: from_time,
                to: to_time,
            },
            Some(FindFeaturesOptions {
                eic_options: Some(eic_opts),
                find_peaks: Some(fp_opts),
                mz_scan_grid: Some(mzr),
                ..Default::default()
            }),
            if cores > 0 { cores as usize } else { 1 },
        );

        let mut arr = Vec::with_capacity(feats.len());
        for f in feats {
            arr.push(serde_json::json!({
                "mz":        f64_ok(f.mz),
                "rt":        f64_ok(f.rt),
                "intensity": f64_ok(f.intensity),
                "from": f64_ok(f.from),
                "to": f64_ok(f.to),
            }));
        }
        let s = serde_json::to_string(&arr).map_err(|_| ERR_PARSE)?;
        write_buf(out_json, s.into_bytes().into_boxed_slice());
        Ok(())
    };

    match catch_unwind(AssertUnwindSafe(run)) {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => ERR_PANIC,
    }
}

fn f64_ok(v: f64) -> f64 {
    if v.is_finite() { v } else { 0.0 }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn calculate_baseline(
    y_ptr: *const f64,
    len: usize,
    baseline_window: c_int,
    baseline_window_factor: c_int,
    out_baseline: *mut Buf,
) -> c_int {
    if y_ptr.is_null() || out_baseline.is_null() {
        return ERR_INVALID_ARGS;
    }
    let res = catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let ys = unsafe { slice::from_raw_parts(y_ptr, len) };
        let def = BaselineOptions::default();
        let opts = BaselineOptions {
            baseline_window: if baseline_window > 0 {
                Some(baseline_window as f64)
            } else {
                def.baseline_window
            },
            baseline_window_factor: if baseline_window_factor > 0 {
                Some(baseline_window_factor as usize)
            } else {
                def.baseline_window_factor
            },
            level: Some(1),
        };
        let base = calculate_baseline_rs(ys, opts);
        let bytes = f64_slice_to_u8_box(&base);
        write_buf(out_baseline, bytes);
        Ok(())
    }));
    match res {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => ERR_PANIC,
    }
}

fn f64_slice_to_u8_box(v: &[f64]) -> Box<[u8]> {
    let n = v.len() * 8;
    let mut out = Vec::<u8>::with_capacity(n);
    unsafe {
        out.set_len(n);
        ptr::copy_nonoverlapping(v.as_ptr() as *const u8, out.as_mut_ptr(), n);
    }
    out.into_boxed_slice()
}

fn write_buf(out: *mut Buf, bytes: Box<[u8]>) {
    let len = bytes.len();
    let ptr_bytes = Box::into_raw(bytes) as *mut u8;
    unsafe {
        ptr::write_unaligned(
            out,
            Buf {
                ptr: ptr_bytes,
                len,
            },
        )
    };
}

fn build_find_peaks_options(options: *const CPeakPOptions) -> FindPeaksOptions {
    if options.is_null() {
        return FindPeaksOptions {
            get_boundaries_options: Some(BoundariesOptions {
                ..Default::default()
            }),
            filter_peaks_options: Some(FilterPeaksOptions {
                ..Default::default()
            }),
            baseline_options: Some(BaselineOptions {
                ..Default::default()
            }),
        };
    }
    let options = unsafe { *options };
    let integral = (options.integral_threshold.is_finite() && options.integral_threshold >= 0.0)
        .then_some(options.integral_threshold);
    let intensity = (options.intensity_threshold.is_finite() && options.intensity_threshold >= 0.0)
        .then_some(options.intensity_threshold);
    let width = (options.width_threshold > 0).then_some(options.width_threshold as usize);
    let noise = (options.noise.is_finite() && options.noise > 0.0).then_some(options.noise);
    let auto_noise = Some(options.auto_noise != 0);
    let auto_baseline = Some(options.auto_baseline != 0);
    let baseline_window = (options.baseline_window > 0).then_some(options.baseline_window as f64);
    let baseline_window_factor =
        (options.baseline_window_factor > 0).then_some(options.baseline_window_factor as usize);
    let allow_overlap = Some(options.allow_overlap != 0);
    let sn_ratio = if options.sn_ratio.is_finite() && options.sn_ratio > 0.0 {
        Some(options.sn_ratio)
    } else {
        Some(1.5)
    };
    let filter = FilterPeaksOptions {
        integral_threshold: integral,
        intensity_threshold: intensity,
        width_threshold: width,
        noise,
        auto_noise,
        auto_baseline,
        allow_overlap,
        sn_ratio,
        ..Default::default()
    };
    FindPeaksOptions {
        get_boundaries_options: Some(BoundariesOptions {
            ..Default::default()
        }),
        filter_peaks_options: Some(filter),
        baseline_options: Some(BaselineOptions {
            baseline_window,
            baseline_window_factor,
            level: Some(1),
        }),
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn collect_ms1_scans(
    bin_ptr: *const u8,
    bin_len: usize,
    from_time: f64,
    to_time: f64,
    out_rt: *mut Buf,
    out_offsets: *mut Buf,
    out_lengths: *mut Buf,
    out_mz: *mut Buf,
    out_intensity: *mut Buf,
) -> c_int {
    if bin_ptr.is_null()
        || out_rt.is_null()
        || out_offsets.is_null()
        || out_lengths.is_null()
        || out_mz.is_null()
        || out_intensity.is_null()
        || !from_time.is_finite()
        || !to_time.is_finite()
        || !(to_time >= from_time)
    {
        return ERR_INVALID_ARGS;
    }

    let run = || -> Result<(), c_int> {
        let bin = unsafe { slice::from_raw_parts(bin_ptr, bin_len) };
        let mzml = decode(bin).map_err(|_| ERR_PARSE)?;
        let (rt, scans) = collect_ms1_scans_rs(
            &mzml,
            FromTo {
                from: from_time,
                to: to_time,
            },
            TimeUnit::Minutes,
        );

        let n_scans = scans.len();
        let mut offsets = Vec::<u32>::with_capacity(n_scans);
        let mut lengths = Vec::<u32>::with_capacity(n_scans);

        let total_points: usize = scans
            .iter()
            .map(|s| s.mz.len().min(s.intensity.len()))
            .sum();

        if total_points > (u32::MAX as usize) {
            return Err(ERR_PARSE);
        }

        let mut mz_all = Vec::<f64>::with_capacity(total_points);
        let mut intensity_all = Vec::<f64>::with_capacity(total_points);

        let mut cursor: usize = 0;
        for s in &scans {
            let len = s.mz.len().min(s.intensity.len());
            offsets.push(cursor as u32);
            lengths.push(len as u32);

            mz_all.extend_from_slice(&s.mz[..len]);
            intensity_all.extend_from_slice(&s.intensity[..len]);

            cursor += len;
        }

        let rt_bytes = f64_slice_to_u8_box(&rt);
        let off_bytes = u32_slice_to_u8_box(&offsets);
        let len_bytes = u32_slice_to_u8_box(&lengths);
        let mz_bytes = f64_slice_to_u8_box(&mz_all);
        let it_bytes = f64_slice_to_u8_box(&intensity_all);

        write_buf(out_rt, rt_bytes);
        write_buf(out_offsets, off_bytes);
        write_buf(out_lengths, len_bytes);
        write_buf(out_mz, mz_bytes);
        write_buf(out_intensity, it_bytes);

        Ok(())
    };

    match catch_unwind(AssertUnwindSafe(run)) {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => ERR_PANIC,
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn find_feature(
    bin_ptr: *const u8,
    bin_len: usize,
    rts_ptr: *const f64,
    mzs_ptr: *const f64,
    windows_ptr: *const f64,
    ids_off_ptr: *const u32,
    ids_len_ptr: *const u32,
    ids_buf_ptr: *const u8,
    ids_buf_len: usize,
    n_items: usize,
    cores: usize,
    scan_eic_ppm_tolerance: f64,
    scan_eic_mz_tolerance: f64,
    eic_ppm_tolerance: f64,
    eic_mz_tolerance: f64,
    peak_opts: *const CPeakPOptions,
    out_json: *mut Buf,
) -> c_int {
    if bin_ptr.is_null()
        || rts_ptr.is_null()
        || mzs_ptr.is_null()
        || windows_ptr.is_null()
        || out_json.is_null()
        || n_items == 0
    {
        return ERR_INVALID_ARGS;
    }

    let run = || -> Result<(), c_int> {
        let bin = unsafe { slice::from_raw_parts(bin_ptr, bin_len) };
        let rts = unsafe { slice::from_raw_parts(rts_ptr, n_items) };
        let mzs = unsafe { slice::from_raw_parts(mzs_ptr, n_items) };
        let wins = unsafe { slice::from_raw_parts(windows_ptr, n_items) };

        let has_ids = !(ids_off_ptr.is_null()
            || ids_len_ptr.is_null()
            || ids_buf_ptr.is_null()
            || ids_buf_len == 0);
        let (offs, lens, ibuf) = if has_ids {
            (
                unsafe { slice::from_raw_parts(ids_off_ptr, n_items) },
                unsafe { slice::from_raw_parts(ids_len_ptr, n_items) },
                Some(unsafe { slice::from_raw_parts(ids_buf_ptr, ids_buf_len) }),
            )
        } else {
            (&[][..], &[][..], None)
        };

        let mzml = decode(bin).map_err(|_| ERR_PARSE)?;

        let mut scan_opts = EicOptions {
            ppm_tolerance: 10.0,
            mz_tolerance: 0.003,
            ..Default::default()
        };
        if scan_eic_ppm_tolerance.is_finite() && scan_eic_ppm_tolerance >= 0.0 {
            scan_opts.ppm_tolerance = scan_eic_ppm_tolerance;
        }
        if scan_eic_mz_tolerance.is_finite() && scan_eic_mz_tolerance >= 0.0 {
            scan_opts.mz_tolerance = scan_eic_mz_tolerance;
        }

        let mut eic_opts = EicOptions {
            ppm_tolerance: 20.0,
            mz_tolerance: 0.005,
            ..Default::default()
        };
        if eic_ppm_tolerance.is_finite() && eic_ppm_tolerance >= 0.0 {
            eic_opts.ppm_tolerance = eic_ppm_tolerance;
        }
        if eic_mz_tolerance.is_finite() && eic_mz_tolerance >= 0.0 {
            eic_opts.mz_tolerance = eic_mz_tolerance;
        }

        let fp_opts = build_find_peaks_options(peak_opts);
        let opts = FindFeatureOptions {
            scan_eic_options: Some(scan_opts),
            eic_options: Some(eic_opts),
            find_peaks: Some(fp_opts),
        };

        let mut rois_owned: Vec<EicRoi> = Vec::with_capacity(n_items);
        for i in 0..n_items {
            let rt = rts[i];
            let mz = mzs[i];
            let window = wins[i];
            let ok = rt.is_finite() && mz.is_finite() && window.is_finite() && window > 0.0;

            let id = if let Some(buf) = ibuf {
                if has_ids {
                    let o = offs[i] as usize;
                    let l = lens[i] as usize;
                    if o.checked_add(l).map_or(true, |e| e > buf.len()) {
                        String::new()
                    } else {
                        std::str::from_utf8(&buf[o..o + l])
                            .unwrap_or("")
                            .to_string()
                    }
                } else {
                    String::new()
                }
            } else {
                String::new()
            };

            if ok {
                rois_owned.push(EicRoi { id, rt, mz, window });
            } else {
                rois_owned.push(EicRoi {
                    id,
                    rt: 0.0,
                    mz: 0.0,
                    window: 0.0,
                });
            }
        }

        let mut roi_refs: Vec<&EicRoi> = Vec::with_capacity(rois_owned.len());
        for r in &rois_owned {
            roi_refs.push(r);
        }

        let results = find_feature_rs(&mzml, roi_refs.as_slice(), cores, Some(opts));

        let mut arr = Vec::with_capacity(results.len());
        for (i, out) in results.iter().enumerate() {
            let jid = &rois_owned[i].id;
            match out {
                Some(f) => {
                    arr.push(serde_json::json!({
                        "id": f.id,
                        "mz": f.mz,
                        "rt": f.rt,
                        "from": f.peak.from,
                        "to":   f.peak.to,
                        "intensity": f.peak.intensity,
                        "integral":  f.peak.integral,
                        "ratio":     f.peak.ratio,
                        "np":        f.peak.np,
                        "noise":     f.peak.noise
                    }));
                }
                None => {
                    arr.push(serde_json::json!({
                        "id": jid,
                        "mz": 0.0,
                        "rt": 0.0,
                        "from": 0.0,
                        "to": 0.0,
                        "intensity": 0.0,
                        "integral": 0.0,
                        "ratio": 0.0,
                        "np": 0,
                        "noise": 0.0
                    }));
                }
            }
        }

        let s = serde_json::to_string(&arr).map_err(|_| ERR_ENCODE)?;
        write_buf(out_json, s.into_bytes().into_boxed_slice());
        Ok(())
    };

    match catch_unwind(AssertUnwindSafe(run)) {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => ERR_PANIC,
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn convert_mzml_to_bin(
    data_ptr: *const u8,
    data_len: usize,
    out_blob: *mut Buf,
    level: u8,
    f32_compress: u8,
) -> c_int {
    if data_ptr.is_null() || out_blob.is_null() {
        return ERR_INVALID_ARGS;
    }
    let res = catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let data = unsafe { slice::from_raw_parts(data_ptr, data_len) };
        let mzml = parse_mzml_rs(data, false).map_err(|_| ERR_PARSE)?;
        let bin = encode(&mzml, level, f32_compress != 0);
        write_buf(out_blob, bin.into_boxed_slice());
        Ok(())
    }));
    match res {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => ERR_PANIC,
    }
}

fn u32_slice_to_u8_box(v: &[u32]) -> Box<[u8]> {
    let n = v.len() * 4;
    let mut out = Vec::<u8>::with_capacity(n);
    unsafe {
        out.set_len(n);
        ptr::copy_nonoverlapping(v.as_ptr() as *const u8, out.as_mut_ptr(), n);
    }
    out.into_boxed_slice()
}
