#[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
use core::ffi::c_int;
#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use core::ffi::{CStr, c_char, c_int};

use memmap2::Mmap;
use serde_json::json;
use std::{
    fs::File,
    panic::{AssertUnwindSafe, catch_unwind},
    path::Path,
    ptr, slice,
    sync::Arc,
};

use ionic::{
    MzML, ScanMeta, SpectrumSource, bin_to_mzml as bin_to_mzml_rs,
    decoder::decode::Ion,
    encoder::{WritingMode, encode},
    ion::DecoderConfig,
    parse_mzml as parse_mzml_rs,
};

pub mod utilities;
use utilities::{
    calculate_baseline::{BaselineOptions, calculate_baseline as calculate_baseline_rs},
    calculate_eic::{
        EicOptions, calculate_eic as calculate_eic_rs, collect_scans as collect_scans_rs,
    },
    find_feature::{FindFeatureOptions, find_feature as find_feature_rs},
    find_features::{FindFeaturesOptions, find_features as find_features_rs},
    find_noise_level::find_noise_level as find_noise_level_rs,
    find_peaks::{FilterPeaksOptions, FindPeaksOptions, find_peaks as find_peaks_rs},
    get_peak::get_peak as get_peak_rs,
    get_peaks_from_chrom::get_peaks_from_chrom as get_peaks_from_chrom_rs,
    get_peaks_from_eic::get_peaks_from_eic as get_peaks_from_eic_rs,
    structs::{ChromRoi, EicRoi},
    structs::{DataXY, FromTo, Roi},
};

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use utilities::get_features::{ConsensusAlignmentConfig, get_features as get_features_rs};

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use crate::utilities::find_features::MzTolerance;

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
pub extern "C" fn alloc(size: usize) -> *mut u8 {
    if size == 0 {
        return core::ptr::null_mut();
    }
    let mut v = Vec::<u8>::with_capacity(size);
    let p = v.as_mut_ptr();
    core::mem::forget(v);
    p
}

/// Free memory returned by `alloc`.
///
/// # Safety
/// `ptr_raw` must be a live pointer from `alloc(size)`.
/// `size` must match the `alloc` call.
/// After this call, `ptr_raw` is invalid.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn free_(ptr_raw: *mut u8, size: usize) {
    if !ptr_raw.is_null() {
        unsafe {
            drop(Vec::<u8>::from_raw_parts(ptr_raw, 0, size));
        }
    }
}

struct BytesBacking(Arc<[u8]>);

impl AsRef<[u8]> for BytesBacking {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

pub struct OwnedIon {
    inner: Ion<'static>,
    _backing: Arc<dyn AsRef<[u8]> + Send + Sync>,
}

impl OwnedIon {
    pub fn from_ion_bytes(data: Arc<[u8]>, max_cache_size: usize) -> Result<Self, String> {
        let static_ref: &'static [u8] =
            unsafe { std::slice::from_raw_parts(data.as_ptr(), data.len()) };
        let inner = Ion::open_with_config(
            static_ref,
            DecoderConfig {
                max_cached_bytes: max_cache_size,
            },
        )?;
        let backing: Arc<dyn AsRef<[u8]> + Send + Sync> = Arc::new(BytesBacking(data));
        Ok(Self {
            _backing: backing,
            inner,
        })
    }

    pub fn from_ion_path(path: &Path, max_cache_size: usize) -> Result<Self, String> {
        let file = File::open(path).map_err(|e| format!("open {}: {}", path.display(), e))?;
        let mmap =
            unsafe { Mmap::map(&file) }.map_err(|e| format!("mmap {}: {}", path.display(), e))?;
        let backing: Arc<dyn AsRef<[u8]> + Send + Sync> = Arc::new(mmap);
        let bytes: &[u8] = backing.as_ref().as_ref();
        let static_ref: &'static [u8] =
            unsafe { std::slice::from_raw_parts(bytes.as_ptr(), bytes.len()) };
        let inner = Ion::open_with_config(
            static_ref,
            DecoderConfig {
                max_cached_bytes: max_cache_size,
            },
        )?;
        Ok(Self {
            _backing: backing,
            inner,
        })
    }

    pub fn from_mzml(mzml: &MzML) -> Result<Self, String> {
        let mut buf = Vec::new();
        encode(mzml, 0, false, WritingMode::Memory, &mut buf)
            .map_err(|e| format!("encode: {e}"))?;
        Self::from_ion_bytes(Arc::from(buf), 0)
    }

    fn ensure_mzml(&mut self) -> Result<MzML, c_int> {
        self.inner.to_mzml().map_err(|_| ERR_PARSE)
    }
}

impl SpectrumSource for OwnedIon {
    fn for_each_scan_in_range(
        &mut self,
        rt_min: f64,
        rt_max: f64,
        ms_level: u8,
        cb: &mut dyn FnMut(f64, &ScanMeta, &[f64], &[f64]),
    ) {
        self.inner
            .for_each_scan_in_range(rt_min, rt_max, ms_level, cb);
    }
}

pub enum ParsedFile {
    Full(Box<MzML>),
    Lazy(Box<OwnedIon>),
}

impl ParsedFile {
    fn with_mzml<T>(&mut self, f: impl FnOnce(&MzML) -> Result<T, c_int>) -> Result<T, c_int> {
        match self {
            ParsedFile::Full(mzml) => f(mzml.as_ref()),
            ParsedFile::Lazy(file) => {
                let mzml = file.ensure_mzml()?;
                f(&mzml)
            }
        }
    }
}

impl SpectrumSource for ParsedFile {
    fn for_each_scan_in_range(
        &mut self,
        rt_min: f64,
        rt_max: f64,
        ms_level: u8,
        cb: &mut dyn FnMut(f64, &ScanMeta, &[f64], &[f64]),
    ) {
        match self {
            ParsedFile::Full(mzml) => mzml.for_each_scan_in_range(rt_min, rt_max, ms_level, cb),
            ParsedFile::Lazy(file) => file.for_each_scan_in_range(rt_min, rt_max, ms_level, cb),
        }
    }
}

/// Free a `ParsedFile` handle.
///
/// # Safety
/// `handle` must be a valid pointer returned by this library.
/// Do not free the same handle twice.
/// After this call, `handle` is invalid.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn free_mzml(handle: *mut ParsedFile) {
    if !handle.is_null() {
        drop(unsafe { Box::from_raw(handle) });
    }
}

/// Parse mzML bytes and store the result in `dest`.
///
/// # Safety
/// `data_ptr` must point to `data_len` readable bytes.
/// `dest` must be a valid writable pointer to `*mut ParsedFile`.
/// On success, `*dest` must be freed with `free_mzml`.
#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
#[unsafe(no_mangle)]
pub unsafe extern "C" fn parse_mzml(
    data_ptr: *const u8,
    data_len: usize,
    dest: *mut *mut ParsedFile,
) -> c_int {
    if data_ptr.is_null() || dest.is_null() {
        return ERR_INVALID_ARGS;
    }
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let parsed = parse_mzml_rs(unsafe { slice::from_raw_parts(data_ptr, data_len) })
            .map_err(|_| ERR_PARSE)?;
        unsafe { *dest = Box::into_raw(Box::new(ParsedFile::Full(Box::new(parsed)))) };
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => ERR_PANIC,
    }
}

#[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
#[unsafe(no_mangle)]
pub unsafe extern "C" fn parse_mzml(
    data_ptr: *const u8,
    data_len: usize,
    dest: *mut *mut ParsedFile,
) -> c_int {
    if data_ptr.is_null() || dest.is_null() {
        return ERR_INVALID_ARGS;
    }
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let parsed = parse_mzml_rs(unsafe { slice::from_raw_parts(data_ptr, data_len) })
            .map_err(|_| ERR_PARSE)?;
        unsafe { *dest = Box::into_raw(Box::new(ParsedFile::Full(Box::new(parsed)))) };
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => ERR_PANIC,
    }
}

/// Parse binary data and store the result in `dest`.
///
/// # Safety
/// `data_ptr` must point to `data_len` readable bytes.
/// `dest` must be a valid writable pointer to `*mut ParsedFile`.
/// On success, `*dest` must be freed with `free_mzml`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn parse_bin(
    data_ptr: *const u8,
    data_len: usize,
    max_cache_size: usize,
    dest: *mut *mut ParsedFile,
) -> c_int {
    if data_ptr.is_null() || dest.is_null() {
        return ERR_INVALID_ARGS;
    }
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let arc: Arc<[u8]> = Arc::from(unsafe { std::slice::from_raw_parts(data_ptr, data_len) });
        let owned = OwnedIon::from_ion_bytes(arc, max_cache_size).map_err(|_| ERR_PARSE)?;
        unsafe { *dest = Box::into_raw(Box::new(ParsedFile::Lazy(Box::new(owned)))) };
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => ERR_PANIC,
    }
}

/// # Safety
/// `h` must be a valid unique `ParsedFile` pointer from this library.
/// `out` must be a valid writable `Buf` pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn bin_to_json(h: *mut ParsedFile, out: *mut Buf) -> c_int {
    if h.is_null() || out.is_null() {
        return ERR_INVALID_ARGS;
    }
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let file = unsafe { &mut *h };
        file.with_mzml(|mzml| {
            write_buf(
                out,
                serde_json::to_string(mzml)
                    .map_err(|_| ERR_PARSE)?
                    .into_bytes()
                    .into_boxed_slice(),
            );
            Ok(())
        })?;
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => ERR_PANIC,
    }
}

/// Convert a parsed file to mzML and write the result to `out`.
///
/// # Safety
/// `h` must be a valid unique `ParsedFile` pointer from this library.
/// `out` must be a valid writable `Buf` pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn bin_to_mzml(h: *mut ParsedFile, out: *mut Buf) -> c_int {
    if h.is_null() || out.is_null() {
        return ERR_INVALID_ARGS;
    }
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let file = unsafe { &mut *h };
        file.with_mzml(|mzml| {
            write_buf(
                out,
                bin_to_mzml_rs(mzml)
                    .map_err(|_| ERR_ENCODE)?
                    .into_bytes()
                    .into_boxed_slice(),
            );
            Ok(())
        })?;
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => ERR_PANIC,
    }
}

/// Convert a parsed file to binary and write the result to `out`.
///
/// # Safety
/// `h` must be a valid unique `ParsedFile` pointer from this library.
/// `out` must be a valid writable `Buf` pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn mzml_to_bin(
    h: *mut ParsedFile,
    out: *mut Buf,
    level: u8,
    f32_compress: u8,
) -> c_int {
    if h.is_null() || out.is_null() {
        return ERR_INVALID_ARGS;
    }
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let file = unsafe { &mut *h };
        file.with_mzml(|mzml| {
            let mut buf = Vec::new();
            encode(
                mzml,
                level,
                f32_compress != 0,
                WritingMode::Memory,
                &mut buf,
            )
            .map_err(|_| ERR_ENCODE)?;
            write_buf(out, buf.into_boxed_slice());
            Ok(())
        })?;
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => ERR_PANIC,
    }
}

/// Find one peak and write the result to `out`.
///
/// # Safety
/// `x_ptr` and `y_ptr` must point to `len` readable `f64` values.
/// `options` must be null or point to a valid `CPeakPOptions`.
/// `out` must be a valid writable `Buf` pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn get_peak(
    x_ptr: *const f64,
    y_ptr: *const f64,
    len: usize,
    rt: f64,
    range: f64,
    options: *const CPeakPOptions,
    out: *mut Buf,
) -> c_int {
    if x_ptr.is_null() || y_ptr.is_null() || out.is_null() || len < 3 {
        return ERR_INVALID_ARGS;
    }
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let data = DataXY {
            x: unsafe { slice::from_raw_parts(x_ptr, len) }.to_vec(),
            y: unsafe { slice::from_raw_parts(y_ptr, len) }.to_vec(),
        };
        let s = match get_peak_rs(&data, &Roi { rt, window: range }, Some(build_fp(options))) {
            Some(p) => json!({"from":p.from,"to":p.to,"rt":p.rt,"integral":p.integral,"intensity":p.intensity,"ratio":p.ratio,"np":p.np}).to_string(),
            None => r#"{"from":0,"to":0,"rt":0,"integral":0,"intensity":0,"ratio":0,"np":0}"#.to_string(),
        };
        write_buf(out, s.into_bytes().into_boxed_slice());
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => ERR_PANIC,
    }
}

/// Find peaks from EIC data and write the result to `out`.
///
/// # Safety
/// `h` must be a valid unique `ParsedFile` pointer from this library.
/// `rts_ptr`, `mzs_ptr`, and `ranges_ptr` must point to `n` readable `f64` values.
/// If IDs are provided, `ids_off` and `ids_len` must point to `n` readable `u32` values,
/// and `ids_buf` must point to `ids_buf_len` readable bytes.
/// `opts` must be null or point to a valid `CPeakPOptions`.
/// `out` must be a valid writable `Buf` pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn get_peaks_from_eic(
    h: *const ParsedFile,
    rts_ptr: *const f64,
    mzs_ptr: *const f64,
    ranges_ptr: *const f64,
    ids_off: *const u32,
    ids_len: *const u32,
    ids_buf: *const u8,
    ids_buf_len: usize,
    n: usize,
    from: f64,
    to: f64,
    opts: *const CPeakPOptions,
    cores: usize,
    out: *mut Buf,
) -> i32 {
    if h.is_null()
        || rts_ptr.is_null()
        || mzs_ptr.is_null()
        || ranges_ptr.is_null()
        || out.is_null()
        || n == 0
    {
        return ERR_INVALID_ARGS;
    }
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), i32> {
        let file = unsafe { &mut *(h as *mut ParsedFile) };
        let items = build_eic_rois(
            unsafe { slice::from_raw_parts(rts_ptr, n) },
            unsafe { slice::from_raw_parts(mzs_ptr, n) },
            unsafe { slice::from_raw_parts(ranges_ptr, n) },
            ids_off,
            ids_len,
            ids_buf,
            ids_buf_len,
        );
        let peaks = get_peaks_from_eic_rs(
            file,
            FromTo { from, to },
            &items,
            Some(build_fp(opts)),
            cores,
        )
        .ok_or(ERR_PARSE)?;
        let arr: Vec<_> = peaks.iter().map(|(id,ort,mz,p)| json!({"id":id,"mz":mz,"ort":ort,"rt":p.rt,"from":p.from,"to":p.to,"intensity":p.intensity,"integral":p.integral,"noise":p.noise})).collect();
        write_buf(
            out,
            serde_json::to_string(&arr)
                .map_err(|_| ERR_PARSE)?
                .into_bytes()
                .into_boxed_slice(),
        );
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => ERR_PANIC,
    }
}

/// Find peaks from chromatogram data and write the result to `out`.
///
/// # Safety
/// `h` must be a valid unique `ParsedFile` pointer from this library.
/// `idxs` must point to `n` readable `u32` values.
/// `rts` and `wins` must point to `n` readable `f64` values.
/// `opts` must be null or point to a valid `CPeakPOptions`.
/// `out` must be a valid writable `Buf` pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn get_peaks_from_chrom(
    h: *mut ParsedFile,
    idxs: *const u32,
    rts: *const f64,
    wins: *const f64,
    n: usize,
    opts: *const CPeakPOptions,
    cores: usize,
    out: *mut Buf,
) -> i32 {
    if h.is_null() || idxs.is_null() || rts.is_null() || wins.is_null() || out.is_null() || n == 0 {
        return ERR_INVALID_ARGS;
    }
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), i32> {
        let file = unsafe { &mut *h };
        file.with_mzml(|mzml| {
            let idxs = unsafe { slice::from_raw_parts(idxs, n) };
            let rts = unsafe { slice::from_raw_parts(rts, n) };
            let wins = unsafe { slice::from_raw_parts(wins, n) };
            let chroms = &mzml
                .run
                .chromatogram_list
                .as_ref()
                .ok_or(ERR_PARSE)?
                .chromatograms;
            let items: Vec<_> = (0..n)
                .map(|i| {
                    let iu = idxs[i];
                    if iu == u32::MAX {
                        return ChromRoi {
                            id: String::new(),
                            idx: usize::MAX,
                            rt: 0.0,
                            window: 0.0,
                        };
                    }
                    let idx = iu as usize;
                    if idx >= chroms.len() {
                        return ChromRoi {
                            id: String::new(),
                            idx,
                            rt: 0.0,
                            window: 0.0,
                        };
                    }
                    ChromRoi {
                        id: chroms[idx].id.clone(),
                        idx,
                        rt: rts[i],
                        window: wins[i],
                    }
                })
                .collect();
            let list =
                get_peaks_from_chrom_rs(mzml, &items, Some(build_fp(opts)), cores).ok_or(ERR_PARSE)?;
            let out_arr: Vec<_> = list
                .iter()
                .map(|row| {
                    json!({
                        "index": row.index,
                        "id": &row.id,
                        "ort": row.target_rt,
                        "rt": row.peak_rt,
                        "from": row.from_rt,
                        "to": row.to_rt,
                        "intensity": row.intensity,
                        "integral": row.area,
                        "total_area": row.total_area,
                        "timestamp": &row.timestamp,
                    })
                })
                .collect();
            write_buf(
                out,
                serde_json::to_string(&out_arr)
                    .map_err(|_| ERR_PARSE)?
                    .into_bytes()
                    .into_boxed_slice(),
            );
            Ok(())
        })?;
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => ERR_PANIC,
    }
}

/// Find peaks and write the result to `out`.
///
/// # Safety
/// `x_ptr` and `y_ptr` must point to `len` readable `f64` values.
/// `opts` must be null or point to a valid `CPeakPOptions`.
/// `out` must be a valid writable `Buf` pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn find_peaks(
    x_ptr: *const f64,
    y_ptr: *const f64,
    len: usize,
    opts: *const CPeakPOptions,
    out: *mut Buf,
) -> c_int {
    if x_ptr.is_null() || y_ptr.is_null() || out.is_null() {
        return ERR_INVALID_ARGS;
    }
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let xs = unsafe { slice::from_raw_parts(x_ptr, len) };
        let ys = unsafe { slice::from_raw_parts(y_ptr, len) };
        if xs.len() != ys.len() || xs.len() < 3 {
            return Err(ERR_INVALID_ARGS);
        }
        let peaks = find_peaks_rs(
            &DataXY {
                x: xs.to_vec(),
                y: ys.to_vec(),
            },
            Some(build_fp(opts)),
        );
        let list: Vec<_> = peaks.iter().map(|p| json!({"from":p.from,"to":p.to,"rt":p.rt,"integral":p.integral,"intensity":p.intensity,"ratio":p.ratio,"np":p.np,"noise":p.noise})).collect();
        write_buf(
            out,
            serde_json::to_string(&list)
                .map_err(|_| ERR_PARSE)?
                .into_bytes()
                .into_boxed_slice(),
        );
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => ERR_PANIC,
    }
}

/// # Safety
/// `y_ptr` must point to `len` readable `f32` values.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn find_noise_level(y_ptr: *const f32, len: usize) -> f32 {
    if y_ptr.is_null() || len == 0 {
        return f32::INFINITY;
    }
    catch_unwind(AssertUnwindSafe(|| {
        find_noise_level_rs(unsafe { slice::from_raw_parts(y_ptr, len) })
    }))
    .map(|n| n.intensity as f32)
    .unwrap_or(f32::INFINITY)
}

/// Calculate an EIC and write `x` and `y` to `ox` and `oy`.
///
/// # Safety
/// `h` must be a valid unique `ParsedFile` pointer from this library.
/// `ox` and `oy` must be valid writable `Buf` pointers.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn calculate_eic(
    h: *mut ParsedFile,
    target: f64,
    from: f64,
    to: f64,
    ppm: f64,
    mz_tol: f64,
    ox: *mut Buf,
    oy: *mut Buf,
) -> c_int {
    if h.is_null() || ox.is_null() || oy.is_null() {
        return ERR_INVALID_ARGS;
    }
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let eic = calculate_eic_rs(
            unsafe { &mut *h },
            target,
            FromTo { from, to },
            EicOptions {
                ppm_tolerance: ppm,
                mz_tolerance: mz_tol,
                ..Default::default()
            },
        );
        write_buf(ox, f64_to_u8(&eic.x));
        write_buf(oy, f64_to_u8(&eic.y));
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => ERR_PANIC,
    }
}

/// Collect scans and write the result to `out`.
///
/// # Safety
/// `h` must be a valid unique `ParsedFile` pointer from this library.
/// `out` must be a valid writable `Buf` pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn collect_scans(
    h: *mut ParsedFile,
    from: f64,
    to: f64,
    level: u8,
    out: *mut Buf,
) -> c_int {
    if h.is_null() || out.is_null() {
        return ERR_INVALID_ARGS;
    }
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let (_, scans) = collect_scans_rs(
            unsafe { &mut *h },
            FromTo { from, to },
            EicOptions::default().time_unit,
            level,
        );
        let arr: Vec<_> = scans.iter().map(|s| json!({"rt":s.rt,"mz":s.mz.as_ref(),"intensity":s.intensity.as_ref(),"metadata":s.metadata})).collect();
        write_buf(
            out,
            serde_json::to_string(&arr)
                .map_err(|_| ERR_PARSE)?
                .into_bytes()
                .into_boxed_slice(),
        );
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => ERR_PANIC,
    }
}

/// Find features from files in `dir` and write the result to `out`.
///
/// # Safety
/// `dir` must point to a valid NUL-terminated C string.
/// `popts` must be null or point to a valid `CPeakPOptions`.
/// `out` must be a valid writable `Buf` pointer.
#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
#[unsafe(no_mangle)]
pub unsafe extern "C" fn get_features(
    dir: *const c_char,
    from: f64,
    to: f64,
    eic_ppm: f64,
    eic_mz: f64,
    gs: f64,
    ge: f64,
    gstep: f64,
    gppm: f64,
    gda: f64,
    grt: f64,
    freq: c_int,
    popts: *const CPeakPOptions,
    cores: c_int,
    use_gpu: c_int,
    batch_size: c_int,
    out: *mut Buf,
) -> c_int {
    if dir.is_null() || out.is_null() {
        return ERR_INVALID_ARGS;
    }
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let path = unsafe { CStr::from_ptr(dir) }
            .to_str()
            .map_err(|_| ERR_INVALID_ARGS)?;

        let mut fc = FindFeaturesOptions::default();
        if eic_ppm.is_finite() && eic_ppm >= 0.0 {
            fc.eic_options.ppm_tolerance = eic_ppm;
        }
        if eic_mz.is_finite() && eic_mz >= 0.0 {
            fc.eic_options.mz_tolerance = eic_mz;
        }
        if gs.is_finite() {
            fc.mz_scan_grid.mz_min = gs;
        }
        if ge.is_finite() {
            fc.mz_scan_grid.mz_max = ge;
        }
        if gstep > 0.0 {
            fc.mz_scan_grid.step_size = gstep;
        }
        fc.find_peaks = build_fp(popts);
        fc.use_gpu = use_gpu != 0;
        fc.batch_size = if batch_size > 0 {
            Some(batch_size as usize)
        } else {
            None
        };

        let fp2 = build_fp(popts);
        let ac = ConsensusAlignmentConfig {
            tolerance: MzTolerance {
                mz_abs: if gda > 0.0 { gda } else { 0.003 },
                ppm: if gppm > 0.0 { gppm } else { 5.0 },
            },
            rt_tolerance: if grt > 0.0 { grt } else { 0.05 },
            frequency: if freq > 0 { freq as usize } else { 1 },
            eic_options: fc.eic_options,
            peak_options: Some(fp2),
        };

        let feats = get_features_rs(
            path,
            FromTo { from, to },
            fc,
            ac,
            if cores > 0 { cores as usize } else { 1 },
        )
        .unwrap_or_default();

        let arr: Vec<_> = feats.iter().map(|f| json!({"mz":f64ok(f.mz),"rt":f64ok(f.rt),"intensity":f64ok(f.intensity),"rintensity":f64ok(f.rintensity),"from":f64ok(f.from),"to":f64ok(f.to),"integral":f64ok(f.integral),"np":f.np,"frequency":f.frequency,"rmz":f.rmz,"rint":f.rint})).collect();
        write_buf(
            out,
            serde_json::to_string(&arr)
                .map_err(|_| ERR_PARSE)?
                .into_bytes()
                .into_boxed_slice(),
        );
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => ERR_PANIC,
    }
}

/// Find features and write the result to `out`.
///
/// # Safety
/// `h` must be a valid unique `ParsedFile` pointer from this library.
/// `popts` must be null or point to a valid `CPeakPOptions`.
/// `out` must be a valid writable `Buf` pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn find_features(
    h: *mut ParsedFile,
    from: f64,
    to: f64,
    eic_ppm: f64,
    eic_mz: f64,
    gs: f64,
    ge: f64,
    gstep: f64,
    popts: *const CPeakPOptions,
    cores: c_int,
    use_gpu: c_int,
    batch_size: c_int,
    out: *mut Buf,
) -> c_int {
    if h.is_null() || out.is_null() || !from.is_finite() || !to.is_finite() || to <= from {
        return ERR_INVALID_ARGS;
    }
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let file = unsafe { &mut *h };

        let mut opts = FindFeaturesOptions::default();
        if eic_ppm.is_finite() && eic_ppm >= 0.0 {
            opts.eic_options.ppm_tolerance = eic_ppm;
        }
        if eic_mz.is_finite() && eic_mz >= 0.0 {
            opts.eic_options.mz_tolerance = eic_mz;
        }
        if gs.is_finite() {
            opts.mz_scan_grid.mz_min = gs;
        }
        if ge.is_finite() {
            opts.mz_scan_grid.mz_max = ge;
        }
        if gstep > 0.0 {
            opts.mz_scan_grid.step_size = gstep;
        }
        opts.find_peaks = build_fp(popts);
        #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
        {
            opts.use_gpu = use_gpu != 0;
            opts.batch_size = if batch_size > 0 {
                Some(batch_size as usize)
            } else {
                None
            };
        }

        let feats = find_features_rs(
            file,
            FromTo { from, to },
            Some(opts),
            if cores > 0 { cores as usize } else { 1 },
        )
        .unwrap_or_default();

        let arr: Vec<_> = feats.iter().map(|f| json!({"mz":f64ok(f.mz),"rt":f64ok(f.rt),"intensity":f64ok(f.intensity),"from":f64ok(f.from),"to":f64ok(f.to),"integral":f64ok(f.integral),"np":f64ok(f.np as f64)})).collect();
        write_buf(
            out,
            serde_json::to_string(&arr)
                .map_err(|_| ERR_PARSE)?
                .into_bytes()
                .into_boxed_slice(),
        );
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => ERR_PANIC,
    }
}

/// Calculate a baseline and write the result to `out`.
///
/// # Safety
/// `y` must point to `len` readable `f64` values.
/// `out` must be a valid writable `Buf` pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn calculate_baseline(
    y: *const f64,
    len: usize,
    bw: c_int,
    bwf: c_int,
    out: *mut Buf,
) -> c_int {
    if y.is_null() || out.is_null() {
        return ERR_INVALID_ARGS;
    }
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let def = BaselineOptions::default();
        let base = calculate_baseline_rs(
            unsafe { slice::from_raw_parts(y, len) },
            BaselineOptions {
                baseline_window: if bw > 0 {
                    Some(bw as f64)
                } else {
                    def.baseline_window
                },
                baseline_window_factor: if bwf > 0 {
                    Some(bwf as usize)
                } else {
                    def.baseline_window_factor
                },
                level: Some(1),
            },
        );
        write_buf(out, f64_to_u8(&base));
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => ERR_PANIC,
    }
}

/// Find features for the given ROIs and write the result to `out`.
///
/// # Safety
/// `h` must be a valid unique `ParsedFile` pointer from this library.
/// `rts`, `mzs`, and `wins` must point to `n` readable `f64` values.
/// If IDs are provided, `ids_off` and `ids_len` must point to `n` readable `u32` values,
/// and `ids_buf` must point to `ids_buf_len` readable bytes.
/// `popts` must be null or point to a valid `CPeakPOptions`.
/// `out` must be a valid writable `Buf` pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn find_feature(
    h: *const ParsedFile,
    rts: *const f64,
    mzs: *const f64,
    wins: *const f64,
    ids_off: *const u32,
    ids_len: *const u32,
    ids_buf: *const u8,
    ids_buf_len: usize,
    n: usize,
    cores: usize,
    sppm: f64,
    smz: f64,
    eppm: f64,
    emz: f64,
    popts: *const CPeakPOptions,
    out: *mut Buf,
) -> c_int {
    if h.is_null() || rts.is_null() || mzs.is_null() || wins.is_null() || out.is_null() || n == 0 {
        return ERR_INVALID_ARGS;
    }
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let file = unsafe { &mut *(h as *mut ParsedFile) };
        let rts = unsafe { slice::from_raw_parts(rts, n) };
        let mzs = unsafe { slice::from_raw_parts(mzs, n) };
        let wins = unsafe { slice::from_raw_parts(wins, n) };
        let mut so = EicOptions {
            ppm_tolerance: 10.0,
            mz_tolerance: 0.003,
            ..Default::default()
        };
        if sppm.is_finite() && sppm >= 0.0 {
            so.ppm_tolerance = sppm;
        }
        if smz.is_finite() && smz >= 0.0 {
            so.mz_tolerance = smz;
        }
        let mut eo = EicOptions {
            ppm_tolerance: 20.0,
            mz_tolerance: 0.005,
            ..Default::default()
        };
        if eppm.is_finite() && eppm >= 0.0 {
            eo.ppm_tolerance = eppm;
        }
        if emz.is_finite() && emz >= 0.0 {
            eo.mz_tolerance = emz;
        }
        let rois = build_eic_rois(rts, mzs, wins, ids_off, ids_len, ids_buf, ids_buf_len);
        let refs: Vec<&EicRoi> = rois.iter().collect();
        let results = find_feature_rs(
            file,
            &refs,
            cores,
            Some(FindFeatureOptions {
                scan_eic_options: Some(so),
                eic_options: Some(eo),
                find_peaks: Some(build_fp(popts)),
            }),
        );
        let arr: Vec<_> = results.iter().enumerate().map(|(i, r)| match r {
            Some(f) => json!({"id":f.id,"mz":f.mz,"rt":f.rt,"from":f.peak.from,"to":f.peak.to,"intensity":f.peak.intensity,"integral":f.peak.integral,"ratio":f.peak.ratio,"np":f.peak.np,"noise":f.peak.noise}),
            None => json!({"id":&rois[i].id,"mz":0.0,"rt":0.0,"from":0.0,"to":0.0,"intensity":0.0,"integral":0.0,"ratio":0.0,"np":0,"noise":0.0}),
        }).collect();
        write_buf(
            out,
            serde_json::to_string(&arr)
                .map_err(|_| ERR_ENCODE)?
                .into_bytes()
                .into_boxed_slice(),
        );
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => ERR_PANIC,
    }
}

fn build_eic_rois(
    rts: &[f64],
    mzs: &[f64],
    wins: &[f64],
    ids_off: *const u32,
    ids_len: *const u32,
    ids_buf: *const u8,
    ids_buf_len: usize,
) -> Vec<EicRoi> {
    let n = rts.len();

    let has = !(ids_off.is_null() || ids_len.is_null() || ids_buf.is_null() || ids_buf_len == 0);
    let (offs, lens, ibuf) = if has {
        (
            unsafe { slice::from_raw_parts(ids_off, n) },
            unsafe { slice::from_raw_parts(ids_len, n) },
            Some(unsafe { slice::from_raw_parts(ids_buf, ids_buf_len) }),
        )
    } else {
        (&[][..], &[][..], None)
    };

    (0..n)
        .map(|i| {
            let (rt, mz, w) = (rts[i], mzs[i], wins[i]);
            let ok = rt.is_finite() && mz.is_finite() && w.is_finite() && w > 0.0;
            let id = ibuf
                .and_then(|buf| {
                    let (o, l) = (offs[i] as usize, lens[i] as usize);
                    o.checked_add(l).filter(|&e| e <= buf.len()).map(|_| {
                        std::str::from_utf8(&buf[o..o + l])
                            .unwrap_or("")
                            .to_string()
                    })
                })
                .unwrap_or_default();

            if ok {
                EicRoi {
                    id,
                    rt,
                    mz,
                    window: w,
                }
            } else {
                EicRoi {
                    id,
                    rt: 0.0,
                    mz: 0.0,
                    window: 0.0,
                }
            }
        })
        .collect()
}

#[inline]
fn f64ok(v: f64) -> f64 {
    if v.is_finite() { v } else { 0.0 }
}

fn f64_to_u8(v: &[f64]) -> Box<[u8]> {
    let n = core::mem::size_of_val(v);
    let mut out = vec![0u8; n];
    unsafe {
        ptr::copy_nonoverlapping(v.as_ptr() as *const u8, out.as_mut_ptr(), n);
    }
    out.into_boxed_slice()
}

fn write_buf(out: *mut Buf, bytes: Box<[u8]>) {
    let len = bytes.len();
    unsafe {
        ptr::write_unaligned(
            out,
            Buf {
                ptr: Box::into_raw(bytes) as *mut u8,
                len,
            },
        )
    };
}

fn build_fp(opts: *const CPeakPOptions) -> FindPeaksOptions {
    if opts.is_null() {
        return FindPeaksOptions {
            get_boundaries_options: Some(Default::default()),
            filter_peaks_options: Some(Default::default()),
            baseline_options: Some(Default::default()),
        };
    }
    let o = unsafe { *opts };
    FindPeaksOptions {
        get_boundaries_options: Some(Default::default()),
        filter_peaks_options: Some(FilterPeaksOptions {
            integral_threshold: (o.integral_threshold.is_finite() && o.integral_threshold >= 0.0)
                .then_some(o.integral_threshold),
            intensity_threshold: (o.intensity_threshold.is_finite()
                && o.intensity_threshold >= 0.0)
                .then_some(o.intensity_threshold),
            width_threshold: (o.width_threshold > 0).then_some(o.width_threshold as usize),
            noise: (o.noise.is_finite() && o.noise > 0.0).then_some(o.noise),
            auto_noise: Some(o.auto_noise != 0),
            auto_baseline: Some(o.auto_baseline != 0),
            allow_overlap: Some(o.allow_overlap != 0),
            sn_ratio: Some(if o.sn_ratio.is_finite() && o.sn_ratio > 0.0 {
                o.sn_ratio
            } else {
                1.5
            }),
        }),
        baseline_options: Some(BaselineOptions {
            baseline_window: (o.baseline_window > 0).then_some(o.baseline_window as f64),
            baseline_window_factor: (o.baseline_window_factor > 0)
                .then_some(o.baseline_window_factor as usize),
            level: Some(1),
        }),
    }
}
