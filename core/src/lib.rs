#[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
use core::ffi::c_int;
#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use core::ffi::{CStr, c_char, c_int};

use serde::Serialize;
use std::{
    panic::{AssertUnwindSafe, catch_unwind},
    ptr, slice,
    sync::Arc,
};
use utilities::structs::ser_finite_f64;

use ionic::{
    ScanSource, ScanSummary, bin_to_mzml as bin_to_mzml_rs,
    encoder::encode,
    ion::{ByteRange, DecoderConfig, Ion, OwnedIon, plan_open_ranges},
    mzml::structs::MzML,
    parse_mzml as parse_mzml_rs,
};

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use ionic::{
    EncodingConfig, IonWriter, SectionChunkMode,
    encoder::FileEncoderOutput,
    encoder::encode::{DEFAULT_MIN_SPLIT_BYTES, DEFAULT_TARGET_SEGMENT_BYTES},
    mzml::MzmlReader,
    stream_to_ion,
};

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
pub mod url_source;
pub mod utilities;

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use url_source::UrlSource;
use utilities::{
    calculate_baseline::{BaselineOptions, calculate_baseline as calculate_baseline_rs},
    calculate_eic::{
        EicOptions, EicReader, FastError, ScanQuery, TimeUnit,
        calculate_eic as calculate_eic_dispatcher, get_scans as get_scans_rs,
    },
    find_feature::{FindFeatureOptions, find_feature as find_feature_rs},
    find_features::{FindFeaturesOptions, find_features as find_features_rs},
    find_noise_level::find_noise_level as find_noise_level_rs,
    find_peaks::{FindPeaksOptions, PeakFilter, find_peaks as find_peaks_rs},
    get_peak::get_peak as get_peak_rs,
    get_peaks_from_chrom::get_peaks_from_chrom as get_peaks_from_chrom_rs,
    get_peaks_from_eic::get_peaks_from_eic as get_peaks_from_eic_rs,
    structs::{ChromRoi, EicRoi},
    structs::{DataXY, FromTo, Roi},
};

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use utilities::get_features::{AlignmentOptions, get_features as get_features_rs};

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use crate::utilities::find_features::MzTolerance;

#[derive(Serialize)]
struct EicPeakOut<'a> {
    id: &'a str,
    #[serde(serialize_with = "ser_finite_f64")]
    mz: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    ort: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    rt: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    from: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    to: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    intensity: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    integral: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    noise: f64,
}

#[derive(Serialize)]
struct ChromPeakRowOut<'a> {
    index: usize,
    id: &'a str,
    #[serde(serialize_with = "ser_finite_f64")]
    ort: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    rt: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    from: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    to: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    intensity: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    integral: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    total_area: f64,
    timestamp: &'a str,
}

#[derive(Serialize)]
struct FoundFeatureOut<'a> {
    id: &'a str,
    #[serde(serialize_with = "ser_finite_f64")]
    mz: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    rt: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    from: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    to: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    intensity: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    integral: f64,
    n_points: usize,
    #[serde(serialize_with = "ser_finite_f64")]
    noise: f64,
}

pub const MSUTILS_ABI_VERSION: u32 = 2;

#[unsafe(no_mangle)]
pub extern "C" fn msutils_abi_version() -> u32 {
    MSUTILS_ABI_VERSION
}

#[unsafe(no_mangle)]
pub extern "C" fn msutils_sizeof_peak_options() -> usize {
    core::mem::size_of::<CPeakOptions>()
}

const OK: c_int = 0;
const ERR_INVALID_ARGS: c_int = 1;
const ERR_PANIC: c_int = 2;
const ERR_PARSE: c_int = 4;
const ERR_ENCODE: c_int = 5;
const ERR_FAST_PATH: c_int = 6;

fn fast_error_to_code(error: FastError) -> c_int {
    match error {
        FastError::InvalidRequest => ERR_INVALID_ARGS,
        _ => ERR_FAST_PATH,
    }
}

#[repr(C)]
pub struct Buf {
    pub ptr: *mut u8,
    pub len: usize,
}

#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct CPeakOptions {
    pub min_integral: f64,
    pub min_intensity: f64,
    pub min_peak_width_points: c_int,
    pub noise: f64,
    pub auto_noise: c_int,
    pub auto_baseline: c_int,
    pub lambda: c_int,
    pub max_iterations: c_int,
    pub allow_overlap: c_int,
    pub min_snr: f64,
}

// Compile-time guard: if any field is added/removed/reordered and the size
// changes, the build fails immediately. Every wrapper hardcodes 64 bytes
// (Python ctypes layout, WASM packing, C++ static_assert). When this fires:
// 1) Bump MSUTILS_ABI_VERSION, 2) Update the wrapper structs to match.
const _: () = assert!(
    core::mem::size_of::<CPeakOptions>() == 64,
    "CPeakOptions size drifted — bump MSUTILS_ABI_VERSION and update all wrappers"
);

pub trait RangeReader {
    fn read(&self, offset: u64, len: u64, dest: &mut [u8]) -> i32;
}

pub fn read_range<R: RangeReader>(
    reader: &R,
    query: ionic::ion::Query,
) -> ionic::ion::IonResult<ionic::ion::QueryPayload> {
    use ionic::ion::{IonError, QueryPayload};

    let len = query.length();
    if len == 0 {
        return Ok(QueryPayload::new(Vec::new()));
    }
    let len_u32 = u32::try_from(len)
        .map_err(|_| IonError::from("range read: length exceeds transport limit"))?;
    let mut buf = vec![0u8; len_u32 as usize];
    let rc = reader.read(query.offset(), len, &mut buf);
    if rc != 0 {
        return Err(IonError::from(format!("range read failed: {rc}")));
    }
    Ok(QueryPayload::new(buf))
}

#[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
#[link(wasm_import_module = "env")]
unsafe extern "C" {
    fn js_log(ptr: *const u8, len: usize);
    fn range_read(
        source_id: u32,
        offset_lo: u32,
        offset_hi: u32,
        len: u32,
        dest_ptr: *mut u8,
    ) -> i32;
}

#[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
struct WasmRangeReader {
    source_id: u32,
}

#[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
impl RangeReader for WasmRangeReader {
    fn read(&self, offset: u64, len: u64, dest: &mut [u8]) -> i32 {
        let len_u32 = match u32::try_from(len) {
            Ok(value) => value,
            Err(_) => return -2,
        };
        unsafe {
            range_read(
                self.source_id,
                offset as u32,
                (offset >> 32) as u32,
                len_u32,
                dest.as_mut_ptr(),
            )
        }
    }
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

pub enum ParsedFile {
    Full(Box<MzML>),
    Lazy(Box<OwnedIon>),
    Remote(Box<Ion>),
}

impl ParsedFile {
    fn with_mzml<T>(&mut self, f: impl FnOnce(&MzML) -> Result<T, c_int>) -> Result<T, c_int> {
        match self {
            ParsedFile::Full(mzml) => f(mzml.as_ref()),
            ParsedFile::Lazy(file) => {
                let mzml = file.to_mzml().map_err(|_| ERR_PARSE)?;
                f(&mzml)
            }
            ParsedFile::Remote(file) => {
                let mzml = file.to_mzml().map_err(|_| ERR_PARSE)?;
                f(&mzml)
            }
        }
    }
}

impl ScanSource for ParsedFile {
    fn for_each_summary(&mut self, cb: &mut dyn FnMut(usize, ScanSummary)) {
        match self {
            ParsedFile::Full(mzml) => mzml.for_each_summary(cb),
            ParsedFile::Lazy(file) => file.for_each_summary(cb),
            ParsedFile::Remote(file) => file.for_each_summary(cb),
        }
    }

    fn load_scan(&mut self, index: usize, mz: &mut Vec<f64>, intensity: &mut Vec<f64>) -> bool {
        match self {
            ParsedFile::Full(mzml) => mzml.load_scan(index, mz, intensity),
            ParsedFile::Lazy(file) => file.load_scan(index, mz, intensity),
            ParsedFile::Remote(file) => file.load_scan(index, mz, intensity),
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

/// Parse mzML bytes and store the result in `dest` (wasm variant).
///
/// # Safety
/// `data_ptr` must point to `data_len` readable bytes.
/// `dest` must be a valid writable pointer to `*mut ParsedFile`.
/// On success, `*dest` must be freed with `free_mzml`.
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

#[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
#[unsafe(no_mangle)]
pub unsafe extern "C" fn parse_ion_url(
    source_id: u32,
    cache_bytes: usize,
    dest: *mut *mut ParsedFile,
) -> c_int {
    if dest.is_null() {
        return ERR_INVALID_ARGS;
    }
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let reader = WasmRangeReader { source_id };
        let ion = Ion::open_with_query(
            move |query| read_range(&reader, query),
            DecoderConfig {
                max_cached_bytes: cache_bytes,
                ..Default::default()
            },
        )
        .map_err(|_| ERR_PARSE)?;

        unsafe { *dest = Box::into_raw(Box::new(ParsedFile::Remote(Box::new(ion)))) };
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
        let owned = Ion::open_bytes(
            arc,
            DecoderConfig {
                max_cached_bytes: max_cache_size,
                ..Default::default()
            },
        )
        .map_err(|_| ERR_PARSE)?;

        unsafe { *dest = Box::into_raw(Box::new(ParsedFile::Lazy(Box::new(owned)))) };
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => ERR_PANIC,
    }
}

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
#[unsafe(no_mangle)]
pub unsafe extern "C" fn parse_ion_path(
    path_ptr: *const c_char,
    max_cache_size: usize,
    dest: *mut *mut ParsedFile,
) -> c_int {
    if path_ptr.is_null() || dest.is_null() {
        return ERR_INVALID_ARGS;
    }
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let path_text = unsafe { CStr::from_ptr(path_ptr) }
            .to_str()
            .map_err(|_| ERR_INVALID_ARGS)?;
        let file_path = std::path::Path::new(path_text);
        let opened_file = Ion::open_file(
            file_path,
            DecoderConfig {
                max_cached_bytes: max_cache_size,
                ..Default::default()
            },
        )
        .map_err(|_| ERR_PARSE)?;
        unsafe { *dest = Box::into_raw(Box::new(ParsedFile::Lazy(Box::new(opened_file)))) };
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => ERR_PANIC,
    }
}

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
#[unsafe(no_mangle)]
pub unsafe extern "C" fn parse_ion_url(
    url_ptr: *const c_char,
    cache_bytes: usize,
    out: *mut *mut ParsedFile,
) -> c_int {
    if url_ptr.is_null() || out.is_null() {
        return ERR_INVALID_ARGS;
    }
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let url = unsafe { CStr::from_ptr(url_ptr) }
            .to_str()
            .map_err(|_| ERR_INVALID_ARGS)?
            .to_string();
        let source = UrlSource::new(url).map_err(|_| ERR_PARSE)?;
        let ion = Ion::open_with_query(
            move |query| source.read(query),
            DecoderConfig {
                max_cached_bytes: cache_bytes,
                ..Default::default()
            },
        )
        .map_err(|_| ERR_PARSE)?;

        unsafe { *out = Box::into_raw(Box::new(ParsedFile::Remote(Box::new(ion)))) };
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => ERR_PANIC,
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn plan_open(
    header_ptr: *const u8,
    header_len: usize,
    out: *mut Buf,
) -> c_int {
    if header_ptr.is_null() || out.is_null() {
        return ERR_INVALID_ARGS;
    }

    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let header = unsafe { slice::from_raw_parts(header_ptr, header_len) };
        let ranges = plan_open_ranges(header).map_err(|_| ERR_FAST_PATH)?;
        let bytes = pack_byte_ranges(&ranges);
        write_buf(out, bytes);
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => ERR_PANIC,
    }
}

/// Serialize a parsed file to JSON and write it to `out`.
///
/// # Safety
/// `h` must be a valid `ParsedFile` pointer from this library.
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
/// `h` must be a valid `ParsedFile` pointer from this library.
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
/// `h` must be a valid `ParsedFile` pointer from this library.
/// `out` must be a valid writable `Buf` pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn mzml_to_bin(
    h: *mut ParsedFile,
    out: *mut Buf,
    compression_level: u8,
    f32_compress: u8,
) -> c_int {
    if h.is_null() || out.is_null() {
        return ERR_INVALID_ARGS;
    }
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let file = unsafe { &mut *h };
        file.with_mzml(|mzml| {
            let mut bytes = Vec::new();
            encode(mzml, compression_level, f32_compress != 0, &mut bytes)
                .map_err(|_| ERR_ENCODE)?;
            write_buf(out, bytes.into_boxed_slice());
            Ok(())
        })?;
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => ERR_PANIC,
    }
}

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
const ION_BLOCK_SIZE_BYTES: usize = 1024 * 1024;

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
#[unsafe(no_mangle)]
pub unsafe extern "C" fn convert_mzml_file_to_ion_file(
    input_path: *const c_char,
    output_path: *const c_char,
    compression_level: u8,
    force_f32: u8,
    section_on_disk: u8,
) -> c_int {
    if input_path.is_null() || output_path.is_null() {
        return ERR_INVALID_ARGS;
    }
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let input_path = unsafe { CStr::from_ptr(input_path) }
            .to_str()
            .map_err(|_| ERR_INVALID_ARGS)?;
        let output_path = unsafe { CStr::from_ptr(output_path) }
            .to_str()
            .map_err(|_| ERR_INVALID_ARGS)?;

        let mut reader =
            MzmlReader::open(std::path::Path::new(input_path)).map_err(|_| ERR_PARSE)?;
        let mut output =
            FileEncoderOutput::open_for_writing(output_path).map_err(|_| ERR_ENCODE)?;

        let config = EncodingConfig {
            compression_level,
            force_f32: force_f32 != 0,
            uncompressed_block_size: ION_BLOCK_SIZE_BYTES,
            parallel: true,
            section_chunk: if section_on_disk != 0 {
                SectionChunkMode::Disk
            } else {
                SectionChunkMode::Memory
            },
            target_segment_bytes: DEFAULT_TARGET_SEGMENT_BYTES,
            min_split_bytes: DEFAULT_MIN_SPLIT_BYTES,
        };

        {
            let mut writer = IonWriter::begin(&mut output, config).map_err(|_| ERR_ENCODE)?;
            stream_to_ion(&mut reader, &mut writer).map_err(|_| ERR_ENCODE)?;
        }
        output.flush().map_err(|_| ERR_ENCODE)?;
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => ERR_PANIC,
    }
}

/// Find one peak and write the result to `out`.
///
/// # Safety
/// `x_ptr` and `y_ptr` must point to `len` readable `f64` values.
/// `options` must be null or point to a valid `CPeakOptions`.
/// `out` must be a valid writable `Buf` pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn get_peak(
    x_ptr: *const f64,
    y_ptr: *const f64,
    len: usize,
    rt: f64,
    range: f64,
    options: *const CPeakOptions,
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
        let peak = get_peak_rs(
            &data,
            &Roi {
                rt,
                half_width: range,
            },
            Some(build_peak_options(options)),
        );
        let s = serde_json::to_string(&peak.unwrap_or_default()).map_err(|_| ERR_ENCODE)?;
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
/// `h` must be a valid `ParsedFile` pointer from this library.
/// `rts_ptr`, `mzs_ptr`, and `ranges_ptr` must point to `n` readable `f64` values.
/// If IDs are provided, `ids_off` and `ids_len` must point to `n` readable `u32` values,
/// and `ids_buf` must point to `ids_buf_len` readable bytes.
/// `opts` must be null or point to a valid `CPeakOptions`.
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
    opts: *const CPeakOptions,
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
        let mut reader = match file {
            ParsedFile::Full(mzml) => EicReader::Mzml(mzml.as_mut()),
            ParsedFile::Lazy(ion) => EicReader::Ion(ion.as_mut()),
            ParsedFile::Remote(ion) => EicReader::Ion(ion.as_mut()),
        };

        let peaks = get_peaks_from_eic_rs(
            &mut reader,
            FromTo { from, to },
            &items,
            Some(build_peak_options(opts)),
            cores,
        )
        .map_err(fast_error_to_code)?;
        let arr: Vec<EicPeakOut> = peaks
            .iter()
            .map(|(id, ort, mz, p)| EicPeakOut {
                id,
                mz: *mz,
                ort: *ort,
                rt: p.rt,
                from: p.from,
                to: p.to,
                intensity: p.intensity,
                integral: p.integral,
                noise: p.noise,
            })
            .collect();
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
/// `h` must be a valid `ParsedFile` pointer from this library.
/// `sample_indices` must point to `n` readable `u32` values.
/// `target_rts` and `half_widths` must point to `n` readable `f64` values.
/// `opts` must be null or point to a valid `CPeakOptions`.
/// `out` must be a valid writable `Buf` pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn get_peaks_from_chrom(
    h: *mut ParsedFile,
    sample_indices: *const u32,
    target_rts: *const f64,
    half_widths: *const f64,
    n: usize,
    opts: *const CPeakOptions,
    cores: usize,
    out: *mut Buf,
) -> i32 {
    if h.is_null()
        || sample_indices.is_null()
        || target_rts.is_null()
        || half_widths.is_null()
        || out.is_null()
        || n == 0
    {
        return ERR_INVALID_ARGS;
    }
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), i32> {
        let file = unsafe { &mut *h };
        file.with_mzml(|mzml| {
            let sample_indices = unsafe { slice::from_raw_parts(sample_indices, n) };
            let target_rts = unsafe { slice::from_raw_parts(target_rts, n) };
            let half_widths = unsafe { slice::from_raw_parts(half_widths, n) };
            let chroms = &mzml
                .run
                .chromatogram_list
                .as_ref()
                .ok_or(ERR_PARSE)?
                .chromatograms;
            let items: Vec<_> = (0..n)
                .map(|i| {
                    let raw_index = sample_indices[i];
                    if raw_index == u32::MAX {
                        return ChromRoi {
                            id: String::new(),
                            sample_index: usize::MAX,
                            rt: 0.0,
                            half_width: 0.0,
                        };
                    }
                    let sample_index = raw_index as usize;
                    if sample_index >= chroms.len() {
                        return ChromRoi {
                            id: String::new(),
                            sample_index,
                            rt: 0.0,
                            half_width: 0.0,
                        };
                    }
                    ChromRoi {
                        id: chroms[sample_index].id.clone(),
                        sample_index,
                        rt: target_rts[i],
                        half_width: half_widths[i],
                    }
                })
                .collect();
            let list = get_peaks_from_chrom_rs(mzml, &items, Some(build_peak_options(opts)), cores)
                .ok_or(ERR_PARSE)?;
            let out_arr: Vec<_> = list
                .iter()
                .map(|row| ChromPeakRowOut {
                    index: row.index,
                    id: &row.id,
                    ort: row.target_rt,
                    rt: row.peak_rt,
                    from: row.from_rt,
                    to: row.to_rt,
                    intensity: row.intensity,
                    integral: row.area,
                    total_area: row.total_area,
                    timestamp: &row.timestamp,
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
/// `opts` must be null or point to a valid `CPeakOptions`.
/// `out` must be a valid writable `Buf` pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn find_peaks(
    x_ptr: *const f64,
    y_ptr: *const f64,
    len: usize,
    opts: *const CPeakOptions,
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
            Some(build_peak_options(opts)),
        );
        write_buf(
            out,
            serde_json::to_string(&peaks)
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

/// Compute the noise intensity for an intensity array.
///
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

fn plan_eic_ranges(
    ion: &mut Ion,
    target_mz: f64,
    from: f64,
    to: f64,
    ppm: f64,
    mz_tol: f64,
) -> Result<Vec<ByteRange>, FastError> {
    if !target_mz.is_finite() || target_mz <= 0.0 {
        return Err(FastError::InvalidRequest);
    }

    let options = EicOptions {
        ppm_tolerance: ppm,
        mz_tolerance: mz_tol,
        ..Default::default()
    };

    let tolerance = utilities::calculate_eic::mz_tolerance_for(target_mz, options);
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(FastError::InvalidRequest);
    }

    ion.require_spectrum_bounds().map_err(FastError::from)?;

    let mz_from = target_mz - tolerance;
    let mz_to = target_mz + tolerance;
    let rt_from = from.min(to);
    let rt_to = from.max(to);

    let mut scan_indices = Vec::new();
    ion.for_each_summary(&mut |scan_index, summary| {
        if summary.rt >= rt_from && summary.rt <= rt_to && summary.ms_level == 1 {
            scan_indices.push(scan_index);
        }
    });

    let mut ranges = Vec::new();
    for scan_index in scan_indices {
        let scan_ranges = ion
            .spectrum_mz_window_block_ranges_strict(scan_index, mz_from, mz_to)
            .map_err(FastError::from)?;
        ranges.extend(scan_ranges);
    }

    sort_and_dedup_ranges(&mut ranges);
    Ok(ranges)
}

/// Calculate an EIC and write `x` and `y` to `out_x` and `out_y`.
///
/// # Safety
/// `h` must be a valid `ParsedFile` pointer from this library.
/// `out_x` and `out_y` must be valid writable `Buf` pointers.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn calculate_eic(
    h: *mut ParsedFile,
    target: f64,
    from: f64,
    to: f64,
    ppm: f64,
    mz_tol: f64,
    out_x: *mut Buf,
    out_y: *mut Buf,
) -> c_int {
    if h.is_null() || out_x.is_null() || out_y.is_null() {
        return ERR_INVALID_ARGS;
    }
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let file = unsafe { &mut *h };

        let mut reader = match file {
            ParsedFile::Full(mzml) => EicReader::Mzml(mzml.as_mut()),
            ParsedFile::Lazy(ion) => EicReader::Ion(ion.as_mut()),
            ParsedFile::Remote(ion) => EicReader::Ion(ion.as_mut()),
        };

        let eic = calculate_eic_dispatcher(
            &mut reader,
            target,
            FromTo { from, to },
            EicOptions {
                ppm_tolerance: ppm,
                mz_tolerance: mz_tol,
                ..Default::default()
            },
        )
        .map_err(fast_error_to_code)?;

        write_buf(out_x, f64_to_u8(&eic.x));
        write_buf(out_y, f64_to_u8(&eic.y));
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => ERR_PANIC,
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn plan_eic(
    h: *mut ParsedFile,
    target: f64,
    from: f64,
    to: f64,
    ppm: f64,
    mz_tol: f64,
    out: *mut Buf,
) -> c_int {
    if h.is_null() || out.is_null() {
        return ERR_INVALID_ARGS;
    }

    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let file = unsafe { &mut *h };

        let ranges = match file {
            ParsedFile::Lazy(ion) => plan_eic_ranges(ion.as_mut(), target, from, to, ppm, mz_tol),
            ParsedFile::Remote(ion) => plan_eic_ranges(ion.as_mut(), target, from, to, ppm, mz_tol),
            ParsedFile::Full(_) => Err(FastError::UnsupportedBackend),
        }
        .map_err(fast_error_to_code)?;

        let bytes = pack_byte_ranges(&ranges);
        write_buf(out, bytes);
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => ERR_PANIC,
    }
}

/// Get scans by query and write the result to `out`.
///
/// `query_type`: 0=RtRange, 1=ClosestRt, 2=MzRange, 3=ClosestMz.
/// `from_value` and `to_value` are the range bounds for range queries (types 0 and 2).
/// For point queries (types 1 and 3), only `from_value` is used; `to_value` is ignored.
///
/// # Safety
/// `h` must be a valid `ParsedFile` pointer from this library.
/// `out` must be a valid writable `Buf` pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn get_scans(
    h: *mut ParsedFile,
    query_type: u8,
    from_value: f64,
    to_value: f64,
    level: u8,
    out: *mut Buf,
) -> c_int {
    if h.is_null() || out.is_null() {
        return ERR_INVALID_ARGS;
    }
    let query = match query_type {
        0 => {
            if !from_value.is_finite() || !to_value.is_finite() {
                return ERR_INVALID_ARGS;
            }
            ScanQuery::RtRange(FromTo {
                from: from_value,
                to: to_value,
            })
        }
        1 => {
            if !from_value.is_finite() {
                return ERR_INVALID_ARGS;
            }
            ScanQuery::ClosestRt(from_value)
        }
        2 => {
            if !from_value.is_finite() || !to_value.is_finite() {
                return ERR_INVALID_ARGS;
            }
            ScanQuery::MzRange(FromTo {
                from: from_value,
                to: to_value,
            })
        }
        3 => {
            if !from_value.is_finite() || from_value <= 0.0 {
                return ERR_INVALID_ARGS;
            }
            ScanQuery::ClosestMz(from_value)
        }
        _ => return ERR_INVALID_ARGS,
    };
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let (_, scans) = get_scans_rs(unsafe { &mut *h }, query, TimeUnit::Minutes, level);
        write_scans_json(out, &scans).map_err(|_| ERR_PARSE)
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => ERR_PANIC,
    }
}

/// Find features across all samples in `dir` and write the result to `out`.
///
/// # Safety
/// `dir` must point to a valid NUL-terminated C string.
/// `peak_opts` must be null or point to a valid `CPeakOptions`.
/// `out` must be a valid writable `Buf` pointer.
#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
#[unsafe(no_mangle)]
pub unsafe extern "C" fn get_features(
    dir: *const c_char,
    from: f64,
    to: f64,
    eic_ppm: f64,
    eic_mz: f64,
    grid_min_mz: f64,
    grid_max_mz: f64,
    grid_step: f64,
    align_ppm: f64,
    align_mz_absolute: f64,
    align_rt_tolerance: f64,
    min_samples: c_int,
    peak_opts: *const CPeakOptions,
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

        let mut feature_opts = FindFeaturesOptions::default();
        if eic_ppm.is_finite() && eic_ppm >= 0.0 {
            feature_opts.final_eic_options.ppm_tolerance = eic_ppm;
        }
        if eic_mz.is_finite() && eic_mz >= 0.0 {
            feature_opts.final_eic_options.mz_tolerance = eic_mz;
        }
        if grid_min_mz.is_finite() {
            feature_opts.mz_scan_grid.min_mz = grid_min_mz;
        }
        if grid_max_mz.is_finite() {
            feature_opts.mz_scan_grid.max_mz = grid_max_mz;
        }
        if grid_step > 0.0 {
            feature_opts.mz_scan_grid.step = grid_step;
        }
        let peak_options = build_peak_options(peak_opts);
        feature_opts.peak_options = peak_options.clone();
        feature_opts.use_gpu = use_gpu != 0;
        feature_opts.batch_size = if batch_size > 0 {
            Some(batch_size as usize)
        } else {
            None
        };

        let alignment_opts = AlignmentOptions {
            mz_tolerance: MzTolerance {
                mz_absolute: if align_mz_absolute > 0.0 {
                    align_mz_absolute
                } else {
                    0.003
                },
                ppm: if align_ppm > 0.0 { align_ppm } else { 5.0 },
            },
            rt_tolerance: if align_rt_tolerance > 0.0 {
                align_rt_tolerance
            } else {
                0.05
            },
            min_samples: if min_samples > 0 {
                min_samples as usize
            } else {
                1
            },
            eic_options: feature_opts.final_eic_options,
            peak_options: Some(peak_options),
        };

        let feats = get_features_rs(
            path,
            FromTo { from, to },
            feature_opts,
            alignment_opts,
            if cores > 0 { cores as usize } else { 1 },
        )
        .map_err(|_| ERR_FAST_PATH)?;

        write_buf(
            out,
            serde_json::to_string(&feats)
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

/// Find features in a single sample and write the result to `out`.
///
/// # Safety
/// `h` must be a valid `ParsedFile` pointer from this library.
/// `peak_opts` must be null or point to a valid `CPeakOptions`.
/// `out` must be a valid writable `Buf` pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn find_features(
    h: *mut ParsedFile,
    from: f64,
    to: f64,
    eic_ppm: f64,
    eic_mz: f64,
    grid_min_mz: f64,
    grid_max_mz: f64,
    grid_step: f64,
    peak_opts: *const CPeakOptions,
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
            opts.final_eic_options.ppm_tolerance = eic_ppm;
        }
        if eic_mz.is_finite() && eic_mz >= 0.0 {
            opts.final_eic_options.mz_tolerance = eic_mz;
        }
        if grid_min_mz.is_finite() {
            opts.mz_scan_grid.min_mz = grid_min_mz;
        }
        if grid_max_mz.is_finite() {
            opts.mz_scan_grid.max_mz = grid_max_mz;
        }
        if grid_step > 0.0 {
            opts.mz_scan_grid.step = grid_step;
        }
        opts.peak_options = build_peak_options(peak_opts);
        #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
        {
            opts.use_gpu = use_gpu != 0;
            opts.batch_size = if batch_size > 0 {
                Some(batch_size as usize)
            } else {
                None
            };
        }

        let mut reader = match file {
            ParsedFile::Full(mzml) => EicReader::Mzml(mzml.as_mut()),
            ParsedFile::Lazy(ion) => EicReader::Ion(ion.as_mut()),
            ParsedFile::Remote(ion) => EicReader::Ion(ion.as_mut()),
        };

        let feats = find_features_rs(
            &mut reader,
            FromTo { from, to },
            Some(opts),
            if cores > 0 { cores as usize } else { 1 },
        )
        .map_err(|_| ERR_FAST_PATH)?;

        write_buf(
            out,
            serde_json::to_string(&feats)
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
    lambda: c_int,
    max_iterations: c_int,
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
                lambda: if lambda > 0 {
                    Some(lambda as f64)
                } else {
                    def.lambda
                },
                max_iterations: if max_iterations > 0 {
                    Some(max_iterations as usize)
                } else {
                    def.max_iterations
                },
                edge_slope_level: Some(1),
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

/// Find a targeted feature for each ROI and write the result to `out`.
///
/// # Safety
/// `h` must be a valid `ParsedFile` pointer from this library.
/// `target_rts`, `target_mzs`, and `half_widths` must point to `n` readable `f64` values.
/// If IDs are provided, `ids_off` and `ids_len` must point to `n` readable `u32` values,
/// and `ids_buf` must point to `ids_buf_len` readable bytes.
/// `peak_opts` must be null or point to a valid `CPeakOptions`.
/// `out` must be a valid writable `Buf` pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn find_feature(
    h: *const ParsedFile,
    target_rts: *const f64,
    target_mzs: *const f64,
    half_widths: *const f64,
    ids_off: *const u32,
    ids_len: *const u32,
    ids_buf: *const u8,
    ids_buf_len: usize,
    n: usize,
    cores: usize,
    seed_ppm: f64,
    seed_mz: f64,
    final_ppm: f64,
    final_mz: f64,
    peak_opts: *const CPeakOptions,
    out: *mut Buf,
) -> c_int {
    if h.is_null()
        || target_rts.is_null()
        || target_mzs.is_null()
        || half_widths.is_null()
        || out.is_null()
        || n == 0
    {
        return ERR_INVALID_ARGS;
    }
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let file = unsafe { &mut *(h as *mut ParsedFile) };
        let target_rts = unsafe { slice::from_raw_parts(target_rts, n) };
        let target_mzs = unsafe { slice::from_raw_parts(target_mzs, n) };
        let half_widths = unsafe { slice::from_raw_parts(half_widths, n) };
        let mut seed_eic = EicOptions {
            ppm_tolerance: 10.0,
            mz_tolerance: 0.003,
            ..Default::default()
        };
        if seed_ppm.is_finite() && seed_ppm >= 0.0 {
            seed_eic.ppm_tolerance = seed_ppm;
        }
        if seed_mz.is_finite() && seed_mz >= 0.0 {
            seed_eic.mz_tolerance = seed_mz;
        }
        let mut final_eic = EicOptions {
            ppm_tolerance: 20.0,
            mz_tolerance: 0.005,
            ..Default::default()
        };
        if final_ppm.is_finite() && final_ppm >= 0.0 {
            final_eic.ppm_tolerance = final_ppm;
        }
        if final_mz.is_finite() && final_mz >= 0.0 {
            final_eic.mz_tolerance = final_mz;
        }
        let rois = build_eic_rois(
            target_rts,
            target_mzs,
            half_widths,
            ids_off,
            ids_len,
            ids_buf,
            ids_buf_len,
        );
        let refs: Vec<&EicRoi> = rois.iter().collect();
        let results = find_feature_rs(
            file,
            &refs,
            cores,
            Some(FindFeatureOptions {
                seed_eic_options: Some(seed_eic),
                final_eic_options: Some(final_eic),
                peak_options: Some(build_peak_options(peak_opts)),
            }),
        );
        let arr: Vec<FoundFeatureOut> = results
            .iter()
            .enumerate()
            .map(|(i, r)| match r {
                Some(f) => FoundFeatureOut {
                    id: &f.id,
                    mz: f.mz,
                    rt: f.rt,
                    from: f.peak.from,
                    to: f.peak.to,
                    intensity: f.peak.intensity,
                    integral: f.peak.integral,
                    n_points: f.peak.n_points,
                    noise: f.peak.noise,
                },
                None => FoundFeatureOut {
                    id: &rois[i].id,
                    mz: 0.0,
                    rt: 0.0,
                    from: 0.0,
                    to: 0.0,
                    intensity: 0.0,
                    integral: 0.0,
                    n_points: 0,
                    noise: 0.0,
                },
            })
            .collect();
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

fn write_scans_json(
    out: *mut Buf,
    scans: &[utilities::calculate_eic::CentroidScan],
) -> Result<(), serde_json::Error> {
    write_buf(
        out,
        serde_json::to_string(scans)?
            .into_bytes()
            .into_boxed_slice(),
    );
    Ok(())
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
                    half_width: w,
                }
            } else {
                EicRoi {
                    id,
                    rt: 0.0,
                    mz: 0.0,
                    half_width: 0.0,
                }
            }
        })
        .collect()
}

#[inline]
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

fn push_u64_le(bytes: &mut Vec<u8>, value: u64) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn pack_byte_ranges(ranges: &[ByteRange]) -> Box<[u8]> {
    let mut bytes = Vec::with_capacity(ranges.len() * 16);
    for range in ranges {
        push_u64_le(&mut bytes, range.offset);
        push_u64_le(&mut bytes, range.length);
    }
    bytes.into_boxed_slice()
}

fn sort_and_dedup_ranges(ranges: &mut Vec<ByteRange>) {
    ranges.sort_unstable_by_key(|range| (range.offset, range.length));
    ranges.dedup_by_key(|range| (range.offset, range.length));
}

fn build_peak_options(opts: *const CPeakOptions) -> FindPeaksOptions {
    if opts.is_null() {
        return FindPeaksOptions {
            boundaries: Some(Default::default()),
            filter: Some(Default::default()),
            baseline: Some(Default::default()),
            artifact_filter: Some(Default::default()),
        };
    }
    let o = unsafe { *opts };
    FindPeaksOptions {
        boundaries: Some(Default::default()),
        filter: Some(PeakFilter {
            min_integral: (o.min_integral.is_finite() && o.min_integral >= 0.0)
                .then_some(o.min_integral),
            min_intensity: (o.min_intensity.is_finite() && o.min_intensity >= 0.0)
                .then_some(o.min_intensity),
            min_peak_width_points: (o.min_peak_width_points > 0)
                .then_some(o.min_peak_width_points as usize),
            noise: (o.noise.is_finite() && o.noise > 0.0 && o.auto_noise == 0).then_some(o.noise),
            auto_noise: Some(o.auto_noise != 0),
            auto_baseline: Some(o.auto_baseline != 0),
            allow_overlap: Some(o.allow_overlap != 0),
            min_snr: Some(if o.min_snr.is_finite() && o.min_snr > 0.0 {
                o.min_snr
            } else {
                1.5
            }),
        }),
        baseline: Some(BaselineOptions {
            lambda: (o.lambda > 0).then_some(o.lambda as f64),
            max_iterations: (o.max_iterations > 0).then_some(o.max_iterations as usize),
            edge_slope_level: Some(1),
        }),
        artifact_filter: Some(Default::default()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ionic::mzml::structs::{
        BinaryData, BinaryDataArray, BinaryDataArrayList, CvParam, NumericType, Run, Scan,
        ScanList, Spectrum, SpectrumList,
    };

    struct FakeReader {
        data: Vec<u8>,
        call_count: std::sync::atomic::AtomicUsize,
        return_code: i32,
    }

    impl RangeReader for FakeReader {
        fn read(&self, _offset: u64, len: u64, dest: &mut [u8]) -> i32 {
            if len == 0 {
                panic!("read called with len=0");
            }
            self.call_count
                .fetch_add(1, std::sync::atomic::Ordering::SeqCst);
            let copy_len = std::cmp::min(len as usize, dest.len().min(self.data.len()));
            dest[..copy_len].copy_from_slice(&self.data[..copy_len]);
            self.return_code
        }
    }

    #[test]
    fn zero_len_returns_empty_and_never_calls_read() {
        use ionic::ion::Query;

        let fake = FakeReader {
            data: vec![1, 2, 3, 4, 5],
            call_count: std::sync::atomic::AtomicUsize::new(0),
            return_code: 0,
        };

        let query = Query::new(100, 0);
        let result = read_range(&fake, query).unwrap();
        let value = result.bytes().to_vec();

        assert_eq!(value.len(), 0);
        assert_eq!(fake.call_count.load(std::sync::atomic::Ordering::SeqCst), 0);
    }

    #[test]
    fn non_zero_read_calls_reader_once() {
        use ionic::ion::Query;

        let fake = FakeReader {
            data: vec![10, 20, 30, 40, 50],
            call_count: std::sync::atomic::AtomicUsize::new(0),
            return_code: 0,
        };

        let query = Query::new(0, 3);
        let _result = read_range(&fake, query).unwrap();

        assert_eq!(fake.call_count.load(std::sync::atomic::Ordering::SeqCst), 1);
    }

    #[test]
    fn negative_return_code_becomes_error() {
        use ionic::ion::Query;

        let fake = FakeReader {
            data: vec![10, 20, 30],
            call_count: std::sync::atomic::AtomicUsize::new(0),
            return_code: -1,
        };

        let query = Query::new(0, 2);
        let result = read_range(&fake, query);

        assert!(result.is_err());
    }

    #[test]
    fn len_exceeding_u32_max_returns_error_without_allocation() {
        use ionic::ion::Query;

        let fake = FakeReader {
            data: vec![1, 2, 3],
            call_count: std::sync::atomic::AtomicUsize::new(0),
            return_code: 0,
        };

        let query = Query::new(0, u64::from(u32::MAX) + 1);
        let result = read_range(&fake, query);

        assert!(result.is_err());
        assert_eq!(fake.call_count.load(std::sync::atomic::Ordering::SeqCst), 0);
    }

    fn decode_byte_ranges(bytes: &[u8]) -> Vec<(u64, u64)> {
        bytes
            .chunks_exact(16)
            .map(|chunk| {
                let offset = u64::from_le_bytes(chunk[0..8].try_into().unwrap());
                let length = u64::from_le_bytes(chunk[8..16].try_into().unwrap());
                (offset, length)
            })
            .collect()
    }

    fn new_buf_out() -> *mut Buf {
        Box::into_raw(Box::new(Buf {
            ptr: std::ptr::null_mut(),
            len: 0,
        }))
    }

    fn drop_buf_out(out: *mut Buf) {
        unsafe {
            if !(*out).ptr.is_null() {
                free_((*out).ptr, (*out).len);
            }
            drop(Box::from_raw(out));
        }
    }

    fn cv_param(accession: &str, value: Option<&str>) -> CvParam {
        CvParam {
            accession: Some(accession.to_string()),
            name: accession.to_string(),
            value: value.map(str::to_string),
            ..Default::default()
        }
    }

    fn cv_param_with_unit(accession: &str, value: &str, unit_accession: &str) -> CvParam {
        CvParam {
            accession: Some(accession.to_string()),
            name: accession.to_string(),
            value: Some(value.to_string()),
            unit_accession: Some(unit_accession.to_string()),
            ..Default::default()
        }
    }

    fn binary_array(role_accession: &str, values: Vec<f64>) -> BinaryDataArray {
        BinaryDataArray {
            array_length: Some(values.len()),
            cv_params: vec![
                cv_param(role_accession, None),
                cv_param("MS:1000523", None),
                cv_param("MS:1000576", None),
            ],
            numeric_type: Some(NumericType::Float64),
            binary: Some(BinaryData::F64(values)),
            ..Default::default()
        }
    }

    fn spectrum(index: usize, rt: f64) -> Spectrum {
        Spectrum {
            id: format!("scan={}", index + 1),
            index: Some(index as u32),
            default_array_length: Some(4),
            ms_level: Some(1),
            cv_params: vec![cv_param("MS:1000511", Some("1"))],
            scan_list: Some(ScanList {
                count: Some(1),
                scans: vec![Scan {
                    cv_params: vec![cv_param_with_unit(
                        "MS:1000016",
                        &rt.to_string(),
                        "UO:0000031",
                    )],
                    ..Default::default()
                }],
                ..Default::default()
            }),
            binary_data_array_list: Some(BinaryDataArrayList {
                count: Some(2),
                binary_data_arrays: vec![
                    binary_array("MS:1000514", vec![100.0, 499.999, 500.001, 900.0]),
                    binary_array("MS:1000515", vec![10.0, 1000.0, 900.0, 20.0]),
                ],
            }),
            ..Default::default()
        }
    }

    fn current_a3_ion_bytes() -> Vec<u8> {
        let mzml = MzML {
            run: Run {
                id: "test".to_string(),
                spectrum_list: Some(SpectrumList {
                    count: Some(3),
                    spectra: vec![spectrum(0, 1.0), spectrum(1, 2.0), spectrum(2, 3.0)],
                    ..Default::default()
                }),
                ..Default::default()
            },
            ..Default::default()
        };

        let mut bytes = Vec::new();
        encode(&mzml, 0, false, &mut bytes).expect("ion encode should succeed");
        bytes
    }

    fn open_current_a3_ion() -> *mut ParsedFile {
        let bytes = current_a3_ion_bytes();
        let bytes_arc = Arc::from(bytes.into_boxed_slice());
        let ion = OwnedIon::open_bytes(bytes_arc, DecoderConfig::default())
            .expect("OwnedIon::open_bytes failed");
        Box::into_raw(Box::new(ParsedFile::Lazy(Box::new(ion))))
    }

    fn open_old_ion_without_a3() -> *mut ParsedFile {
        let bytes = include_bytes!("../tests/fixtures/no_a3_test.ion");
        let bytes_arc = Arc::from(bytes.as_slice());
        let ion = OwnedIon::open_bytes(bytes_arc, DecoderConfig::default())
            .expect("OwnedIon::open_bytes failed");
        Box::into_raw(Box::new(ParsedFile::Lazy(Box::new(ion))))
    }

    #[test]
    fn plan_open_rejects_null_header_ptr() {
        let out = new_buf_out();
        let code = unsafe { plan_open(std::ptr::null(), 1024, out) };
        assert_eq!(code, ERR_INVALID_ARGS);
        drop_buf_out(out);
    }

    #[test]
    fn plan_open_rejects_null_out() {
        let code = unsafe { plan_open(b"header".as_ptr(), 6, std::ptr::null_mut()) };
        assert_eq!(code, ERR_INVALID_ARGS);
    }

    #[test]
    fn plan_open_rejects_short_header() {
        let out = new_buf_out();
        let code = unsafe { plan_open(b"x".as_ptr(), 1, out) };
        assert_eq!(code, ERR_FAST_PATH);
        drop_buf_out(out);
    }

    #[test]
    fn plan_open_returns_buffer_divisible_by_16() {
        let out = new_buf_out();
        let header = vec![0u8; 2048];
        let _code = unsafe { plan_open(header.as_ptr(), header.len(), out) };
        let len = unsafe { (*out).len };
        assert_eq!(len % 16, 0);
        drop_buf_out(out);
    }

    #[test]
    fn pack_byte_ranges_writes_little_endian_pairs() {
        let ranges = vec![
            ByteRange {
                offset: 100,
                length: 32,
            },
            ByteRange {
                offset: 5000,
                length: 64,
            },
        ];

        let bytes = pack_byte_ranges(&ranges);

        assert_eq!(bytes.len(), 32);
        assert_eq!(decode_byte_ranges(&bytes), vec![(100, 32), (5000, 64)]);
    }

    #[test]
    fn plan_eic_rejects_null_handle() {
        let out = new_buf_out();
        let code = unsafe { plan_eic(std::ptr::null_mut(), 500.0, 0.0, 10.0, 20.0, 0.005, out) };
        assert_eq!(code, ERR_INVALID_ARGS);
        drop_buf_out(out);
    }

    #[test]
    fn plan_eic_rejects_null_out() {
        let code = unsafe {
            plan_eic(
                std::ptr::null_mut(),
                500.0,
                0.0,
                10.0,
                20.0,
                0.005,
                std::ptr::null_mut(),
            )
        };
        assert_eq!(code, ERR_INVALID_ARGS);
    }

    #[test]
    fn plan_eic_rejects_invalid_target_mz() {
        let handle = open_current_a3_ion();
        let out = new_buf_out();
        let code = unsafe { plan_eic(handle, -1.0, 0.0, 10.0, 20.0, 0.005, out) };
        assert_eq!(code, ERR_INVALID_ARGS);
        drop_buf_out(out);
        unsafe { free_mzml(handle) };
    }

    #[test]
    fn plan_eic_rejects_zero_tolerance() {
        let handle = open_current_a3_ion();
        let out = new_buf_out();
        let code = unsafe { plan_eic(handle, 500.0, 0.0, 10.0, 0.0, 0.0, out) };
        assert_eq!(code, ERR_INVALID_ARGS);
        drop_buf_out(out);
        unsafe { free_mzml(handle) };
    }

    #[test]
    fn plan_eic_returns_ranges_for_current_a3_ion() {
        let handle = open_current_a3_ion();
        let out = new_buf_out();

        let code = unsafe { plan_eic(handle, 500.0, 0.0, 10.0, 20.0, 0.005, out) };

        assert_eq!(code, OK);

        let bytes = unsafe { slice::from_raw_parts((*out).ptr, (*out).len).to_vec() };
        assert_eq!(bytes.len() % 16, 0);

        let ranges = decode_byte_ranges(&bytes);
        assert!(!ranges.is_empty());
        assert!(ranges.iter().all(|(_, length)| *length > 0));

        drop_buf_out(out);
        unsafe { free_mzml(handle) };
    }

    #[test]
    fn plan_eic_requires_spectrum_bounds_for_ion() {
        let handle = open_old_ion_without_a3();
        let out = new_buf_out();

        let code = unsafe { plan_eic(handle, 500.0, 0.0, 9999.0, 20.0, 0.005, out) };

        assert_eq!(code, ERR_FAST_PATH);

        drop_buf_out(out);
        unsafe { free_mzml(handle) };
    }
}
