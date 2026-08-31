#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use core::ffi::{CStr, c_char, c_int, c_void};
#[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
use core::ffi::{c_int, c_void};
use std::{
    panic::{AssertUnwindSafe, catch_unwind},
    ptr, slice,
    sync::Arc,
};

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use ionic::{
    DEFAULT_MZ_WINDOW, IonWriter,
    ion::{FileWriter, SectionStorage},
    mzml::MzmlReader,
};
use ionic::{
    ScanSource, ScanSummary, WriteOptions, bin_to_mzml as bin_to_mzml_rs, coalesce_byte_ranges,
    ion::{ByteRange, BytesSource, CallbackSource, IonReader, ReadBytes, ReadOptions, open_ranges},
    mzml::structs::MzML,
    parse_mzml as parse_mzml_rs, write_mzml_to_ion,
};

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use crate::utilities::find_features::MzTolerance;
#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use crate::utilities::get_features::{AlignmentOptions, get_features as get_features_rs};
use crate::utilities::{
    calculate_baseline::{BaselineOptions, calculate_baseline as calculate_baseline_rs},
    calculate_eic::{
        EicOptions, EicReader, FastError, ScanQuery, TimeUnit,
        calculate_eic as calculate_eic_dispatcher, get_scans as get_scans_rs, plan_eic_ranges,
        plan_scan_ranges,
    },
    find_feature::{FindFeatureOptions, find_feature as find_feature_rs},
    find_features::{FindFeaturesOptions, find_features as find_features_rs},
    find_noise_level::find_noise_level as find_noise_level_rs,
    find_peaks::{ArtifactFilter, FindPeaksOptions, PeakFilter, find_peaks as find_peaks_rs},
    bridge::*,
    fit_peak::{
        PeakParameters, PeakSeed, PeakShape, draw_peak as draw_peak_rs, fit_peak as fit_peak_rs,
    },
    get_peak::get_peak as get_peak_rs,
    get_peaks_from_chrom::get_peaks_from_chrom as get_peaks_from_chrom_rs,
    get_peaks_from_eic::get_peaks_from_eic as get_peaks_from_eic_rs,
    mz_estimator::MzEstimatorKind,
    structs::{DataXY, FromTo, Roi},
};

struct ChromPeakRowOut<'a> {
    index: usize,
    id: &'a str,
    ort: f64,
    rt: f64,
    from: f64,
    to: f64,
    intensity: f64,
    integral: f64,
    total_area: f64,
    timestamp: &'a str,
}

struct FoundFeatureOut<'a> {
    id: &'a str,
    mz: f64,
    rt: f64,
    from: f64,
    to: f64,
    intensity: f64,
    integral: f64,
    n_points: usize,
    noise: f64,
}

pub const QUANTION_ABI_VERSION: u32 = 1;

#[unsafe(no_mangle)]
pub extern "C" fn quantion_abi_version() -> u32 {
    QUANTION_ABI_VERSION
}

#[unsafe(no_mangle)]
pub extern "C" fn quantion_sizeof_peak_options() -> usize {
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
    pub shape: c_int,
    pub noise: f64,
    pub auto_noise: c_int,
    pub auto_baseline: c_int,
    pub lambda: c_int,
    pub max_iterations: c_int,
    pub allow_overlap: c_int,
    pub min_snr: f64,
    pub min_r2: f64,
    pub kernel_size: c_int,
}

const _: () = assert!(
    core::mem::size_of::<CPeakOptions>() == 80,
    "CPeakOptions size drifted — bump QUANTION_ABI_VERSION and update all wrappers"
);

pub trait RangeReader {
    fn read(&self, offset: u64, len: u64, dest: &mut [u8]) -> i32;
}

pub fn read_range<R: RangeReader>(
    reader: &R,
    range: ionic::ion::ByteRange,
) -> ionic::ion::IonResult<Vec<u8>> {
    use ionic::ion::IonError;

    let len = range.length;
    if len == 0 {
        return Ok(Vec::new());
    }
    let len_u32 = u32::try_from(len)
        .map_err(|_| IonError::from("range read: length exceeds transport limit"))?;
    let mut buf = vec![0u8; len_u32 as usize];
    let rc = reader.read(range.offset, len, &mut buf);
    if rc != 0 {
        return Err(IonError::from(format!("range read failed: {rc}")));
    }
    Ok(buf)
}

/// A function the caller supplies so the library can ask for a slice of the file.
///
/// It must write `length` bytes taken from `offset` into `dest`, and return 0.
/// Any other value means the read failed.
/// It is always called on the thread that called the entry point.
pub type ReadRangeCallback =
    unsafe extern "C" fn(context: *mut c_void, offset: u64, length: u64, dest: *mut u8) -> i32;

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
struct CallerRangeReader {
    read: ReadRangeCallback,
    context: *mut c_void,
}

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
unsafe impl Send for CallerRangeReader {}

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
unsafe impl Sync for CallerRangeReader {}

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
impl RangeReader for CallerRangeReader {
    fn read(&self, offset: u64, len: u64, dest: &mut [u8]) -> i32 {
        unsafe { (self.read)(self.context, offset, len, dest.as_mut_ptr()) }
    }
}

/// Open an ion file whose bytes the caller supplies on demand.
///
/// The library never reads a file or a socket itself: every time it needs a slice
/// it calls `read`. Use this for a file that lives somewhere the library cannot
/// reach, such as a URL.
///
/// # Safety
/// `read` must stay callable, and `context` must stay valid, until the returned
/// handle is freed with `free_mzml`.
#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
#[unsafe(no_mangle)]
pub unsafe extern "C" fn parse_ion_source(
    read: Option<
        unsafe extern "C" fn(context: *mut c_void, offset: u64, length: u64, dest: *mut u8) -> i32,
    >,
    context: *mut c_void,
    max_cache_size: usize,
    out: *mut *mut ParsedFile,
) -> c_int {
    let Some(read) = read else {
        return ERR_INVALID_ARGS;
    };
    if out.is_null() {
        return ERR_INVALID_ARGS;
    }
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let reader = CallerRangeReader { read, context };
        let ion = IonReader::open_source(
            Arc::new(CallbackSource::new(move |range| read_range(&reader, range)))
                as Arc<dyn ReadBytes>,
            ReadOptions {
                max_cached_bytes: max_cache_size,
                ..Default::default()
            },
        )
        .map_err(|_| ERR_PARSE)?;

        unsafe { *out = Box::into_raw(Box::new(ParsedFile::new(FileSource::Remote(Box::new(ion))))) };
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => ERR_PANIC,
    }
}

#[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
#[link(wasm_import_module = "env")]
unsafe extern "C" {
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

#[unsafe(no_mangle)]
pub extern "C" fn alloc(size: usize) -> *mut u8 {
    if size == 0 {
        return core::ptr::null_mut();
    }
    catch_unwind(AssertUnwindSafe(|| {
        let mut v = Vec::<u8>::with_capacity(size);
        let p = v.as_mut_ptr();
        core::mem::forget(v);
        p
    }))
    .unwrap_or(core::ptr::null_mut())
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

pub enum FileSource {
    Full(Box<MzML>),
    Lazy(Box<IonReader>),
    Remote(Box<IonReader>),
}

pub struct ParsedFile {
    source: FileSource,
    broken: bool,
}

impl ParsedFile {
    pub(crate) fn new(source: FileSource) -> Self {
        Self {
            source,
            broken: false,
        }
    }

    fn source_mut(&mut self) -> &mut FileSource {
        &mut self.source
    }

    fn with_mzml<T>(&mut self, f: impl FnOnce(&MzML) -> Result<T, c_int>) -> Result<T, c_int> {
        match &mut self.source {
            FileSource::Full(mzml) => f(mzml.as_ref()),
            FileSource::Lazy(file) => {
                let mzml = file.to_mzml().map_err(|_| ERR_PARSE)?;
                f(&mzml)
            }
            FileSource::Remote(ion) => {
                let mzml = ion.to_mzml().map_err(|_| ERR_PARSE)?;
                f(&mzml)
            }
        }
    }
}

fn file_is_broken(h: *const ParsedFile) -> bool {
    !h.is_null() && unsafe { (*h).broken }
}

fn mark_file_broken(h: *const ParsedFile) {
    if !h.is_null() {
        unsafe { (*(h as *mut ParsedFile)).broken = true };
    }
}

impl ScanSource for ParsedFile {
    fn for_each_summary(&mut self, cb: &mut dyn FnMut(usize, ScanSummary)) {
        match &mut self.source {
            FileSource::Full(mzml) => mzml.for_each_summary(cb),
            FileSource::Lazy(file) => file.for_each_summary(cb),
            FileSource::Remote(ion) => ion.for_each_summary(cb),
        }
    }

    fn load_scan(&mut self, index: usize, mz: &mut Vec<f64>, intensity: &mut Vec<f64>) -> bool {
        match &mut self.source {
            FileSource::Full(mzml) => mzml.load_scan(index, mz, intensity),
            FileSource::Lazy(file) => file.load_scan(index, mz, intensity),
            FileSource::Remote(ion) => ion.load_scan(index, mz, intensity),
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
        let _ = catch_unwind(AssertUnwindSafe(|| {
            drop(unsafe { Box::from_raw(handle) });
        }));
    }
}

/// Parse mzML bytes and store the result in `dest`.
///
/// # Safety
/// `data_ptr` must point to `data_len` readable bytes.
/// `dest` must be a valid writable pointer to `*mut ParsedFile`.
/// On success, `*dest` must be freed with `free_mzml`.
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
        unsafe { *dest = Box::into_raw(Box::new(ParsedFile::new(FileSource::Full(Box::new(parsed))))) };
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
        let ion = IonReader::open_source(
            Arc::new(CallbackSource::new(move |range| read_range(&reader, range)))
                as Arc<dyn ReadBytes>,
            ReadOptions {
                max_cached_bytes: cache_bytes,
                ..Default::default()
            },
        )
        .map_err(|_| ERR_PARSE)?;

        unsafe { *dest = Box::into_raw(Box::new(ParsedFile::new(FileSource::Remote(Box::new(ion))))) };
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
        let owned = IonReader::open_source(
            Arc::new(BytesSource::new(arc)) as Arc<dyn ReadBytes>,
            ReadOptions {
                max_cached_bytes: max_cache_size,
                ..Default::default()
            },
        )
        .map_err(|_| ERR_PARSE)?;

        unsafe { *dest = Box::into_raw(Box::new(ParsedFile::new(FileSource::Lazy(Box::new(owned))))) };
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
        let opened_file = IonReader::open_file(
            file_path,
            ReadOptions {
                max_cached_bytes: max_cache_size,
                ..Default::default()
            },
        )
        .map_err(|_| ERR_PARSE)?;
        unsafe { *dest = Box::into_raw(Box::new(ParsedFile::new(FileSource::Lazy(Box::new(opened_file))))) };
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => ERR_PANIC,
    }
}

const LARGEST_OPEN_GAP: u64 = 131072;

fn open_gap_for(ranges: &[ByteRange]) -> u64 {
    let total = ranges
        .iter()
        .map(|range| range.offset + range.length)
        .max()
        .unwrap_or(0);
    LARGEST_OPEN_GAP.min(total / 8)
}

/// Plan the byte ranges an open needs.
///
/// # Safety
/// `header_ptr` must point to at least `header_len` readable bytes, and `header_len` must be at
/// least 1024 so the header can be parsed.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn plan_open(
    header_ptr: *const u8,
    header_len: usize,
    out: *mut Buf,
) -> c_int {
    if header_ptr.is_null() || out.is_null() {
        return ERR_INVALID_ARGS;
    }
    clear_buf(out);

    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let header = unsafe { slice::from_raw_parts(header_ptr, header_len) };
        let mut ranges = open_ranges(header).map_err(|_| ERR_FAST_PATH)?;
        let gap = open_gap_for(&ranges);
        coalesce_byte_ranges(&mut ranges, gap);
        let bytes = pack_byte_ranges(&ranges);
        write_buf(out, bytes);
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
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
    if file_is_broken(h) {
        return ERR_PANIC;
    }
    clear_buf(out);
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
        Err(_) => {
            mark_file_broken(h);
            ERR_PANIC
        }
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
    if file_is_broken(h) {
        return ERR_PANIC;
    }
    clear_buf(out);
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let file = unsafe { &mut *h };
        file.with_mzml(|mzml| {
            let mut bytes = Vec::new();
            write_mzml_to_ion(
                mzml,
                WriteOptions {
                    compression_level,
                    force_f32: f32_compress != 0,
                    ..Default::default()
                },
                &mut bytes,
            )
            .map_err(|_| ERR_ENCODE)?;
            write_buf(out, bytes.into_boxed_slice());
            Ok(())
        })?;
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => {
            mark_file_broken(h);
            ERR_PANIC
        }
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
        let mut output = FileWriter::open(output_path).map_err(|_| ERR_ENCODE)?;

        let config = WriteOptions {
            compression_level,
            force_f32: force_f32 != 0,
            block_size: ION_BLOCK_SIZE_BYTES,
            parallel: true,
            section_storage: SectionStorage::Disk,
            mz_window: DEFAULT_MZ_WINDOW,
        };

        {
            let mut writer = IonWriter::create(&mut output, config).map_err(|_| ERR_ENCODE)?;
            writer.write_stream(&mut reader).map_err(|_| ERR_ENCODE)?;
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
    clear_buf(out);
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let data = DataXY {
            x: unsafe { slice::from_raw_parts(x_ptr, len) }.to_vec(),
            y: unsafe { slice::from_raw_parts(y_ptr, len) }.to_vec(),
        };
        let peak = get_peak_rs(
            &data,
            &Roi::peak(rt, range),
            Some(build_peak_options(options)),
        );
        let bytes = build_peaks_bridge(&[peak]).ok_or(ERR_ENCODE)?;
        write_buf(out, bytes);
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => ERR_PANIC,
    }
}

fn peak_shape_from_code(shape: c_int) -> PeakShape {
    if shape == 0 {
        PeakShape::Gaussian
    } else {
        PeakShape::EMG
    }
}

/// Fit a Gaussian or EMG model to a peak and write the parameters JSON to `out`.
///
/// `shape` is 0 for Gaussian, 1 for EMG.
///
/// # Safety
/// `x_ptr` and `y_ptr` must point to `len` readable `f64` values.
/// `out` must be a valid writable `Buf` pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn fit_peak(
    x_ptr: *const f64,
    y_ptr: *const f64,
    len: usize,
    rt: f64,
    intensity: f64,
    shape: c_int,
    out: *mut Buf,
) -> c_int {
    if x_ptr.is_null() || y_ptr.is_null() || out.is_null() || len < 5 {
        return ERR_INVALID_ARGS;
    }
    clear_buf(out);
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let data = DataXY {
            x: unsafe { slice::from_raw_parts(x_ptr, len) }.to_vec(),
            y: unsafe { slice::from_raw_parts(y_ptr, len) }.to_vec(),
        };
        let params = fit_peak_rs(
            &data,
            &PeakSeed { rt, intensity },
            peak_shape_from_code(shape),
        );
        let bytes = build_fit_bridge(&params).ok_or(ERR_ENCODE)?;
        write_buf(out, bytes);
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => ERR_PANIC,
    }
}

/// Render a fitted peak model over `x` and write the y curve to `out_y`.
///
/// `shape` is 0 for Gaussian, 1 for EMG. The output has the same length as `x`.
///
/// # Safety
/// `x_ptr` must point to `len` readable `f64` values.
/// `out_y` must be a valid writable `Buf` pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn draw_peak(
    x_ptr: *const f64,
    len: usize,
    shape: c_int,
    height: f64,
    center: f64,
    fwhm: f64,
    tail: f64,
    out_y: *mut Buf,
) -> c_int {
    if x_ptr.is_null() || out_y.is_null() || len == 0 {
        return ERR_INVALID_ARGS;
    }
    clear_buf(out_y);
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let x = unsafe { slice::from_raw_parts(x_ptr, len) }.to_vec();
        let params = PeakParameters {
            shape: peak_shape_from_code(shape),
            height,
            center,
            fwhm,
            tail,
            r2: 0.0,
        };
        let data = DataXY {
            x: x.clone(),
            y: vec![0.0; len],
        };
        let drawn = draw_peak_rs(&data, &params);
        write_buf(out_y, f64_to_u8(&drawn.y));
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
    if file_is_broken(h) {
        return ERR_PANIC;
    }
    clear_buf(out);
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
        let mut reader = match file.source_mut() {
            FileSource::Full(mzml) => EicReader::Mzml(mzml.as_mut()),
            FileSource::Lazy(ion) => EicReader::Ion(ion.as_mut()),
            FileSource::Remote(ion) => EicReader::Ion(ion.as_mut()),
        };

        let peaks = get_peaks_from_eic_rs(
            &mut reader,
            FromTo { from, to },
            &items,
            Some(build_peak_options(opts)),
            cores,
        )
        .map_err(fast_error_to_code)?;
        let bytes = build_eic_peaks_bridge(&peaks).ok_or(ERR_ENCODE)?;
        write_buf(out, bytes);
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => {
            mark_file_broken(h);
            ERR_PANIC
        }
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
    ranges: *const f64,
    n: usize,
    opts: *const CPeakOptions,
    cores: usize,
    out: *mut Buf,
) -> i32 {
    if h.is_null()
        || sample_indices.is_null()
        || target_rts.is_null()
        || ranges.is_null()
        || out.is_null()
        || n == 0
    {
        return ERR_INVALID_ARGS;
    }
    if file_is_broken(h) {
        return ERR_PANIC;
    }
    clear_buf(out);
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), i32> {
        let file = unsafe { &mut *h };
        file.with_mzml(|mzml| {
            let sample_indices = unsafe { slice::from_raw_parts(sample_indices, n) };
            let target_rts = unsafe { slice::from_raw_parts(target_rts, n) };
            let ranges = unsafe { slice::from_raw_parts(ranges, n) };
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
                        return Roi::chrom("", usize::MAX, 0.0, 0.0);
                    }
                    let sample_index = raw_index as usize;
                    if sample_index >= chroms.len() {
                        return Roi::chrom("", sample_index, 0.0, 0.0);
                    }
                    Roi::chrom(
                        chroms[sample_index].id.clone(),
                        sample_index,
                        target_rts[i],
                        ranges[i],
                    )
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
            let bytes = build_chrom_peaks_bridge(&out_arr).ok_or(ERR_ENCODE)?;
            write_buf(out, bytes);
            Ok(())
        })?;
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => {
            mark_file_broken(h);
            ERR_PANIC
        }
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
    clear_buf(out);
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
        let bytes = build_peaks_bridge(&peaks).ok_or(ERR_ENCODE)?;
        write_buf(out, bytes);
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => ERR_PANIC,
    }
}

/// Write the noise window width and the noise intensity for an intensity array.
///
/// # Safety
/// `y_ptr` must point to `len` readable `f32` values.
/// `out_width` and `out_intensity` must be valid writable pointers.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn find_noise_level(
    y_ptr: *const f32,
    len: usize,
    out_width: *mut usize,
    out_intensity: *mut f64,
) -> c_int {
    if y_ptr.is_null() || out_width.is_null() || out_intensity.is_null() || len == 0 {
        return ERR_INVALID_ARGS;
    }
    match catch_unwind(AssertUnwindSafe(|| {
        find_noise_level_rs(unsafe { slice::from_raw_parts(y_ptr, len) })
    })) {
        Ok(noise) => {
            unsafe {
                *out_width = noise.width;
                *out_intensity = noise.intensity;
            }
            OK
        }
        Err(_) => ERR_PANIC,
    }
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
    if file_is_broken(h) {
        return ERR_PANIC;
    }
    clear_buf(out_x);
    clear_buf(out_y);
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let file = unsafe { &mut *h };

        let mut reader = match file.source_mut() {
            FileSource::Full(mzml) => EicReader::Mzml(mzml.as_mut()),
            FileSource::Lazy(ion) => EicReader::Ion(ion.as_mut()),
            FileSource::Remote(ion) => EicReader::Ion(ion.as_mut()),
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

        let x_bytes = f64_to_u8(&eic.x);
        let y_bytes = f64_to_u8(&eic.y);
        write_buf(out_x, x_bytes);
        write_buf(out_y, y_bytes);
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => {
            mark_file_broken(h);
            ERR_PANIC
        }
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
    if file_is_broken(h) {
        return ERR_PANIC;
    }
    clear_buf(out);

    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let file = unsafe { &mut *h };

        let ranges = match file.source_mut() {
            FileSource::Lazy(ion) => plan_eic_ranges(ion.as_mut(), target, from, to, ppm, mz_tol),
            FileSource::Remote(ion) => plan_eic_ranges(ion.as_mut(), target, from, to, ppm, mz_tol),
            FileSource::Full(_) => Err(FastError::UnsupportedBackend),
        }
        .map_err(fast_error_to_code)?;

        let bytes = pack_byte_ranges(&ranges);
        write_buf(out, bytes);
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => {
            mark_file_broken(h);
            ERR_PANIC
        }
    }
}

fn scan_query_from_parts(query_type: u8, from_value: f64, to_value: f64) -> Option<ScanQuery> {
    match query_type {
        0 => {
            if !from_value.is_finite() || !to_value.is_finite() {
                return None;
            }
            Some(ScanQuery::RtRange(FromTo {
                from: from_value,
                to: to_value,
            }))
        }
        1 => {
            if !from_value.is_finite() {
                return None;
            }
            Some(ScanQuery::ClosestRt(from_value))
        }
        2 => {
            if !from_value.is_finite() || !to_value.is_finite() {
                return None;
            }
            Some(ScanQuery::MzRange(FromTo {
                from: from_value,
                to: to_value,
            }))
        }
        3 => {
            if !from_value.is_finite() || from_value <= 0.0 {
                return None;
            }
            Some(ScanQuery::ClosestMz(from_value))
        }
        _ => None,
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn plan_scans(
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
    if file_is_broken(h) {
        return ERR_PANIC;
    }
    clear_buf(out);
    let query = match scan_query_from_parts(query_type, from_value, to_value) {
        Some(query) => query,
        None => return ERR_INVALID_ARGS,
    };

    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let file = unsafe { &mut *h };
        let ranges = match file.source_mut() {
            FileSource::Lazy(ion) => {
                plan_scan_ranges(ion.as_mut(), query, TimeUnit::Minutes, level)
            }
            FileSource::Remote(ion) => {
                plan_scan_ranges(ion.as_mut(), query, TimeUnit::Minutes, level)
            }
            FileSource::Full(_) => Err(FastError::UnsupportedBackend),
        }
        .map_err(fast_error_to_code)?;

        let bytes = pack_byte_ranges(&ranges);
        write_buf(out, bytes);
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => {
            mark_file_broken(h);
            ERR_PANIC
        }
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn plan_image(
    h: *mut ParsedFile,
    target: f64,
    tolerance: f64,
    level: u8,
    out: *mut Buf,
) -> c_int {
    if h.is_null() || out.is_null() {
        return ERR_INVALID_ARGS;
    }
    if file_is_broken(h) {
        return ERR_PANIC;
    }
    clear_buf(out);

    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let file = unsafe { &mut *h };
        let ranges = match file.source_mut() {
            FileSource::Lazy(ion) => crate::utilities::ion_image::plan_ion_image_ranges(
                ion.as_mut(),
                target,
                tolerance,
                level,
            ),
            FileSource::Remote(ion) => crate::utilities::ion_image::plan_ion_image_ranges(
                ion.as_mut(),
                target,
                tolerance,
                level,
            ),
            FileSource::Full(_) => Err(FastError::UnsupportedBackend),
        }
        .map_err(fast_error_to_code)?;

        let bytes = pack_byte_ranges(&ranges);
        write_buf(out, bytes);
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => {
            mark_file_broken(h);
            ERR_PANIC
        }
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
    if file_is_broken(h) {
        return ERR_PANIC;
    }
    let query = match scan_query_from_parts(query_type, from_value, to_value) {
        Some(query) => query,
        None => return ERR_INVALID_ARGS,
    };
    clear_buf(out);
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let (_, scans) = get_scans_rs(unsafe { &mut *h }, query, TimeUnit::Minutes, level);
        let bytes = build_scans_bridge(scans).ok_or(ERR_ENCODE)?;
        write_buf(out, bytes);
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => {
            mark_file_broken(h);
            ERR_PANIC
        }
    }
}

/// Build a 2D ion image for a target m/z by summing intensity in `[target - tolerance, target + tolerance]`
/// per spectrum and scattering the mean into a position_x/position_y grid. Writes an ion image bridge to `out`.
///
/// # Safety
/// `h` must be a valid `ParsedFile` pointer from this library.
/// `out` must be a valid writable `Buf` pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn get_ion_image(
    h: *mut ParsedFile,
    target: f64,
    tolerance: f64,
    level: u8,
    out: *mut Buf,
) -> c_int {
    if h.is_null() || out.is_null() {
        return ERR_INVALID_ARGS;
    }
    if file_is_broken(h) {
        return ERR_PANIC;
    }
    clear_buf(out);
    if !target.is_finite() || !tolerance.is_finite() || tolerance < 0.0 {
        return ERR_INVALID_ARGS;
    }
    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let file = unsafe { &mut *h };

        let mut reader = match file.source_mut() {
            FileSource::Full(mzml) => EicReader::Mzml(mzml.as_mut()),
            FileSource::Lazy(ion) => EicReader::Ion(ion.as_mut()),
            FileSource::Remote(ion) => EicReader::Ion(ion.as_mut()),
        };

        let image =
            crate::utilities::ion_image::compute_ion_image(&mut reader, target, tolerance, level);
        let bytes = build_image_bridge(&image).ok_or(ERR_ENCODE)?;
        write_buf(out, bytes);
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => {
            mark_file_broken(h);
            ERR_PANIC
        }
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn image_begin(
    h: *mut ParsedFile,
    target: f64,
    tolerance: f64,
    level: u8,
    out_session: *mut usize,
) -> c_int {
    if h.is_null() || out_session.is_null() {
        return ERR_INVALID_ARGS;
    }
    if file_is_broken(h) {
        return ERR_PANIC;
    }
    if !target.is_finite() || !tolerance.is_finite() || tolerance < 0.0 {
        return ERR_INVALID_ARGS;
    }

    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let file = unsafe { &mut *h };
        let mut reader = match file.source_mut() {
            FileSource::Full(mzml) => EicReader::Mzml(mzml.as_mut()),
            FileSource::Lazy(ion) => EicReader::Ion(ion.as_mut()),
            FileSource::Remote(ion) => EicReader::Ion(ion.as_mut()),
        };
        let session =
            crate::utilities::ion_image::image_session_begin(&mut reader, target, tolerance, level);
        let pointer = Box::into_raw(Box::new(session)) as usize;
        unsafe { *out_session = pointer };
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => {
            mark_file_broken(h);
            ERR_PANIC
        }
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn image_scan_count(
    session: *mut crate::utilities::ion_image::ImageSession,
    out_count: *mut usize,
) -> c_int {
    if session.is_null() || out_count.is_null() {
        return ERR_INVALID_ARGS;
    }
    match catch_unwind(AssertUnwindSafe(|| {
        let session = unsafe { &*session };
        unsafe { *out_count = session.scan_count() };
    })) {
        Ok(()) => OK,
        Err(_) => ERR_PANIC,
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn image_ranges(
    h: *mut ParsedFile,
    session: *mut crate::utilities::ion_image::ImageSession,
    from: usize,
    count: usize,
    out: *mut Buf,
) -> c_int {
    if h.is_null() || session.is_null() || out.is_null() {
        return ERR_INVALID_ARGS;
    }
    if session_is_broken(session) {
        return ERR_PANIC;
    }
    if file_is_broken(h) {
        return ERR_PANIC;
    }
    clear_buf(out);

    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let file = unsafe { &mut *h };
        let session = unsafe { &*session };
        let ranges = match file.source_mut() {
            FileSource::Lazy(ion) => session.ranges(ion.as_mut(), from, count),
            FileSource::Remote(ion) => session.ranges(ion.as_mut(), from, count),
            FileSource::Full(_) => Err(FastError::UnsupportedBackend),
        }
        .map_err(fast_error_to_code)?;

        let bytes = pack_byte_ranges(&ranges);
        write_buf(out, bytes);
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => {
            mark_file_broken(h);
            ERR_PANIC
        }
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn image_fold(
    h: *mut ParsedFile,
    session: *mut crate::utilities::ion_image::ImageSession,
    from: usize,
    count: usize,
) -> c_int {
    if h.is_null() || session.is_null() {
        return ERR_INVALID_ARGS;
    }
    if session_is_broken(session) {
        return ERR_PANIC;
    }
    if file_is_broken(h) {
        return ERR_PANIC;
    }

    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let file = unsafe { &mut *h };
        let session = unsafe { &mut *session };
        let mut reader = match file.source_mut() {
            FileSource::Full(mzml) => EicReader::Mzml(mzml.as_mut()),
            FileSource::Lazy(ion) => EicReader::Ion(ion.as_mut()),
            FileSource::Remote(ion) => EicReader::Ion(ion.as_mut()),
        };
        session.fold(&mut reader, from, count);
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => {
            mark_file_broken(h);
            ERR_PANIC
        }
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn image_finish(
    session: *mut crate::utilities::ion_image::ImageSession,
    out: *mut Buf,
) -> c_int {
    if session.is_null() || out.is_null() {
        return ERR_INVALID_ARGS;
    }
    if session_is_broken(session) {
        return ERR_PANIC;
    }
    clear_buf(out);

    match catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let session = unsafe { &*session };
        let image = session.finish();
        let bytes = build_image_bridge(&image).ok_or(ERR_ENCODE)?;
        write_buf(out, bytes);
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => {
            mark_session_broken(session);
            ERR_PANIC
        }
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn image_free(session: *mut crate::utilities::ion_image::ImageSession) {
    if !session.is_null() {
        let _ = catch_unwind(AssertUnwindSafe(|| {
            drop(unsafe { Box::from_raw(session) });
        }));
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
    out: *mut Buf,
) -> c_int {
    if dir.is_null() || out.is_null() {
        return ERR_INVALID_ARGS;
    }
    clear_buf(out);
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
            mz_estimator: MzEstimatorKind::MedianMzApex,
            ..AlignmentOptions::default()
        };

        let feats = get_features_rs(
            path,
            FromTo { from, to },
            feature_opts,
            alignment_opts,
            if cores > 0 { cores as usize } else { 1 },
        )
        .map_err(|_| ERR_FAST_PATH)?;

        let bytes = build_consensus_bridge(&feats).ok_or(ERR_ENCODE)?;
        write_buf(out, bytes);
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
    out: *mut Buf,
) -> c_int {
    if h.is_null() || out.is_null() || !from.is_finite() || !to.is_finite() || to <= from {
        return ERR_INVALID_ARGS;
    }
    if file_is_broken(h) {
        return ERR_PANIC;
    }
    clear_buf(out);
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

        let mut reader = match file.source_mut() {
            FileSource::Full(mzml) => EicReader::Mzml(mzml.as_mut()),
            FileSource::Lazy(ion) => EicReader::Ion(ion.as_mut()),
            FileSource::Remote(ion) => EicReader::Ion(ion.as_mut()),
        };

        let feats = find_features_rs(
            &mut reader,
            FromTo { from, to },
            Some(opts),
            if cores > 0 { cores as usize } else { 1 },
        )
        .map_err(|_| ERR_FAST_PATH)?;

        let bytes = build_features_bridge(&feats).ok_or(ERR_ENCODE)?;
        write_buf(out, bytes);
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => {
            mark_file_broken(h);
            ERR_PANIC
        }
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
    clear_buf(out);
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
    if file_is_broken(h) {
        return ERR_PANIC;
    }
    clear_buf(out);
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
        let refs: Vec<&Roi> = rois.iter().collect();
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
                    id: rois[i].id(),
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
        let bytes = build_found_features_bridge(&arr).ok_or(ERR_ENCODE)?;
        write_buf(out, bytes);
        Ok(())
    })) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => {
            mark_file_broken(h);
            ERR_PANIC
        }
    }
}

fn build_scans_bridge(
    scans: Vec<crate::utilities::calculate_eic::CentroidScan>,
) -> Option<Box<[u8]>> {
    let scan_count = scans.len() as u64;
    let mut total_points: u64 = 0;
    for scan in scans.iter() {
        total_points = total_points.checked_add(scan.mz.len() as u64)?;
    }

    let mut builder = BridgeBuilder::new(QUANTION_PAYLOAD_SCANS, scan_count);
    builder.add_section(
        QUANTION_SECTION_POINT_STARTS,
        QUANTION_ELEMENT_U64,
        scan_count + 1,
    );
    builder.add_section(QUANTION_SECTION_MZ, QUANTION_ELEMENT_F64, total_points);
    builder.add_section(QUANTION_SECTION_INTENSITY, QUANTION_ELEMENT_F64, total_points);
    for id in SCAN_METADATA_SECTIONS {
        builder.add_section(id, QUANTION_ELEMENT_F64, scan_count);
    }
    let mut bridge = builder.build()?;

    let mut point_starts = Vec::with_capacity(scans.len() + 1);
    point_starts.push(0u64);
    let mut next_point: u64 = 0;

    for (index, scan) in scans.into_iter().enumerate() {
        let position = index as u64;
        bridge.write_f64_at(QUANTION_SECTION_MZ, next_point, &scan.mz);
        bridge.write_f64_at(QUANTION_SECTION_INTENSITY, next_point, &scan.intensity);

        let summary = &scan.metadata;
        bridge.write_f64_at(QUANTION_SECTION_RT, position, &[scan.rt]);
        bridge.write_f64_at(QUANTION_SECTION_RT_SECONDS, position, &[summary.rt_seconds]);
        bridge.write_f64_at(
            QUANTION_SECTION_BASE_PEAK_MZ,
            position,
            &[summary.base_peak_mz],
        );
        bridge.write_f64_at(
            QUANTION_SECTION_SELECTED_ION_MZ,
            position,
            &[summary.selected_ion_mz],
        );
        bridge.write_f64_at(
            QUANTION_SECTION_BASE_PEAK_INT,
            position,
            &[summary.base_peak_int],
        );
        bridge.write_f64_at(
            QUANTION_SECTION_TOTAL_ION_CURRENT,
            position,
            &[summary.total_ion_current],
        );
        bridge.write_f64_at(
            QUANTION_SECTION_MS_LEVEL,
            position,
            &[f64::from(summary.ms_level)],
        );
        bridge.write_f64_at(
            QUANTION_SECTION_POLARITY,
            position,
            &[f64::from(summary.polarity)],
        );
        bridge.write_f64_at(
            QUANTION_SECTION_POSITION_X,
            position,
            &[f64::from(summary.position_x)],
        );
        bridge.write_f64_at(
            QUANTION_SECTION_POSITION_Y,
            position,
            &[f64::from(summary.position_y)],
        );
        bridge.write_f64_at(
            QUANTION_SECTION_POSITION_Z,
            position,
            &[f64::from(summary.position_z)],
        );

        next_point += scan.mz.len() as u64;
        point_starts.push(next_point);
    }

    bridge.write_u64_section(QUANTION_SECTION_POINT_STARTS, &point_starts);
    Some(bridge.into_bytes())
}

fn build_image_bridge(image: &crate::utilities::ion_image::IonImage) -> Option<Box<[u8]>> {
    let cell_count = image.data.len() as u64;
    let mut builder = BridgeBuilder::new(QUANTION_PAYLOAD_ION_IMAGE, cell_count);
    builder.add_section(QUANTION_SECTION_IMAGE_SHAPE, QUANTION_ELEMENT_U32, 6);
    builder.add_section(QUANTION_SECTION_IMAGE_DATA, QUANTION_ELEMENT_F64, cell_count);
    builder.add_section(
        QUANTION_SECTION_IMAGE_COUNTS,
        QUANTION_ELEMENT_U32,
        image.counts.len() as u64,
    );
    let mut bridge = builder.build()?;

    let shape = [
        image.width,
        image.height,
        image.min_x,
        image.min_y,
        image.min_z,
        image.max_z,
    ];
    bridge.write_u32_section(QUANTION_SECTION_IMAGE_SHAPE, &shape);
    bridge.write_f64_section(QUANTION_SECTION_IMAGE_DATA, &image.data);
    bridge.write_u32_section(QUANTION_SECTION_IMAGE_COUNTS, &image.counts);
    Some(bridge.into_bytes())
}

const DEFAULT_EIC_HALF_WIDTH: f64 = 0.5;

fn build_eic_rois(
    rts: &[f64],
    mzs: &[f64],
    wins: &[f64],
    ids_off: *const u32,
    ids_len: *const u32,
    ids_buf: *const u8,
    ids_buf_len: usize,
) -> Vec<Roi> {
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
            let (rt, mz, window) = (rts[i], mzs[i], wins[i]);
            let has_target = rt.is_finite() && mz.is_finite();
            let range = if window.is_finite() && window > 0.0 {
                window
            } else {
                DEFAULT_EIC_HALF_WIDTH
            };
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

            if has_target {
                Roi::eic(id, rt, mz, range)
            } else {
                Roi::eic(id, 0.0, 0.0, 0.0)
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


fn count_as_u32(value: usize) -> u32 {
    u32::try_from(value).unwrap_or(u32::MAX)
}

fn build_peaks_bridge(peaks: &[crate::utilities::structs::Peak]) -> Option<Box<[u8]>> {
    let from: Vec<f64> = peaks.iter().map(|peak| peak.from).collect();
    let to: Vec<f64> = peaks.iter().map(|peak| peak.to).collect();
    let rt: Vec<f64> = peaks.iter().map(|peak| peak.rt).collect();
    let integral: Vec<f64> = peaks.iter().map(|peak| peak.integral).collect();
    let intensity: Vec<f64> = peaks.iter().map(|peak| peak.intensity).collect();
    let point_count: Vec<u32> = peaks.iter().map(|peak| count_as_u32(peak.n_points)).collect();
    let noise: Vec<f64> = peaks.iter().map(|peak| peak.noise).collect();
    let r2: Vec<f64> = peaks.iter().map(|peak| peak.r2.unwrap_or(f64::NAN)).collect();

    build_record_bridge(
        QUANTION_PAYLOAD_PEAKS,
        peaks.len() as u64,
        &[
            Column::Numbers { id: QUANTION_SECTION_PEAK_FROM, values: &from },
            Column::Numbers { id: QUANTION_SECTION_PEAK_TO, values: &to },
            Column::Numbers { id: QUANTION_SECTION_PEAK_RT, values: &rt },
            Column::Numbers { id: QUANTION_SECTION_PEAK_INTEGRAL, values: &integral },
            Column::Numbers { id: QUANTION_SECTION_PEAK_INTENSITY, values: &intensity },
            Column::Counts { id: QUANTION_SECTION_PEAK_POINT_COUNT, values: &point_count },
            Column::Numbers { id: QUANTION_SECTION_PEAK_NOISE, values: &noise },
            Column::Numbers { id: QUANTION_SECTION_PEAK_R2, values: &r2 },
        ],
    )
}

fn build_fit_bridge(params: &Option<PeakParameters>) -> Option<Box<[u8]>> {
    let (shape, height, center, fwhm, tail, r2) = match params {
        Some(found) => (
            vec![match found.shape {
                PeakShape::Gaussian => 0.0,
                PeakShape::EMG => 1.0,
            }],
            vec![found.height],
            vec![found.center],
            vec![found.fwhm],
            vec![found.tail],
            vec![found.r2],
        ),
        None => (vec![], vec![], vec![], vec![], vec![], vec![]),
    };

    build_record_bridge(
        QUANTION_PAYLOAD_FIT_RESULT,
        shape.len() as u64,
        &[
            Column::Numbers { id: QUANTION_SECTION_FIT_SHAPE, values: &shape },
            Column::Numbers { id: QUANTION_SECTION_FIT_HEIGHT, values: &height },
            Column::Numbers { id: QUANTION_SECTION_FIT_CENTER, values: &center },
            Column::Numbers { id: QUANTION_SECTION_FIT_FWHM, values: &fwhm },
            Column::Numbers { id: QUANTION_SECTION_FIT_TAIL, values: &tail },
            Column::Numbers { id: QUANTION_SECTION_FIT_R2, values: &r2 },
        ],
    )
}

fn build_eic_peaks_bridge(rows: &[(&str, f64, f64, crate::utilities::structs::Peak)]) -> Option<Box<[u8]>> {
    let ids: Vec<&str> = rows.iter().map(|row| row.0).collect();
    let ort: Vec<f64> = rows.iter().map(|row| row.1).collect();
    let mz: Vec<f64> = rows.iter().map(|row| row.2).collect();
    let rt: Vec<f64> = rows.iter().map(|row| row.3.rt).collect();
    let from: Vec<f64> = rows.iter().map(|row| row.3.from).collect();
    let to: Vec<f64> = rows.iter().map(|row| row.3.to).collect();
    let intensity: Vec<f64> = rows.iter().map(|row| row.3.intensity).collect();
    let integral: Vec<f64> = rows.iter().map(|row| row.3.integral).collect();
    let noise: Vec<f64> = rows.iter().map(|row| row.3.noise).collect();

    build_record_bridge(
        QUANTION_PAYLOAD_EIC_PEAKS,
        rows.len() as u64,
        &[
            Column::Text {
                starts_id: QUANTION_SECTION_EIC_PEAK_ID_STARTS,
                bytes_id: QUANTION_SECTION_EIC_PEAK_ID_BYTES,
                values: &ids,
            },
            Column::Numbers { id: QUANTION_SECTION_EIC_PEAK_MZ, values: &mz },
            Column::Numbers { id: QUANTION_SECTION_EIC_PEAK_ORT, values: &ort },
            Column::Numbers { id: QUANTION_SECTION_EIC_PEAK_RT, values: &rt },
            Column::Numbers { id: QUANTION_SECTION_EIC_PEAK_FROM, values: &from },
            Column::Numbers { id: QUANTION_SECTION_EIC_PEAK_TO, values: &to },
            Column::Numbers { id: QUANTION_SECTION_EIC_PEAK_INTENSITY, values: &intensity },
            Column::Numbers { id: QUANTION_SECTION_EIC_PEAK_INTEGRAL, values: &integral },
            Column::Numbers { id: QUANTION_SECTION_EIC_PEAK_NOISE, values: &noise },
        ],
    )
}

fn build_features_bridge(features: &[crate::utilities::find_features::Feature]) -> Option<Box<[u8]>> {
    let mz: Vec<f64> = features.iter().map(|item| item.mz).collect();
    let rt: Vec<f64> = features.iter().map(|item| item.rt).collect();
    let from: Vec<f64> = features.iter().map(|item| item.from).collect();
    let to: Vec<f64> = features.iter().map(|item| item.to).collect();
    let intensity: Vec<f64> = features.iter().map(|item| item.intensity).collect();
    let integral: Vec<f64> = features.iter().map(|item| item.integral).collect();
    let point_count: Vec<u32> = features.iter().map(|item| count_as_u32(item.n_points)).collect();
    let noise: Vec<f64> = features.iter().map(|item| item.noise).collect();

    build_record_bridge(
        QUANTION_PAYLOAD_FEATURES,
        features.len() as u64,
        &[
            Column::Numbers { id: QUANTION_SECTION_FEATURE_MZ, values: &mz },
            Column::Numbers { id: QUANTION_SECTION_FEATURE_RT, values: &rt },
            Column::Numbers { id: QUANTION_SECTION_FEATURE_FROM, values: &from },
            Column::Numbers { id: QUANTION_SECTION_FEATURE_TO, values: &to },
            Column::Numbers { id: QUANTION_SECTION_FEATURE_INTENSITY, values: &intensity },
            Column::Numbers { id: QUANTION_SECTION_FEATURE_INTEGRAL, values: &integral },
            Column::Counts { id: QUANTION_SECTION_FEATURE_POINT_COUNT, values: &point_count },
            Column::Numbers { id: QUANTION_SECTION_FEATURE_NOISE, values: &noise },
        ],
    )
}


#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
fn build_consensus_bridge(
    features: &[crate::utilities::get_features::ConsensusFeature],
) -> Option<Box<[u8]>> {
    let mz: Vec<f64> = features.iter().map(|item| item.mz).collect();
    let rt: Vec<f64> = features.iter().map(|item| item.rt).collect();
    let from: Vec<f64> = features.iter().map(|item| item.from).collect();
    let to: Vec<f64> = features.iter().map(|item| item.to).collect();
    let intensity: Vec<f64> = features.iter().map(|item| item.intensity).collect();
    let integral: Vec<f64> = features.iter().map(|item| item.integral).collect();
    let frequency: Vec<f64> = features.iter().map(|item| item.frequency).collect();

    build_record_bridge(
        QUANTION_PAYLOAD_CONSENSUS_FEATURES,
        features.len() as u64,
        &[
            Column::Numbers { id: QUANTION_SECTION_CONSENSUS_MZ, values: &mz },
            Column::Numbers { id: QUANTION_SECTION_CONSENSUS_RT, values: &rt },
            Column::Numbers { id: QUANTION_SECTION_CONSENSUS_FROM, values: &from },
            Column::Numbers { id: QUANTION_SECTION_CONSENSUS_TO, values: &to },
            Column::Numbers { id: QUANTION_SECTION_CONSENSUS_INTENSITY, values: &intensity },
            Column::Numbers { id: QUANTION_SECTION_CONSENSUS_INTEGRAL, values: &integral },
            Column::Numbers { id: QUANTION_SECTION_CONSENSUS_FREQUENCY, values: &frequency },
        ],
    )
}

fn build_found_features_bridge(rows: &[FoundFeatureOut]) -> Option<Box<[u8]>> {
    let ids: Vec<&str> = rows.iter().map(|row| row.id).collect();
    let mz: Vec<f64> = rows.iter().map(|row| row.mz).collect();
    let rt: Vec<f64> = rows.iter().map(|row| row.rt).collect();
    let from: Vec<f64> = rows.iter().map(|row| row.from).collect();
    let to: Vec<f64> = rows.iter().map(|row| row.to).collect();
    let intensity: Vec<f64> = rows.iter().map(|row| row.intensity).collect();
    let integral: Vec<f64> = rows.iter().map(|row| row.integral).collect();
    let point_count: Vec<u32> = rows.iter().map(|row| count_as_u32(row.n_points)).collect();
    let noise: Vec<f64> = rows.iter().map(|row| row.noise).collect();

    build_record_bridge(
        QUANTION_PAYLOAD_FOUND_FEATURES,
        rows.len() as u64,
        &[
            Column::Text {
                starts_id: QUANTION_SECTION_FOUND_ID_STARTS,
                bytes_id: QUANTION_SECTION_FOUND_ID_BYTES,
                values: &ids,
            },
            Column::Numbers { id: QUANTION_SECTION_FOUND_MZ, values: &mz },
            Column::Numbers { id: QUANTION_SECTION_FOUND_RT, values: &rt },
            Column::Numbers { id: QUANTION_SECTION_FOUND_FROM, values: &from },
            Column::Numbers { id: QUANTION_SECTION_FOUND_TO, values: &to },
            Column::Numbers { id: QUANTION_SECTION_FOUND_INTENSITY, values: &intensity },
            Column::Numbers { id: QUANTION_SECTION_FOUND_INTEGRAL, values: &integral },
            Column::Counts { id: QUANTION_SECTION_FOUND_POINT_COUNT, values: &point_count },
            Column::Numbers { id: QUANTION_SECTION_FOUND_NOISE, values: &noise },
        ],
    )
}

fn build_chrom_peaks_bridge(rows: &[ChromPeakRowOut]) -> Option<Box<[u8]>> {
    let index: Vec<u32> = rows.iter().map(|row| count_as_u32(row.index)).collect();
    let ids: Vec<&str> = rows.iter().map(|row| row.id).collect();
    let timestamps: Vec<&str> = rows.iter().map(|row| row.timestamp).collect();
    let target_rt: Vec<f64> = rows.iter().map(|row| row.ort).collect();
    let rt: Vec<f64> = rows.iter().map(|row| row.rt).collect();
    let from: Vec<f64> = rows.iter().map(|row| row.from).collect();
    let to: Vec<f64> = rows.iter().map(|row| row.to).collect();
    let intensity: Vec<f64> = rows.iter().map(|row| row.intensity).collect();
    let integral: Vec<f64> = rows.iter().map(|row| row.integral).collect();
    let total_area: Vec<f64> = rows.iter().map(|row| row.total_area).collect();

    build_record_bridge(
        QUANTION_PAYLOAD_CHROM_PEAKS,
        rows.len() as u64,
        &[
            Column::Counts { id: QUANTION_SECTION_CHROM_INDEX, values: &index },
            Column::Numbers { id: QUANTION_SECTION_CHROM_TARGET_RT, values: &target_rt },
            Column::Numbers { id: QUANTION_SECTION_CHROM_RT, values: &rt },
            Column::Numbers { id: QUANTION_SECTION_CHROM_FROM, values: &from },
            Column::Numbers { id: QUANTION_SECTION_CHROM_TO, values: &to },
            Column::Numbers { id: QUANTION_SECTION_CHROM_INTENSITY, values: &intensity },
            Column::Numbers { id: QUANTION_SECTION_CHROM_INTEGRAL, values: &integral },
            Column::Numbers { id: QUANTION_SECTION_CHROM_TOTAL_AREA, values: &total_area },
            Column::Text {
                starts_id: QUANTION_SECTION_CHROM_ID_STARTS,
                bytes_id: QUANTION_SECTION_CHROM_ID_BYTES,
                values: &ids,
            },
            Column::Text {
                starts_id: QUANTION_SECTION_CHROM_TIMESTAMP_STARTS,
                bytes_id: QUANTION_SECTION_CHROM_TIMESTAMP_BYTES,
                values: &timestamps,
            },
        ],
    )
}

const SCAN_METADATA_SECTIONS: [u32; 11] = [
    QUANTION_SECTION_RT,
    QUANTION_SECTION_RT_SECONDS,
    QUANTION_SECTION_BASE_PEAK_MZ,
    QUANTION_SECTION_SELECTED_ION_MZ,
    QUANTION_SECTION_BASE_PEAK_INT,
    QUANTION_SECTION_TOTAL_ION_CURRENT,
    QUANTION_SECTION_MS_LEVEL,
    QUANTION_SECTION_POLARITY,
    QUANTION_SECTION_POSITION_X,
    QUANTION_SECTION_POSITION_Y,
    QUANTION_SECTION_POSITION_Z,
];

fn clear_buf(out: *mut Buf) {
    unsafe {
        ptr::write_unaligned(
            out,
            Buf {
                ptr: ptr::null_mut(),
                len: 0,
            },
        )
    };
}

fn session_is_broken(session: *const crate::utilities::ion_image::ImageSession) -> bool {
    !session.is_null() && unsafe { (*session).is_broken() }
}

fn mark_session_broken(session: *mut crate::utilities::ion_image::ImageSession) {
    if !session.is_null() {
        unsafe { (*session).mark_broken() };
    }
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
    let default_artifact = ArtifactFilter::default();
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
            noise_method: None,
            kernel_size: (o.kernel_size > 0).then_some(o.kernel_size as usize),
        }),
        baseline: Some(BaselineOptions {
            lambda: (o.lambda > 0).then_some(o.lambda as f64),
            max_iterations: (o.max_iterations > 0).then_some(o.max_iterations as usize),
            edge_slope_level: Some(1),
        }),
        artifact_filter: Some(ArtifactFilter {
            min_r2: if o.min_r2.is_finite() && o.min_r2 >= 0.0 {
                o.min_r2
            } else {
                default_artifact.min_r2
            },
            shape: peak_shape_from_code(o.shape),
        }),
    }
}

#[cfg(test)]
mod tests {
    use ionic::{
        mzml::structs::{
            BinaryDataArray, BinaryDataArrayList, CvParam, NumericArray, NumericType, Run, Scan,
            ScanList, Spectrum, SpectrumList,
        },
        write_mzml_to_ion,
    };

    use super::*;
    use crate::utilities::calculate_eic::plan_window_ranges;

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
        use ionic::ion::ByteRange;

        let fake = FakeReader {
            data: vec![1, 2, 3, 4, 5],
            call_count: std::sync::atomic::AtomicUsize::new(0),
            return_code: 0,
        };

        let range = ByteRange {
            offset: 100,
            length: 0,
        };
        let value = read_range(&fake, range).unwrap();

        assert_eq!(value.len(), 0);
        assert_eq!(fake.call_count.load(std::sync::atomic::Ordering::SeqCst), 0);
    }

    #[test]
    fn non_zero_read_calls_reader_once() {
        use ionic::ion::ByteRange;

        let fake = FakeReader {
            data: vec![10, 20, 30, 40, 50],
            call_count: std::sync::atomic::AtomicUsize::new(0),
            return_code: 0,
        };

        let range = ByteRange {
            offset: 0,
            length: 3,
        };
        let _result = read_range(&fake, range).unwrap();

        assert_eq!(fake.call_count.load(std::sync::atomic::Ordering::SeqCst), 1);
    }

    #[test]
    fn negative_return_code_becomes_error() {
        use ionic::ion::ByteRange;

        let fake = FakeReader {
            data: vec![10, 20, 30],
            call_count: std::sync::atomic::AtomicUsize::new(0),
            return_code: -1,
        };

        let range = ByteRange {
            offset: 0,
            length: 2,
        };
        let result = read_range(&fake, range);

        assert!(result.is_err());
    }

    #[test]
    fn len_exceeding_u32_max_returns_error_without_allocation() {
        use ionic::ion::ByteRange;

        let fake = FakeReader {
            data: vec![1, 2, 3],
            call_count: std::sync::atomic::AtomicUsize::new(0),
            return_code: 0,
        };

        let range = ByteRange {
            offset: 0,
            length: u64::from(u32::MAX) + 1,
        };
        let result = read_range(&fake, range);

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
            binary: Some(NumericArray::F64(values)),
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

    fn sample_mzml() -> MzML {
        MzML {
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
        }
    }

    fn build_ion(mz_window: f64, block_size: usize) -> Vec<u8> {
        let options = WriteOptions {
            compression_level: 0,
            force_f32: false,
            block_size,
            parallel: true,
            section_storage: SectionStorage::Memory,
            mz_window,
        };
        let mut bytes = Vec::new();
        write_mzml_to_ion(&sample_mzml(), options, &mut bytes).expect("ion encode should succeed");
        bytes
    }

    fn open_ion(bytes: Vec<u8>) -> *mut ParsedFile {
        let bytes_arc = Arc::from(bytes.into_boxed_slice());
        let ion = IonReader::open_source(
            Arc::new(ionic::ion::BytesSource::new(bytes_arc)) as Arc<dyn ionic::ion::ReadBytes>,
            ReadOptions::default(),
        )
        .expect("IonReader::open_bytes failed");
        Box::into_raw(Box::new(ParsedFile::new(FileSource::Lazy(Box::new(ion)))))
    }

    fn windowed_ion_bytes() -> Vec<u8> {
        build_ion(DEFAULT_MZ_WINDOW, ION_BLOCK_SIZE_BYTES)
    }

    fn open_windowed_ion() -> *mut ParsedFile {
        open_ion(windowed_ion_bytes())
    }

    static mut SERVED_BYTES: usize = 0;

    unsafe extern "C" fn read_from_slice(
        context: *mut c_void,
        offset: u64,
        length: u64,
        dest: *mut u8,
    ) -> i32 {
        let bytes = unsafe { &*(context as *const Vec<u8>) };
        let start = offset as usize;
        let end = start + length as usize;
        if end > bytes.len() {
            return -1;
        }
        unsafe {
            std::ptr::copy_nonoverlapping(bytes.as_ptr().add(start), dest, length as usize);
            SERVED_BYTES += length as usize;
        }
        0
    }

    #[test]
    fn parse_ion_source_reads_only_what_it_needs() {
        let fixture =
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/api/api.ion");
        let bytes = std::fs::read(&fixture).expect("read the api fixture");
        let total = bytes.len();
        unsafe { SERVED_BYTES = 0 };

        let mut handle: *mut ParsedFile = std::ptr::null_mut();
        let code = unsafe {
            parse_ion_source(
                Some(read_from_slice),
                &bytes as *const Vec<u8> as *mut c_void,
                0,
                &mut handle,
            )
        };
        assert_eq!(code, OK, "parse_ion_source failed");
        assert!(!handle.is_null());

        let served = unsafe { SERVED_BYTES };
        assert!(
            served < total,
            "opening pulled {served} of {total} bytes: it is not reading lazily"
        );
        println!("opening read {served} of {total} bytes");

        unsafe { free_mzml(handle) };
    }

    #[test]
    fn parse_ion_source_rejects_a_missing_callback() {
        let mut handle: *mut ParsedFile = std::ptr::null_mut();
        let code = unsafe { parse_ion_source(None, std::ptr::null_mut(), 0, &mut handle) };
        assert_eq!(code, ERR_INVALID_ARGS);
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
        let handle = open_windowed_ion();
        let out = new_buf_out();
        let code = unsafe { plan_eic(handle, -1.0, 0.0, 10.0, 20.0, 0.005, out) };
        assert_eq!(code, ERR_INVALID_ARGS);
        drop_buf_out(out);
        unsafe { free_mzml(handle) };
    }

    #[test]
    fn plan_eic_rejects_zero_tolerance() {
        let handle = open_windowed_ion();
        let out = new_buf_out();
        let code = unsafe { plan_eic(handle, 500.0, 0.0, 10.0, 0.0, 0.0, out) };
        assert_eq!(code, ERR_INVALID_ARGS);
        drop_buf_out(out);
        unsafe { free_mzml(handle) };
    }

    #[test]
    fn plan_eic_returns_ranges_for_windowed_ion() {
        let handle = open_windowed_ion();
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
    fn plan_eic_works_for_wide_window_large_block_ion() {
        let handle = open_ion(build_ion(250.0, 8 * 1024 * 1024));
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

    struct RecordingReader {
        data: Vec<u8>,
        reads: std::sync::Mutex<Vec<(u64, u64)>>,
        fail: std::sync::atomic::AtomicBool,
    }

    impl RecordingReader {
        fn new(data: Vec<u8>) -> Self {
            Self {
                data,
                reads: std::sync::Mutex::new(Vec::new()),
                fail: std::sync::atomic::AtomicBool::new(false),
            }
        }

        fn clear_reads(&self) {
            self.reads.lock().unwrap().clear();
        }

        fn recorded(&self) -> Vec<(u64, u64)> {
            self.reads.lock().unwrap().clone()
        }

        fn set_fail(&self, value: bool) {
            self.fail.store(value, std::sync::atomic::Ordering::SeqCst);
        }
    }

    impl RangeReader for RecordingReader {
        fn read(&self, offset: u64, len: u64, dest: &mut [u8]) -> i32 {
            self.reads.lock().unwrap().push((offset, len));
            if self.fail.load(std::sync::atomic::Ordering::SeqCst) {
                return -1;
            }
            let start = offset as usize;
            let end = match start.checked_add(len as usize) {
                Some(value) => value,
                None => return -1,
            };
            if end > self.data.len() || dest.len() < len as usize {
                return -1;
            }
            dest[..len as usize].copy_from_slice(&self.data[start..end]);
            0
        }
    }

    #[test]
    fn eic_compute_reads_are_within_plan() {
        let bytes = windowed_ion_bytes();
        let file_len = bytes.len() as u64;
        let reader = std::sync::Arc::new(RecordingReader::new(bytes));
        let read_source = reader.clone();
        let mut ion = ionic::ion::IonReader::open_source(
            Arc::new(ionic::ion::CallbackSource::new(move |range| {
                read_range(&*read_source, range)
            })) as Arc<dyn ionic::ion::ReadBytes>,
            ReadOptions::default(),
        )
        .expect("open_remote failed");

        let plan = plan_eic_ranges(&mut ion, 500.0, 0.0, 10.0, 20.0, 0.005)
            .expect("plan_eic_ranges failed");
        assert!(!plan.is_empty());

        let planned_bytes: u64 = plan.iter().map(|range| range.length).sum();
        assert!(planned_bytes < file_len);

        reader.clear_reads();

        let eic = calculate_eic_dispatcher(
            &mut EicReader::Ion(&mut ion),
            500.0,
            FromTo {
                from: 0.0,
                to: 10.0,
            },
            EicOptions {
                ppm_tolerance: 20.0,
                mz_tolerance: 0.005,
                ..Default::default()
            },
        )
        .expect("calculate_eic failed");
        assert!(!eic.y.is_empty());

        let compute_reads = reader.recorded();
        assert!(
            !compute_reads.is_empty(),
            "compute recorded no reads — prefetch check would be vacuous"
        );

        for (offset, len) in &compute_reads {
            let contained = plan
                .iter()
                .any(|r| r.offset <= *offset && offset.saturating_add(*len) <= r.offset + r.length);
            assert!(
                contained,
                "compute read (offset={offset}, len={len}) is not contained in any planned range"
            );
        }
    }

    #[test]
    fn open_remote_fails_cleanly_when_range_read_fails() {
        let reader = std::sync::Arc::new(RecordingReader::new(windowed_ion_bytes()));
        reader.set_fail(true);
        let read_source = reader.clone();
        let result = ionic::ion::IonReader::open_source(
            Arc::new(ionic::ion::CallbackSource::new(move |range| {
                read_range(&*read_source, range)
            })) as Arc<dyn ionic::ion::ReadBytes>,
            ReadOptions::default(),
        );
        assert!(result.is_err());
    }

    #[test]
    fn eic_compute_errors_cleanly_when_data_read_fails() {
        let reader = std::sync::Arc::new(RecordingReader::new(windowed_ion_bytes()));
        let read_source = reader.clone();
        let mut ion = ionic::ion::IonReader::open_source(
            Arc::new(ionic::ion::CallbackSource::new(move |range| {
                read_range(&*read_source, range)
            })) as Arc<dyn ionic::ion::ReadBytes>,
            ReadOptions::default(),
        )
        .expect("open_remote failed");

        plan_eic_ranges(&mut ion, 500.0, 0.0, 10.0, 20.0, 0.005).expect("plan_eic_ranges failed");

        reader.set_fail(true);

        let result = calculate_eic_dispatcher(
            &mut EicReader::Ion(&mut ion),
            500.0,
            FromTo {
                from: 0.0,
                to: 10.0,
            },
            EicOptions {
                ppm_tolerance: 20.0,
                mz_tolerance: 0.005,
                ..Default::default()
            },
        );
        assert!(result.is_err());
    }

    struct PlannedOnly {
        data: Vec<u8>,
        allowed: Vec<(u64, u64)>,
    }

    impl RangeReader for PlannedOnly {
        fn read(&self, offset: u64, len: u64, dest: &mut [u8]) -> i32 {
            let inside = self
                .allowed
                .iter()
                .any(|(start, size)| offset >= *start && offset + len <= start + size);
            if !inside {
                return -1;
            }
            let start = offset as usize;
            dest[..len as usize].copy_from_slice(&self.data[start..start + len as usize]);
            0
        }
    }

    fn open_serving_only(allowed: Vec<(u64, u64)>, data: Vec<u8>) -> bool {
        let reader = PlannedOnly { data, allowed };
        ionic::ion::IonReader::open_source(
            Arc::new(ionic::ion::CallbackSource::new(move |range| {
                read_range(&reader, range)
            })) as Arc<dyn ionic::ion::ReadBytes>,
            ReadOptions::default(),
        )
        .is_ok()
    }

    fn planned_ranges_for(bytes: &[u8]) -> Vec<(u64, u64)> {
        let out = new_buf_out();
        let code = unsafe { plan_open(bytes.as_ptr(), 1024, out) };
        assert_eq!(code, OK);
        let packed = unsafe { slice::from_raw_parts((*out).ptr, (*out).len).to_vec() };
        drop_buf_out(out);
        decode_byte_ranges(&packed)
    }

    #[test]
    fn planned_ranges_plus_the_trailer_are_enough_to_open() {
        let bytes = windowed_ion_bytes();
        let mut planned = planned_ranges_for(&bytes);
        planned.push((bytes.len() as u64 - 8, 8));
        assert!(
            open_serving_only(planned, bytes),
            "serving the planned ranges plus the trailer still did not open, so something else is unplanned"
        );
    }

    fn spectrum_in_seconds(index: usize, rt: f64) -> Spectrum {
        let mut spectrum = spectrum(index, rt);
        spectrum.scan_list = Some(ScanList {
            count: Some(1),
            scans: vec![Scan {
                cv_params: vec![cv_param_with_unit(
                    "MS:1000016",
                    &rt.to_string(),
                    "UO:0000010",
                )],
                ..Default::default()
            }],
            ..Default::default()
        });
        spectrum
    }

    fn ion_in_seconds() -> Vec<u8> {
        let mzml = MzML {
            run: Run {
                id: "test".to_string(),
                spectrum_list: Some(SpectrumList {
                    count: Some(3),
                    spectra: vec![
                        spectrum_in_seconds(0, 60.0),
                        spectrum_in_seconds(1, 120.0),
                        spectrum_in_seconds(2, 180.0),
                    ],
                    ..Default::default()
                }),
                ..Default::default()
            },
            ..Default::default()
        };
        let options = WriteOptions {
            compression_level: 0,
            force_f32: false,
            block_size: ION_BLOCK_SIZE_BYTES,
            parallel: true,
            section_storage: SectionStorage::Memory,
            mz_window: DEFAULT_MZ_WINDOW,
        };
        let mut bytes = Vec::new();
        write_mzml_to_ion(&mzml, options, &mut bytes).expect("ion encode should succeed");
        bytes
    }

    fn reader_for(bytes: Vec<u8>) -> IonReader {
        let bytes_arc = Arc::from(bytes.into_boxed_slice());
        IonReader::open_source(
            Arc::new(ionic::ion::BytesSource::new(bytes_arc)) as Arc<dyn ionic::ion::ReadBytes>,
            ReadOptions::default(),
        )
        .expect("open_source failed")
    }

    #[test]
    fn plan_window_ranges_is_single_range_on_new_files() {
        let mut ion = reader_for(windowed_ion_bytes());
        let plan = plan_window_ranges(&mut ion, 0.0, 10.0, 499.0, 501.0)
            .expect("plan_window_ranges failed");
        assert_eq!(
            plan.len(),
            1,
            "window-major plan should be one run: {plan:?}"
        );
    }

    #[test]
    fn plan_window_ranges_converts_rt_units() {
        let mut ion = reader_for(ion_in_seconds());
        let from_minutes =
            plan_window_ranges(&mut ion, 1.0, 3.0, 499.0, 501.0).expect("minute plan failed");
        let from_seconds = ion
            .eic_byte_ranges(
                ionic::ion::Range {
                    from: 499.0,
                    to: 501.0,
                },
                Some(ionic::ion::Range {
                    from: 60.0,
                    to: 180.0,
                }),
                0,
            )
            .expect("second plan failed");
        assert!(!from_minutes.is_empty());
        assert_eq!(
            from_minutes
                .iter()
                .map(|r| (r.offset, r.length))
                .collect::<Vec<_>>(),
            from_seconds
                .iter()
                .map(|r| (r.offset, r.length))
                .collect::<Vec<_>>()
        );
    }

    fn covers(planned: &[(u64, u64)], range: &ionic::ion::ByteRange) -> bool {
        planned.iter().any(|(offset, length)| {
            *offset <= range.offset && offset + length >= range.offset + range.length
        })
    }

    #[test]
    fn plan_open_coalesces() {
        let small = windowed_ion_bytes();
        let planned = planned_ranges_for(&small);
        assert_eq!(
            planned,
            vec![(0, small.len() as u64)],
            "a file under the open gap collapses to one range"
        );
        for range in open_ranges(&small[..1024]).expect("open_ranges") {
            assert!(covers(&planned, &range), "{range:?} is not covered");
        }

        let large = std::fs::read("tests/fixtures/test.ion").expect("fixture");
        let planned = planned_ranges_for(&large);
        assert_eq!(planned.len(), 2, "expected header and tail: {planned:?}");
        assert_eq!(planned[0], (0, 1024));
        assert_eq!(planned[1].0 + planned[1].1, large.len() as u64);
        for range in open_ranges(&large[..1024]).expect("open_ranges") {
            assert!(covers(&planned, &range), "{range:?} is not covered");
        }
    }

    #[test]
    fn open_succeeds_using_only_the_ranges_plan_open_asked_for() {
        let bytes = windowed_ion_bytes();
        let planned = planned_ranges_for(&bytes);
        assert!(
            open_serving_only(planned, bytes),
            "plan_open did not list every range the open needs"
        );
    }

    fn scan_queries() -> Vec<(&'static str, ScanQuery)> {
        vec![
            (
                "rt_range",
                ScanQuery::RtRange(FromTo {
                    from: 0.0,
                    to: 10.0,
                }),
            ),
            ("closest_rt", ScanQuery::ClosestRt(5.0)),
            (
                "mz_range",
                ScanQuery::MzRange(FromTo {
                    from: 0.0,
                    to: 1000.0,
                }),
            ),
            ("closest_mz", ScanQuery::ClosestMz(500.0)),
        ]
    }

    fn open_recording_ion(reader: &Arc<RecordingReader>) -> ionic::ion::IonReader {
        let read_source = reader.clone();
        ionic::ion::IonReader::open_source(
            Arc::new(ionic::ion::CallbackSource::new(move |range| {
                read_range(&*read_source, range)
            })) as Arc<dyn ionic::ion::ReadBytes>,
            ReadOptions::default(),
        )
        .expect("open_source failed")
    }

    #[test]
    fn scan_reads_are_within_plan() {
        for (name, query) in scan_queries() {
            let reader = Arc::new(RecordingReader::new(windowed_ion_bytes()));
            let mut ion = open_recording_ion(&reader);

            let plan = plan_scan_ranges(&mut ion, query, TimeUnit::Minutes, 0)
                .expect("plan_scan_ranges failed");

            reader.clear_reads();
            let (_, scans) = get_scans_rs(&mut ion, query, TimeUnit::Minutes, 0);
            if scans.is_empty() {
                assert!(
                    plan.is_empty(),
                    "{name}: nothing was loaded but the plan asked for bytes"
                );
                continue;
            }
            assert!(
                !plan.is_empty(),
                "{name}: scans loaded but the plan is empty"
            );

            let compute_reads = reader.recorded();
            assert!(
                !compute_reads.is_empty(),
                "{name}: compute recorded no reads, the check would be vacuous"
            );

            for (offset, len) in &compute_reads {
                let contained = plan.iter().any(|range| {
                    range.offset <= *offset
                        && offset.saturating_add(*len) <= range.offset + range.length
                });
                assert!(
                    contained,
                    "{name}: read (offset={offset}, len={len}) is outside every planned range"
                );
            }
        }
    }

    struct PlanGate {
        data: Vec<u8>,
        allowed: std::sync::Mutex<Vec<(u64, u64)>>,
        gated: std::sync::atomic::AtomicBool,
    }

    impl PlanGate {
        fn new(data: Vec<u8>) -> Self {
            Self {
                data,
                allowed: std::sync::Mutex::new(Vec::new()),
                gated: std::sync::atomic::AtomicBool::new(false),
            }
        }

        fn allow_only(&self, ranges: Vec<(u64, u64)>) {
            *self.allowed.lock().unwrap() = ranges;
            self.gated.store(true, std::sync::atomic::Ordering::SeqCst);
        }
    }

    impl RangeReader for PlanGate {
        fn read(&self, offset: u64, len: u64, dest: &mut [u8]) -> i32 {
            if self.gated.load(std::sync::atomic::Ordering::SeqCst) {
                let allowed = self.allowed.lock().unwrap();
                let inside = allowed
                    .iter()
                    .any(|(start, size)| offset >= *start && offset + len <= start + size);
                if !inside {
                    return -1;
                }
            }
            let start = offset as usize;
            dest[..len as usize].copy_from_slice(&self.data[start..start + len as usize]);
            0
        }
    }

    #[test]
    fn scans_load_using_only_the_planned_ranges() {
        for (name, query) in scan_queries() {
            let bytes = windowed_ion_bytes();
            let reader = Arc::new(RecordingReader::new(bytes.clone()));
            let mut ion = open_recording_ion(&reader);
            let (_, expected) = get_scans_rs(&mut ion, query, TimeUnit::Minutes, 0);
            if expected.is_empty() {
                continue;
            }

            let gate = Arc::new(PlanGate::new(bytes));
            let read_source = gate.clone();
            let mut served = ionic::ion::IonReader::open_source(
                Arc::new(ionic::ion::CallbackSource::new(move |range| {
                    read_range(&*read_source, range)
                })) as Arc<dyn ionic::ion::ReadBytes>,
                ReadOptions::default(),
            )
            .expect("open_source failed");

            let plan = plan_scan_ranges(&mut served, query, TimeUnit::Minutes, 0)
                .expect("plan_scan_ranges failed");
            gate.allow_only(plan.iter().map(|r| (r.offset, r.length)).collect());

            let (_, got) = get_scans_rs(&mut served, query, TimeUnit::Minutes, 0);
            assert_eq!(
                got.len(),
                expected.len(),
                "{name}: serving only the planned ranges truncated the result"
            );
        }
    }

    #[test]
    fn plan_scans_rejects_a_bad_query() {
        let out = new_buf_out();
        let handle = open_windowed_ion();
        assert_eq!(
            unsafe { plan_scans(handle, 9, 0.0, 10.0, 0, out) },
            ERR_INVALID_ARGS
        );
        assert_eq!(
            unsafe { plan_scans(std::ptr::null_mut(), 0, 0.0, 10.0, 0, out) },
            ERR_INVALID_ARGS
        );
        assert_eq!(
            unsafe { plan_scans(handle, 0, 0.0, 10.0, 0, std::ptr::null_mut()) },
            ERR_INVALID_ARGS
        );
        unsafe { free_mzml(handle) };
        drop_buf_out(out);
    }

    fn noisy_signal() -> Vec<f32> {
        let mut signal = vec![0.0f32; 200];
        for (index, value) in signal.iter_mut().enumerate() {
            *value = 10.0 + ((index % 7) as f32);
        }
        for index in 90..100 {
            signal[index] = 900.0;
        }
        signal
    }

    #[test]
    fn find_noise_level_reports_width_and_intensity() {
        let signal = noisy_signal();
        let expected = find_noise_level_rs(signal.as_slice());
        let mut width = 0usize;
        let mut intensity = 0.0f64;
        let code =
            unsafe { find_noise_level(signal.as_ptr(), signal.len(), &mut width, &mut intensity) };
        assert_eq!(code, OK);
        assert_eq!(width, expected.width);
        assert_eq!(intensity, expected.intensity);
    }

    #[test]
    fn find_noise_level_rejects_missing_outputs() {
        let signal = noisy_signal();
        let mut width = 0usize;
        let mut intensity = 0.0f64;
        assert_eq!(
            unsafe {
                find_noise_level(
                    signal.as_ptr(),
                    signal.len(),
                    std::ptr::null_mut(),
                    &mut intensity,
                )
            },
            ERR_INVALID_ARGS
        );
        assert_eq!(
            unsafe {
                find_noise_level(
                    signal.as_ptr(),
                    signal.len(),
                    &mut width,
                    std::ptr::null_mut(),
                )
            },
            ERR_INVALID_ARGS
        );
        assert_eq!(
            unsafe { find_noise_level(std::ptr::null(), signal.len(), &mut width, &mut intensity) },
            ERR_INVALID_ARGS
        );
        assert_eq!(
            unsafe { find_noise_level(signal.as_ptr(), 0, &mut width, &mut intensity) },
            ERR_INVALID_ARGS
        );
    }
}
