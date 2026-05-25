from __future__ import annotations

import ctypes
import math
from ctypes import POINTER, c_double, c_float, c_int32, c_size_t, c_uint8, c_uint32
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

from msutils._bridge import (
    PeakOptions,
    _ABI,
    _Buf,
    _check,
    buf_to_bytes,
    buf_to_f64,
    buf_to_json,
    buf_to_str,
)
from msutils._pack import (
    encode_target_ids,
    pack_peak_options,
    unpack_targets,
)
from msutils._shared import to_cores
from msutils.file import SampleFile

_abi: Optional[_ABI] = None


def _get_abi() -> _ABI:
    if _abi is None:
        raise RuntimeError("msutils backend not loaded.")
    return _abi


def _set_abi(abi: _ABI) -> None:
    global _abi
    _abi = abi


class _EicDefaults:
    mz_tolerance: float  = 0.005
    ppm_tolerance: float = 10.0

class _GridDefaults:
    start: float     = 40.0
    end: float       = 1000.0
    step_size: float = 0.005

class _GroupingDefaults:
    ppm_tolerance: float = 5.0
    mz_tolerance: float  = 0.0025
    rt_tolerance: float  = 0.05
    frequency: int       = 1

class _FindPeakDefaults:
    min_intensity: float        = 150.0
    min_peak_width_points: int  = 5
    auto_noise: bool            = True
    auto_baseline: bool         = True
    min_snr: float              = 1.0


def _as_f64_ptr(array: np.ndarray) -> POINTER(c_double):
    return array.ctypes.data_as(POINTER(c_double))

def _as_u32_ptr(array: np.ndarray) -> POINTER(c_uint32):
    return array.ctypes.data_as(POINTER(c_uint32))

def _as_u8_ptr(data: bytes):
    arr = (c_uint8 * len(data)).from_buffer_copy(data)
    return ctypes.cast(arr, POINTER(c_uint8)), arr

def _as_opts_ptr(peak_options: Optional[PeakOptions]):
    if peak_options is None:
        return ctypes.cast(None, POINTER(PeakOptions))
    return ctypes.byref(peak_options)

def _to_f64_array(sequence: Sequence) -> np.ndarray:
    return np.asarray(sequence, dtype=np.float64)


def _get_float(mapping: Optional[Dict], key: str, default: float) -> float:
    if mapping is None:
        return default
    value = mapping.get(key, default)
    parsed = float(value)
    return parsed if math.isfinite(parsed) and parsed >= 0 else math.nan


def _get_tolerance(mapping: Optional[Dict], key: str, default: float) -> float:
    if mapping is None:
        return default
    value = mapping.get(key, default)
    parsed = float(value)
    if not math.isfinite(parsed) or parsed < 0:
        raise TypeError(f"{key} must be a finite non-negative number")
    return parsed


def _get_numeric(mapping: Optional[Dict], key: str, default: float) -> float:
    if mapping is None:
        return default
    value = mapping.get(key, default)
    parsed = float(value)
    return parsed if math.isfinite(parsed) else default


def _ion_to_json_raw(file: SampleFile) -> Any:
    abi = file._abi
    buf = _Buf()
    _check("bin_to_json", abi.bin_to_json(file.handle, ctypes.byref(buf)))
    return buf_to_json(abi, buf)


def _ion_to_mzml_raw(file: SampleFile) -> str:
    abi = file._abi
    buf = _Buf()
    _check("bin_to_mzml", abi.bin_to_mzml(file.handle, ctypes.byref(buf)))
    return buf_to_str(abi, buf)


def _mzml_to_ion_raw(file: SampleFile, level: int, f32_compress: bool) -> bytes:
    abi = file._abi
    if not (0 <= level <= 22):
        raise ValueError("level must be in [0, 22]")
    buf = _Buf()
    _check(
        "mzml_to_bin",
        abi.mzml_to_bin(
            file.handle,
            ctypes.byref(buf),
            c_uint8(level),
            c_uint8(1 if f32_compress else 0),
        ),
    )
    return buf_to_bytes(abi, buf)


def parse_mzml(data: bytes) -> SampleFile:
    """Load an mzML file.

    Args:
        data: Raw bytes from one mzML file.

    Returns: A SampleFile for use in other msutils functions.
    """
    abi = _get_abi()
    arr = (c_uint8 * len(data)).from_buffer_copy(data)
    ptr = ctypes.c_void_p()
    _check("parse_mzml", abi.parse_mzml(
        ctypes.cast(arr, POINTER(c_uint8)), len(data),
        ctypes.byref(ptr),
    ))
    return SampleFile(ptr, abi)


def parse_ion(data: bytes, max_cache_size: int = 0) -> SampleFile:
    """Load an ion file.

    Args:
        data: Raw bytes from one ion file.
        max_cache_size: Cache size in bytes. 0 sets no limit. Larger valuesspeed up repeated reads at the cost of RAM (default 0).

    Returns: A SampleFile for use in other msutils functions.
    """
    abi = _get_abi()
    if not isinstance(max_cache_size, int) or max_cache_size < 0:
        raise ValueError("parse_ion: max_cache_size must be a non-negative integer")
    arr = (c_uint8 * len(data)).from_buffer_copy(data)
    ptr = ctypes.c_void_p()
    _check("parse_bin", abi.parse_bin(
        ctypes.cast(arr, POINTER(c_uint8)),
        c_size_t(len(data)),
        c_size_t(max_cache_size),
        ctypes.byref(ptr),
    ))
    return SampleFile(ptr, abi)


def ion_to_json(file: SampleFile) -> Any:
    """Get JSON from a sample.

    Args:
        file: SampleFile from parse_mzml() or parse_ion().

    Returns: Parsed JSON object with all run metadata, spectra, and chromatograms.
    """
    return _ion_to_json_raw(file)


def ion_to_mzml(file: SampleFile) -> str:
    """Get mzML file from a sample.

    Args:
        file: SampleFile from parse_mzml() or parse_ion().

    Returns: A string with valid mzML content.
    """
    return _ion_to_mzml_raw(file)


def mzml_to_ion(
    file: SampleFile,
    level: int = 12,
    f32_compress: bool = False,
) -> bytes:
    """Convert a sample to ion binary bytes.

    Args:
        file: SampleFile from parse_mzml() or parse_ion().
        level: Compression level, integer 0 (none) to 22 (max). Default 12.
        f32_compress: If True, compress intensity values to 32-bit float.

    Returns: Raw ion binary bytes.
    """
    if not isinstance(level, int) or not (0 <= level <= 22):
        raise ValueError("mzml_to_ion: level must be an integer in [0, 22]")
    if not isinstance(f32_compress, bool):
        raise TypeError("mzml_to_ion: f32_compress must be a bool")
    return _mzml_to_ion_raw(file, level, f32_compress)


def calculate_eic(
    file: SampleFile,
    target_mz: float,
    from_rt: float,
    to_rt: float,
    ppm_tol: float = 20.0,
    mz_tol: float  = 0.005,
) -> Dict[str, np.ndarray]:
    """Get an extracted ion chromatogram (EIC).

    Extracts an EIC for one target m/z over a given retention-time range.

    Args:
        file: SampleFile from parse_mzml() or parse_ion().
        target_mz: Target m/z value (Da).
        from_rt: Start retention time (minutes).
        to_rt: End retention time (minutes).
        ppm_tol: Mass tolerance in ppm. Default 20.
        mz_tol: Absolute mass tolerance in Da. Default 0.005.

    Returns: Dict with keys 'x' (retention time array) and 'y' (intensity array).
    """
    abi = _get_abi()
    buf_x, buf_y = _Buf(), _Buf()
    _check("calculate_eic", abi.calculate_eic(
        file.handle,
        c_double(target_mz), c_double(from_rt), c_double(to_rt),
        c_double(ppm_tol), c_double(mz_tol),
        ctypes.byref(buf_x), ctypes.byref(buf_y),
    ))
    return {"x": buf_to_f64(abi, buf_x), "y": buf_to_f64(abi, buf_y)}


def find_peaks(
    x: Sequence,
    y: Sequence,
    options: Optional[Dict[str, Any]] = None,
) -> List[Dict]:
    """Find all peaks in a trace.

    Detects peaks in one numeric trace (x = time, y = intensity).
    Returns all peaks that pass the given filters.

    Args:
        x: Retention-time sequence. Minimum 3 points.
        y: Intensity sequence, same length as x.
        options: Optional dict of peak filter settings. Keys:
            min_intensity (float),
            min_peak_width_points (int),
            noise (float),
            auto_noise (bool, default True),
            auto_baseline (bool, default True),
            lambda_ (int),
            max_iterations (int),
            allow_overlap (bool),
            min_snr (float).

    Returns:
        List of dicts, one per peak. Each dict has keys: rt, from, to, intensity and integral.
    """
    abi = _get_abi()
    x_array, y_array = _to_f64_array(x), _to_f64_array(y)
    if len(x_array) != len(y_array) or len(x_array) < 3:
        raise ValueError("x and y must have equal length >= 3")
    peak_options = pack_peak_options(options)
    buf = _Buf()
    _check("find_peaks", abi.find_peaks(
        _as_f64_ptr(x_array), _as_f64_ptr(y_array), c_size_t(len(x_array)),
        _as_opts_ptr(peak_options), ctypes.byref(buf),
    ))
    return buf_to_json(abi, buf)


def get_peak(
    x: Sequence,
    y: Sequence,
    rt: float,
    range_: float,
    options: Optional[Dict[str, Any]] = None,
) -> Dict:
    """Find one peak near a target RT.

    Finds the best peak within a given retention-time window around one
    target RT. Use this when you already know approximately where to look.

    Args:
        x: Retention-time sequence. Minimum 3 points.
        y: Intensity sequence, same length as x.
        rt: Target retention time.
        range_: Half-width of the search window around rt.
        options: Optional dict of peak filter settings (same keys as find_peaks).

    Returns:
        Dict with peak fields: rt, from, to, intensity, integral.
        Empty dict if no peak was found.
    """
    abi = _get_abi()
    x_array, y_array = _to_f64_array(x), _to_f64_array(y)
    if len(x_array) != len(y_array) or len(x_array) < 3:
        raise ValueError("x and y must have equal length >= 3")
    peak_options = pack_peak_options(options)
    buf = _Buf()
    _check("get_peak", abi.get_peak(
        _as_f64_ptr(x_array),
        _as_f64_ptr(y_array),
        c_size_t(len(x_array)),
        c_double(rt),
        c_double(range_),
        _as_opts_ptr(peak_options),
        ctypes.byref(buf),
    ))
    return buf_to_json(abi, buf)


def find_noise_level(y: Sequence) -> float:
    """Estimate the noise level of a signal.

    Args:
        y: Intensity sequence.

    Returns:
        Estimated noise level as a float.
    """
    abi = _get_abi()
    arr = np.asarray(y, dtype=np.float32)
    ptr = arr.ctypes.data_as(POINTER(c_float))
    return float(abi.find_noise_level(ptr, len(arr)))


def calculate_baseline(
    y: Sequence,
    lambda_: int = 0,
    max_iterations: int = 0,
) -> np.ndarray:
    """Estimate a signal baseline.

    Args:
        y: Numeric signal sequence.
        lambda_: Smoothness parameter for the asymmetric least-squares baseline. 0 uses a default.
        max_iterations: Maximum ALS iterations. 0 uses a default.

    Returns:
        Numpy array, same length as y, with the estimated baseline values.
    """
    abi = _get_abi()
    y_array = _to_f64_array(y)
    buf = _Buf()
    _check("calculate_baseline", abi.calculate_baseline(
        _as_f64_ptr(y_array), c_size_t(len(y_array)),
        c_int32(int(lambda_)), c_int32(int(max_iterations)),
        ctypes.byref(buf),
    ))
    return buf_to_f64(abi, buf)


_QUERY_RT_RANGE   = 0
_QUERY_CLOSEST_RT = 1
_QUERY_MZ_RANGE   = 2
_QUERY_CLOSEST_MZ = 3


def get_scans(
    file: SampleFile,
    *,
    rt_range: Optional[tuple] = None,
    rt: Optional[float] = None,
    mz_range: Optional[tuple] = None,
    mz: Optional[float] = None,
    level: int = 1,
) -> "List[Dict] | Optional[Dict]":
    """Get scans by RT or m/z query.

    Fetches scans from a loaded file. Exactly one query mode must be provided.

    Args:
        file: SampleFile from parse_mzml() or parse_ion().
        rt_range: Tuple (from, to) in minutes — all scans in RT window.
        rt: Target RT in minutes — returns the single closest scan.
        mz_range: Tuple (from, to) in Da — all scans with precursor m/z in window.
        mz: Target precursor m/z in Da — returns the single closest scan.
        level: MS level to query (0–255). Default 1.

    Returns:
        Range query: list of scan dicts, each with rt, mz, intensity, metadata.
        Point query: single scan dict, or None if no scan was found.
    """
    abi = _get_abi()
    if not (0 <= level <= 255):
        raise ValueError("get_scans: level must be in [0, 255]")

    provided = sum(x is not None for x in (rt_range, rt, mz_range, mz))
    if provided != 1:
        raise ValueError("get_scans: exactly one of rt_range, rt, mz_range, mz must be provided")

    if rt_range is not None:
        a, b = float(rt_range[0]), float(rt_range[1])
        if not math.isfinite(a) or not math.isfinite(b):
            raise ValueError("get_scans: rt_range values must be finite")
        query_type = _QUERY_RT_RANGE
    elif rt is not None:
        a = float(rt)
        if not math.isfinite(a):
            raise ValueError("get_scans: rt must be finite")
        query_type, b = _QUERY_CLOSEST_RT, math.nan
    elif mz_range is not None:
        a, b = float(mz_range[0]), float(mz_range[1])
        if not math.isfinite(a) or not math.isfinite(b):
            raise ValueError("get_scans: mz_range values must be finite")
        query_type = _QUERY_MZ_RANGE
    else:
        a = float(mz)  # type: ignore[arg-type]
        if not math.isfinite(a) or a <= 0:
            raise ValueError("get_scans: mz must be a positive finite number")
        query_type, b = _QUERY_CLOSEST_MZ, math.nan

    buf = _Buf()
    _check("get_scans", abi.get_scans(
        file.handle,
        c_uint8(query_type), c_double(a), c_double(b), c_uint8(level),
        ctypes.byref(buf),
    ))
    scans = buf_to_json(abi, buf)
    if query_type in (_QUERY_CLOSEST_RT, _QUERY_CLOSEST_MZ):
        return scans[0] if scans else None
    return scans


def get_peaks_from_eic(
    file: SampleFile,
    targets: Sequence[Dict[str, Any]],
    from_rt: float = 0.5,
    to_rt: float   = 5.0,
    options: Optional[Dict[str, Any]] = None,
    cores: int = 1,
) -> List[Dict]:
    """Find peaks for many targets using EICs.

    Builds an EIC for each target and finds the best peak near the given RT.
    Use this for batch targeted peak extraction from one file.

    Args:
        file: SampleFile from parse_mzml() or parse_ion().
        targets: List of dicts, each with keys: id (str), rt (float),
            mz (float), range (float, RT half-window).
        from_rt: EIC extraction start time (minutes). Default 0.5.
        to_rt: EIC extraction end time (minutes). Default 5.0.
        options: Optional dict of peak filter settings (same keys as find_peaks).
        cores: Number of CPU cores to use. Default 1.

    Returns:
        List of dicts, one per target. Each dict has keys: id, mz, rt,
        from, to, intensity, integral.
    """
    abi = _get_abi()
    rts, mzs, ranges, ids = unpack_targets(targets)
    offsets, lengths, id_buf = encode_target_ids(ids)
    peak_options = pack_peak_options(options)
    buf = _Buf()

    if id_buf:
        id_ptr, id_arr = _as_u8_ptr(id_buf)
        rc = abi.get_peaks_from_eic(
            file.handle,
            _as_f64_ptr(rts), _as_f64_ptr(mzs), _as_f64_ptr(ranges),
            _as_u32_ptr(offsets), _as_u32_ptr(lengths),
            id_ptr, c_size_t(len(id_buf)),
            c_size_t(len(targets)), c_double(from_rt), c_double(to_rt),
            _as_opts_ptr(peak_options), c_size_t(to_cores(cores)),
            ctypes.byref(buf),
        )
    else:
        rc = abi.get_peaks_from_eic(
            file.handle,
            _as_f64_ptr(rts), _as_f64_ptr(mzs), _as_f64_ptr(ranges),
            None, None, None, c_size_t(0),
            c_size_t(len(targets)), c_double(from_rt), c_double(to_rt),
            _as_opts_ptr(peak_options), c_size_t(to_cores(cores)),
            ctypes.byref(buf),
        )
    _check("get_peaks_from_eic", rc)
    return buf_to_json(abi, buf)


def get_peaks_from_chrom(
    file: SampleFile,
    items: Sequence[Dict[str, Any]],
    options: Optional[Dict[str, Any]] = None,
    cores: int = 1,
) -> List[Dict]:
    """Find peaks from stored chromatograms.

    Finds peaks from file-stored chromatograms using index, RT, and RT window.
    Use this when you already know which chromatogram index to read.

    Args:
        file: SampleFile from parse_mzml() or parse_ion().
        items: List of dicts with keys: idx (int, chromatogram index),
            rt (float, target retention time), window or range (float, RT half-window).
        options: Optional dict of peak filter settings (same keys as find_peaks).
        cores: Number of CPU cores to use. Default 1.

    Returns:
        List of dicts, one per input item. Each dict has keys: index, rt,
        from, to, intensity, integral.
    """
    abi = _get_abi()
    n = len(items)
    indices = np.empty(n, dtype=np.uint32)
    rts     = np.empty(n, dtype=np.float64)
    windows = np.empty(n, dtype=np.float64)
    for i, item in enumerate(items):
        raw_idx = item.get("idx", item.get("index", -1))
        idx = int(raw_idx) if raw_idx is not None else -1
        indices[i] = idx if idx >= 0 else 0xFFFFFFFF
        rts[i]     = float(item.get("rt", 0.0))
        windows[i] = float(item.get("window", item.get("range", 0.0)))
    peak_options = pack_peak_options(options)
    buf = _Buf()
    _check("get_peaks_from_chrom", abi.get_peaks_from_chrom(
        file.handle,
        _as_u32_ptr(indices), _as_f64_ptr(rts), _as_f64_ptr(windows),
        c_size_t(n), _as_opts_ptr(peak_options), c_size_t(to_cores(cores)),
        ctypes.byref(buf),
    ))
    return buf_to_json(abi, buf)


def find_features(
    file: SampleFile,
    from_rt: float,
    to_rt: float,
    *,
    eic: Optional[Dict[str, float]] = None,
    grid: Optional[Dict[str, float]] = None,
    options: Optional[Dict[str, Any]] = None,
    cores: int = 1,
    use_gpu: bool = False,
    batch_size: int = 0,
) -> List[Dict]:
    """Find all features in one file.

    Scans the full m/z grid and finds all chromatographic features.
    No target list is needed.

    Args:
        file: SampleFile from parse_mzml() or parse_ion().
        from_rt: Start retention time (minutes).
        to_rt: End retention time (minutes).
        eic: Optional dict with EIC extraction settings. Keys:
            ppm_tolerance (float, default 10.0), mz_tolerance (float, default 0.005).
        grid: Optional dict with m/z grid settings. Keys:
            start (float, default 40.0), end (float, default 1000.0),
            step_size (float, default 0.005).
        options: Optional dict of peak filter settings (same keys as find_peaks).
        cores: Number of CPU cores to use. Default 1.
        use_gpu: Use GPU acceleration if available. Default False.
        batch_size: GPU batch size. 0 selects automatically.

    Returns:
        List of dicts, one per detected feature. Each dict has keys: mz, rt,
        from, to, intensity, integral, n_points.
    """
    abi = _get_abi()

    eic_ppm  = _get_float(eic,  "ppm_tolerance", _EicDefaults.ppm_tolerance)
    eic_mz   = _get_float(eic,  "mz_tolerance",  _EicDefaults.mz_tolerance)
    g_start  = _get_float(grid, "start",          _GridDefaults.start)
    g_end    = _get_float(grid, "end",            _GridDefaults.end)
    g_step   = _get_float(grid, "step_size",      _GridDefaults.step_size)
    if g_step <= 0:
        g_step = math.nan

    peak_defaults = {
        "min_intensity":        _FindPeakDefaults.min_intensity,
        "min_peak_width_points": _FindPeakDefaults.min_peak_width_points,
        "auto_noise":           _FindPeakDefaults.auto_noise,
        "auto_baseline":        _FindPeakDefaults.auto_baseline,
        "min_snr":              _FindPeakDefaults.min_snr,
    }
    if options:
        peak_defaults.update(options)
    peak_options = pack_peak_options(peak_defaults)

    buf = _Buf()
    _check("find_features", abi.find_features(
        file.handle,
        c_double(from_rt), c_double(to_rt),
        c_double(eic_ppm), c_double(eic_mz),
        c_double(g_start), c_double(g_end), c_double(g_step),
        _as_opts_ptr(peak_options), c_int32(to_cores(cores)),
        c_int32(1 if use_gpu else 0), c_int32(int(batch_size)),
        ctypes.byref(buf),
    ))
    return buf_to_json(abi, buf)


def find_feature(
    file: SampleFile,
    targets: Sequence[Dict[str, Any]],
    *,
    scan_eic: Optional[Dict[str, float]] = None,
    eic: Optional[Dict[str, float]]      = None,
    options: Optional[Dict[str, Any]]    = None,
    cores: int = 1,
) -> List[Dict]:
    """Find targeted features by m/z and RT.

    Extracts peaks for one or more targets defined by m/z and RT.
    Use this when you already know what you are looking for.

    Args:
        file: SampleFile from parse_mzml() or parse_ion().
        targets: List of dicts with keys: id (str, optional), rt (float),
            mz (float), range (float, RT half-window).
        scan_eic: Optional dict for spectrum scan tolerance. Keys:
            ppm_tolerance (float), mz_tolerance (float).
        eic: Optional dict for EIC extraction tolerance. Keys:
            ppm_tolerance (float), mz_tolerance (float).
        options: Optional dict of peak filter settings (same keys as find_peaks).
        cores: Number of CPU cores to use. Default 1.

    Returns:
        List of dicts, one per target. Each dict has keys: id, mz, rt,
        from, to, intensity, integral.
    """
    abi = _get_abi()

    scan_ppm = _get_tolerance(scan_eic, "ppm_tolerance", 10.0)
    scan_mz  = _get_tolerance(scan_eic, "mz_tolerance",  0.003)
    eic_ppm  = _get_tolerance(eic,      "ppm_tolerance", 20.0)
    eic_mz   = _get_tolerance(eic,      "mz_tolerance",  0.005)

    peak_defaults = {
        "min_intensity":        _FindPeakDefaults.min_intensity,
        "min_peak_width_points": _FindPeakDefaults.min_peak_width_points,
        "auto_noise":           _FindPeakDefaults.auto_noise,
        "auto_baseline":        _FindPeakDefaults.auto_baseline,
        "min_snr":              _FindPeakDefaults.min_snr,
    }
    if options:
        peak_defaults.update(options)
    peak_options = pack_peak_options(peak_defaults)

    rts, mzs, ranges, ids = unpack_targets(targets, default_range=0.5)
    offsets, lengths, id_buf = encode_target_ids(ids)
    buf = _Buf()

    if id_buf:
        id_ptr, id_arr = _as_u8_ptr(id_buf)
        rc = abi.find_feature(
            file.handle,
            _as_f64_ptr(rts), _as_f64_ptr(mzs), _as_f64_ptr(ranges),
            _as_u32_ptr(offsets), _as_u32_ptr(lengths),
            id_ptr, c_size_t(len(id_buf)),
            c_size_t(len(targets)), c_size_t(to_cores(cores)),
            c_double(scan_ppm), c_double(scan_mz),
            c_double(eic_ppm),  c_double(eic_mz),
            _as_opts_ptr(peak_options), ctypes.byref(buf),
        )
    else:
        rc = abi.find_feature(
            file.handle,
            _as_f64_ptr(rts), _as_f64_ptr(mzs), _as_f64_ptr(ranges),
            None, None, None, c_size_t(0),
            c_size_t(len(targets)), c_size_t(to_cores(cores)),
            c_double(scan_ppm), c_double(scan_mz),
            c_double(eic_ppm),  c_double(eic_mz),
            _as_opts_ptr(peak_options), ctypes.byref(buf),
        )
    _check("find_feature", rc)
    return buf_to_json(abi, buf)


def get_features(
    directory_path: str,
    from_rt: float,
    to_rt: float,
    *,
    eic: Optional[Dict[str, float]]      = None,
    grid: Optional[Dict[str, float]]     = None,
    grouping: Optional[Dict[str, Any]]   = None,
    options: Optional[Dict[str, Any]]    = None,
    cores: int = 1,
    use_gpu: bool = False,
    batch_size: int = 0,
) -> List[Dict]:
    """Find and align features across many files.

    Detects features in every file in a folder and groups them by m/z and RT.
    Returns a consensus feature table aligned across all samples.

    Args:
        directory_path: Path to a folder with mzML or ion files.
        from_rt: Start retention time (minutes).
        to_rt: End retention time (minutes).
        eic: Optional dict for EIC extraction settings. Keys:
            ppm_tolerance (float, default 10.0), mz_tolerance (float, default 0.005).
        grid: Optional dict for m/z grid settings. Keys:
            start (float, default 40.0), end (float, default 1000.0),
            step_size (float, default 0.005).
        grouping: Optional dict for cross-sample alignment settings. Keys:
            ppm_tolerance (float), mz_tolerance (float), rt_tolerance (float),
            frequency (int, min samples a feature must appear in).
        options: Optional dict of peak filter settings (same keys as find_peaks).
        cores: Number of CPU cores to use. Default 1.
        use_gpu: Use GPU acceleration if available. Default False.
        batch_size: GPU batch size. 0 selects automatically.

    Returns:
        List of dicts, one per consensus feature. Each dict has the
        following keys:

        - ``mz`` (float): median m/z of the feature (Da).
        - ``rt`` (float): median retention time (minutes).
        - ``from`` (float): start of the retention-time window (minutes).
        - ``to`` (float): end of the retention-time window (minutes).
        - ``intensity`` (float): median apex intensity.
        - ``integral`` (float): median peak area.
        - ``frequency`` (float): fraction of input samples in which this
          feature was detected (``n_samples / total_samples``, in
          ``[0.0, 1.0]``).
    """
    abi = _get_abi()

    if not directory_path:
        raise ValueError("get_features: directory_path is required")

    eic_ppm    = _get_numeric(eic,      "ppm_tolerance", _EicDefaults.ppm_tolerance)
    eic_mz     = _get_numeric(eic,      "mz_tolerance",  _EicDefaults.mz_tolerance)
    grid_start = _get_numeric(grid,     "start",          _GridDefaults.start)
    grid_end   = _get_numeric(grid,     "end",            _GridDefaults.end)
    grid_step  = _get_numeric(grid,     "step_size",      _GridDefaults.step_size)
    group_ppm  = _get_numeric(grouping, "ppm_tolerance",  _GroupingDefaults.ppm_tolerance)
    group_mz   = _get_numeric(grouping, "mz_tolerance",   _GroupingDefaults.mz_tolerance)
    group_rt   = _get_numeric(grouping, "rt_tolerance",   _GroupingDefaults.rt_tolerance)
    prevalence = int((grouping or {}).get(
        "frequency", (grouping or {}).get("prevalence", _GroupingDefaults.frequency)
    ))

    peak_defaults = {
        "min_intensity":        _FindPeakDefaults.min_intensity,
        "min_peak_width_points": _FindPeakDefaults.min_peak_width_points,
        "auto_noise":           _FindPeakDefaults.auto_noise,
        "auto_baseline":        _FindPeakDefaults.auto_baseline,
        "min_snr":              _FindPeakDefaults.min_snr,
    }
    if options:
        peak_defaults.update(options)
    peak_options = pack_peak_options(peak_defaults)

    buf = _Buf()
    _check("get_features", abi.get_features(
        directory_path.encode("utf-8"),
        c_double(from_rt), c_double(to_rt),
        c_double(eic_ppm), c_double(eic_mz),
        c_double(grid_start), c_double(grid_end), c_double(grid_step),
        c_double(group_ppm), c_double(group_mz), c_double(group_rt),
        c_int32(prevalence),
        _as_opts_ptr(peak_options), c_int32(to_cores(cores)),
        c_int32(1 if use_gpu else 0), c_int32(int(batch_size)),
        ctypes.byref(buf),
    ))
    return buf_to_json(abi, buf)
