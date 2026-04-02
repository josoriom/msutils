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
from msutils.file import MzMlFile

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
    intensity_threshold: float = 150.0
    width_threshold: int       = 5
    auto_noise: bool           = True
    auto_baseline: bool        = True
    sn_ratio: float            = 1.0


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


def _bin_to_json_raw(file: MzMlFile) -> Any:
    abi = file._abi
    buf = _Buf()
    _check("bin_to_json", abi.bin_to_json(file.handle, ctypes.byref(buf)))
    return buf_to_json(abi, buf)


def _bin_to_mzml_raw(file: MzMlFile) -> str:
    abi = file._abi
    buf = _Buf()
    _check("bin_to_mzml", abi.bin_to_mzml(file.handle, ctypes.byref(buf)))
    return buf_to_str(abi, buf)


def _mzml_to_bin_raw(file: MzMlFile, level: int, f32_compress: bool) -> bytes:
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


def parse_mzml(data: bytes) -> MzMlFile:
    abi = _get_abi()
    arr = (c_uint8 * len(data)).from_buffer_copy(data)
    ptr = ctypes.c_void_p()
    _check("parse_mzml", abi.parse_mzml(
        ctypes.cast(arr, POINTER(c_uint8)), len(data),
        ctypes.byref(ptr),
    ))
    return MzMlFile(ptr, abi)


def parse_bin(data: bytes, max_cache_size: int = 0) -> MzMlFile:
    abi = _get_abi()
    if not isinstance(max_cache_size, int) or max_cache_size < 0:
        raise ValueError("parse_bin: max_cache_size must be a non-negative integer")
    arr = (c_uint8 * len(data)).from_buffer_copy(data)
    ptr = ctypes.c_void_p()
    _check("parse_bin", abi.parse_bin(
        ctypes.cast(arr, POINTER(c_uint8)),
        c_size_t(len(data)),
        c_size_t(max_cache_size),
        ctypes.byref(ptr),
    ))
    return MzMlFile(ptr, abi)


def bin_to_json(file: MzMlFile) -> Any:
    return _bin_to_json_raw(file)


def convert_bin_to_mzml(file: MzMlFile) -> str:
    return _bin_to_mzml_raw(file)


def mzml_to_bin(
    file: MzMlFile,
    level: int = 12,
    f32_compress: bool = False,
) -> bytes:
    if not isinstance(level, int) or not (0 <= level <= 22):
        raise ValueError("mzml_to_bin: level must be an integer in [0, 22]")
    if not isinstance(f32_compress, bool):
        raise TypeError("mzml_to_bin: f32_compress must be a bool")
    return _mzml_to_bin_raw(file, level, f32_compress)


def calculate_eic(
    file: MzMlFile,
    target_mz: float,
    from_rt: float,
    to_rt: float,
    ppm_tol: float = 20.0,
    mz_tol: float  = 0.005,
) -> Dict[str, np.ndarray]:
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
    abi = _get_abi()
    arr = np.asarray(y, dtype=np.float32)
    ptr = arr.ctypes.data_as(POINTER(c_float))
    return float(abi.find_noise_level(ptr, len(arr)))


def calculate_baseline(
    y: Sequence,
    baseline_window: int = 0,
    baseline_window_factor: int = 0,
) -> np.ndarray:
    abi = _get_abi()
    y_array = _to_f64_array(y)
    buf = _Buf()
    _check("calculate_baseline", abi.calculate_baseline(
        _as_f64_ptr(y_array), c_size_t(len(y_array)),
        c_int32(int(baseline_window)), c_int32(int(baseline_window_factor)),
        ctypes.byref(buf),
    ))
    return buf_to_f64(abi, buf)


def collect_scans(
    file: MzMlFile,
    from_rt: float,
    to_rt: float,
    level: int = 1,
) -> List[Dict]:
    abi = _get_abi()
    if not (0 <= level <= 255):
        raise ValueError("collect_scans: level must be in [0, 255]")
    buf = _Buf()
    _check("collect_scans", abi.collect_scans(
        file.handle,
        c_double(from_rt), c_double(to_rt), c_uint8(level),
        ctypes.byref(buf),
    ))
    return buf_to_json(abi, buf)


def get_peaks_from_eic(
    file: MzMlFile,
    targets: Sequence[Dict[str, Any]],
    from_rt: float = 0.5,
    to_rt: float   = 5.0,
    options: Optional[Dict[str, Any]] = None,
    cores: int = 1,
) -> List[Dict]:
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
    file: MzMlFile,
    items: Sequence[Dict[str, Any]],
    options: Optional[Dict[str, Any]] = None,
    cores: int = 1,
) -> List[Dict]:
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
    file: MzMlFile,
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
    abi = _get_abi()

    eic_ppm  = _get_float(eic,  "ppm_tolerance", _EicDefaults.ppm_tolerance)
    eic_mz   = _get_float(eic,  "mz_tolerance",  _EicDefaults.mz_tolerance)
    g_start  = _get_float(grid, "start",          _GridDefaults.start)
    g_end    = _get_float(grid, "end",            _GridDefaults.end)
    g_step   = _get_float(grid, "step_size",      _GridDefaults.step_size)
    if g_step <= 0:
        g_step = math.nan

    peak_defaults = {
        "intensity_threshold": _FindPeakDefaults.intensity_threshold,
        "width_threshold":     _FindPeakDefaults.width_threshold,
        "auto_noise":          _FindPeakDefaults.auto_noise,
        "auto_baseline":       _FindPeakDefaults.auto_baseline,
        "sn_ratio":            _FindPeakDefaults.sn_ratio,
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
    file: MzMlFile,
    targets: Sequence[Dict[str, Any]],
    *,
    scan_eic: Optional[Dict[str, float]] = None,
    eic: Optional[Dict[str, float]]      = None,
    options: Optional[Dict[str, Any]]    = None,
    cores: int = 1,
) -> List[Dict]:
    abi = _get_abi()

    scan_ppm = _get_tolerance(scan_eic, "ppm_tolerance", 10.0)
    scan_mz  = _get_tolerance(scan_eic, "mz_tolerance",  0.003)
    eic_ppm  = _get_tolerance(eic,      "ppm_tolerance", 20.0)
    eic_mz   = _get_tolerance(eic,      "mz_tolerance",  0.005)

    peak_defaults = {
        "intensity_threshold": _FindPeakDefaults.intensity_threshold,
        "width_threshold":     _FindPeakDefaults.width_threshold,
        "auto_noise":          _FindPeakDefaults.auto_noise,
        "auto_baseline":       _FindPeakDefaults.auto_baseline,
        "sn_ratio":            _FindPeakDefaults.sn_ratio,
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
        "intensity_threshold": _FindPeakDefaults.intensity_threshold,
        "width_threshold":     _FindPeakDefaults.width_threshold,
        "auto_noise":          _FindPeakDefaults.auto_noise,
        "auto_baseline":       _FindPeakDefaults.auto_baseline,
        "sn_ratio":            _FindPeakDefaults.sn_ratio,
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