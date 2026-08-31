from __future__ import annotations

import ctypes
import math
import os
from ctypes import POINTER, c_double, c_float, c_int32, c_size_t, c_uint8, c_uint32
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

from quantion._bridge import (
    PeakOptions,
    READ_RANGE,
    _ABI,
    _Buf,
    _check,
    buf_to_bytes,
    buf_to_f64,
    buf_to_records,
    buf_to_scans,
    buf_to_str,
)
from quantion._pack import (
    encode_target_ids,
    pack_peak_options,
    unpack_targets,
)
from quantion._shared import to_cores
from quantion.file import SampleFile

_abi: Optional[_ABI] = None


def _get_abi() -> _ABI:
    if _abi is None:
        raise RuntimeError("quantion backend not loaded.")
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
    min_intensity: float        = 500.0
    min_peak_width_points: int  = 3
    auto_noise: bool            = True
    auto_baseline: bool         = True
    min_snr: float              = 2.0
    min_r2: float               = 0.0
    shape: str                  = "emg"


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
    return np.ascontiguousarray(sequence, dtype=np.float64)


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

    Returns: A SampleFile for use in other quantion functions.
    """
    abi = _get_abi()
    arr = (c_uint8 * len(data)).from_buffer_copy(data)
    ptr = ctypes.c_void_p()
    _check("parse_mzml", abi.parse_mzml(
        ctypes.cast(arr, POINTER(c_uint8)), len(data),
        ctypes.byref(ptr),
    ))
    return SampleFile(ptr, abi)


def parse_ion(source, max_cache_size: int = 0) -> SampleFile:
    """Load an ion file from raw bytes, a path, or a URL.

    Picks the right one for you: bytes go to parse_ion_raw, a string that looks
    like a URL goes to parse_ion_remote, any other string goes to parse_ion_path.
    """
    if isinstance(source, (str, os.PathLike)):
        text = os.fspath(source)
        if text.startswith("http://") or text.startswith("https://"):
            return parse_ion_remote(text, max_cache_size)
        return parse_ion_path(text, max_cache_size)
    return parse_ion_raw(source, max_cache_size)


def parse_ion_raw(source, max_cache_size: int = 0) -> SampleFile:
    """Load an ion file from raw bytes already in memory."""
    if not isinstance(source, (bytes, bytearray, memoryview)):
        raise TypeError("parse_ion_raw: source must be ion bytes")
    if not isinstance(max_cache_size, int) or max_cache_size < 0:
        raise ValueError("parse_ion_raw: max_cache_size must be a non-negative integer")
    abi = _get_abi()
    arr = (c_uint8 * len(source)).from_buffer_copy(source)
    ptr = ctypes.c_void_p()
    _check("parse_bin", abi.parse_bin(
        ctypes.cast(arr, POINTER(c_uint8)),
        c_size_t(len(source)),
        c_size_t(max_cache_size),
        ctypes.byref(ptr),
    ))
    return SampleFile(ptr, abi)


def parse_ion_path(path, max_cache_size: int = 0) -> SampleFile:
    """Load an ion file from a path, reading only what it needs from disk."""
    if not isinstance(max_cache_size, int) or max_cache_size < 0:
        raise ValueError("parse_ion_path: max_cache_size must be a non-negative integer")
    file_path = os.fspath(path)
    if not file_path:
        raise ValueError("parse_ion_path: path must be a non-empty string")
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"parse_ion_path: no file at {file_path}")
    abi = _get_abi()
    ptr = ctypes.c_void_p()
    _check("parse_ion_path", abi.parse_ion_path(
        os.fsencode(file_path),
        c_size_t(max_cache_size),
        ctypes.byref(ptr),
    ))
    return SampleFile(ptr, abi)


LARGEST_GAP = 131072


def gap_for(total: int) -> int:
    return min(LARGEST_GAP, total // 8)


class _RemoteSource:
    HEADER_BYTES = 1024

    def __init__(self, url: str) -> None:
        self.url = url
        self.cache: dict[tuple[int, int], bytes] = {}
        self.bytes_fetched = 0
        self.total = 0
        header = self.fetch(0, self.HEADER_BYTES)
        self.cache[(0, len(header))] = header
        planned = _plan_open(header)
        self.total = max((offset + length for offset, length in planned), default=0)
        self.prefetch(planned)

    def fetch(self, offset: int, length: int) -> bytes:
        from urllib.request import Request, urlopen

        last = offset + length - 1
        request = Request(self.url, headers={"Range": f"bytes={offset}-{last}"})
        with urlopen(request) as response:
            payload = response.read()
        self.bytes_fetched += len(payload)
        return payload

    def read_from_cache(self, offset: int, length: int) -> Optional[bytes]:
        exact = self.cache.get((offset, length))
        if exact is not None:
            return exact
        for (start, size), payload in self.cache.items():
            if offset >= start and offset + length <= start + size:
                begin = offset - start
                return payload[begin:begin + length]
        return None

    def read(self, offset: int, length: int) -> bytes:
        held = self.read_from_cache(offset, length)
        if held is not None:
            return held
        payload = self.fetch(offset, length)
        self.cache[(offset, length)] = payload
        return payload

    def prefetch(self, ranges: list[tuple[int, int]]) -> None:
        missing = [
            (offset, length) for offset, length in ranges
            if length > 0 and self.read_from_cache(offset, length) is None
        ]
        for offset, length in _merge_ranges(missing, gap_for(self.total)):
            if self.read_from_cache(offset, length) is None:
                self.cache[(offset, length)] = self.fetch(offset, length)


def _merge_ranges(
    ranges: list[tuple[int, int]], gap: int = LARGEST_GAP
) -> list[tuple[int, int]]:
    wanted = sorted((offset, length) for offset, length in ranges if length)
    merged: list[list[int]] = []
    for offset, length in wanted:
        if merged and offset - (merged[-1][0] + merged[-1][1]) <= gap:
            end = max(merged[-1][0] + merged[-1][1], offset + length)
            merged[-1][1] = end - merged[-1][0]
        else:
            merged.append([offset, length])
    return [(offset, length) for offset, length in merged]


def _plan_open(header: bytes) -> list[tuple[int, int]]:
    abi = _get_abi()
    buffer = _Buf()
    array = (c_uint8 * len(header)).from_buffer_copy(header)
    _check("plan_open", abi.plan_open(
        ctypes.cast(array, POINTER(c_uint8)),
        c_size_t(len(header)),
        ctypes.byref(buffer),
    ))
    packed = buf_to_bytes(abi, buffer)
    ranges = []
    for index in range(0, len(packed), 16):
        offset = int.from_bytes(packed[index:index + 8], "little")
        length = int.from_bytes(packed[index + 8:index + 16], "little")
        ranges.append((offset, length))
    return ranges


def _plan_eic(file: SampleFile, target_mz: float, from_rt: float, to_rt: float,
              ppm_tol: float, mz_tol: float) -> list[tuple[int, int]]:
    buffer = _Buf()
    _check("plan_eic", file._abi.plan_eic(
        file.handle,
        c_double(target_mz), c_double(from_rt), c_double(to_rt),
        c_double(ppm_tol), c_double(mz_tol),
        ctypes.byref(buffer),
    ))
    packed = buf_to_bytes(file._abi, buffer)
    return [
        (int.from_bytes(packed[i:i + 8], "little"),
         int.from_bytes(packed[i + 8:i + 16], "little"))
        for i in range(0, len(packed), 16)
    ]


def _plan_scans(file: SampleFile, query_type: int, from_value: float,
                to_value: float, level: int) -> list[tuple[int, int]]:
    buffer = _Buf()
    _check("plan_scans", file._abi.plan_scans(
        file.handle,
        c_uint8(query_type), c_double(from_value), c_double(to_value),
        c_uint8(level),
        ctypes.byref(buffer),
    ))
    packed = buf_to_bytes(file._abi, buffer)
    return [
        (int.from_bytes(packed[i:i + 8], "little"),
         int.from_bytes(packed[i + 8:i + 16], "little"))
        for i in range(0, len(packed), 16)
    ]


def parse_ion_remote(url: str, max_cache_size: int = 0) -> SampleFile:
    """Load an ion file over HTTP, reading only the bytes it needs.

    Fetches the header, asks the library for the file's map, prefetches those
    ranges in one pass, then serves everything else on demand.

    Args:
        url: HTTP or HTTPS URL pointing to an ion file.
        max_cache_size: Cache size in bytes. 0 sets no limit.

    Returns: A SampleFile for use in other quantion functions.
    """
    if not isinstance(url, str) or not url:
        raise ValueError("parse_ion_remote: url must be a non-empty string")
    if not (url.startswith("http://") or url.startswith("https://")):
        raise ValueError("parse_ion_remote: url must start with http:// or https://")
    if not isinstance(max_cache_size, int) or max_cache_size < 0:
        raise ValueError("parse_ion_remote: max_cache_size must be a non-negative integer")

    source = _RemoteSource(url)

    def serve(context, offset, length, dest):
        try:
            payload = source.read(int(offset), int(length))
        except Exception:
            return -1
        if len(payload) != length:
            return -1
        ctypes.memmove(dest, payload, length)
        return 0

    callback = READ_RANGE(serve)
    abi = _get_abi()
    ptr = ctypes.c_void_p()
    _check("parse_ion_source", abi.parse_ion_source(
        callback,
        None,
        c_size_t(max_cache_size),
        ctypes.byref(ptr),
    ))
    file = SampleFile(ptr, abi)
    file._range_callback = callback
    file._remote_source = source
    return file


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


def mzml_to_ion_file(
    input_path: str,
    output_path: str,
    level: int = 12,
    f32_compress: bool = False,
) -> None:
    """Stream-convert an mzML file on disk into an ion file on disk.

    Reads the mzML file one scan at a time and writes the ion file directly,
    so the whole file is never held in memory. Use this for large files.

    Args:
        input_path: Path to the source mzML file.
        output_path: Path where the ion file is written.
        level: Compression level, integer 0 (none) to 22 (max). Default 12.
        f32_compress: If True, store intensity values as 32-bit float.
    """
    abi = _get_abi()
    if not input_path:
        raise ValueError("mzml_to_ion_file: input_path is required")
    if not output_path:
        raise ValueError("mzml_to_ion_file: output_path is required")
    if not isinstance(level, int) or not (0 <= level <= 22):
        raise ValueError("mzml_to_ion_file: level must be an integer in [0, 22]")
    if not isinstance(f32_compress, bool):
        raise TypeError("mzml_to_ion_file: f32_compress must be a bool")
    _check(
        "convert_mzml_file_to_ion_file",
        abi.convert_mzml_file_to_ion_file(
            input_path.encode("utf-8"),
            output_path.encode("utf-8"),
            c_uint8(level),
            c_uint8(1 if f32_compress else 0),
        ),
    )


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
    source = getattr(file, "_remote_source", None)
    if source is not None:
        source.prefetch(_plan_eic(file, target_mz, from_rt, to_rt, ppm_tol, mz_tol))
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
    return buf_to_records(abi, buf)


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
    rows = buf_to_records(abi, buf)
    return rows[0] if rows else {}


def _shape_code(shape: str) -> int:
    name = str(shape).lower()
    if name == "gaussian":
        return 0
    if name == "emg":
        return 1
    raise ValueError(f"shape must be 'gaussian' or 'emg', got {shape!r}")


def fit_peak(
    x: Sequence,
    y: Sequence,
    rt: float,
    intensity: float,
    shape: str = "emg",
) -> Optional[Dict]:
    """Fit a Gaussian or EMG model to one peak.

    The apex (rt, intensity) from peak picking seeds the fit, so the optimizer
    only has to find the shape.

    Args:
        x: Retention-time sequence. Minimum 5 points.
        y: Intensity sequence, same length as x.
        rt: Apex retention time.
        intensity: Apex intensity.
        shape: 'gaussian' or 'emg'. Default 'emg'.

    Returns:
        Dict with keys shape, height, center, fwhm, tail, r2.
        None if the peak could not be fit.
    """
    abi = _get_abi()
    x_array, y_array = _to_f64_array(x), _to_f64_array(y)
    if len(x_array) != len(y_array) or len(x_array) < 5:
        raise ValueError("x and y must have equal length >= 5")
    buf = _Buf()
    _check("fit_peak", abi.fit_peak(
        _as_f64_ptr(x_array),
        _as_f64_ptr(y_array),
        c_size_t(len(x_array)),
        c_double(rt),
        c_double(intensity),
        c_int32(_shape_code(shape)),
        ctypes.byref(buf),
    ))
    rows = buf_to_records(abi, buf)
    return rows[0] if rows else None


def draw_peak(x: Sequence, params: Dict) -> np.ndarray:
    """Render a fitted peak model over x.

    Args:
        x: Retention-time sequence to draw on.
        params: Result of fit_peak (shape, height, center, fwhm, tail).

    Returns:
        Intensity array with the same length and x as the input.
    """
    abi = _get_abi()
    x_array = _to_f64_array(x)
    buf = _Buf()
    _check("draw_peak", abi.draw_peak(
        _as_f64_ptr(x_array),
        c_size_t(len(x_array)),
        c_int32(_shape_code(params["shape"])),
        c_double(params["height"]),
        c_double(params["center"]),
        c_double(params["fwhm"]),
        c_double(params.get("tail", 0.0)),
        ctypes.byref(buf),
    ))
    return buf_to_f64(abi, buf)


def find_noise_level(y: Sequence) -> Dict[str, Any]:
    """Estimate the noise level of a signal.

    Args:
        y: Intensity sequence.

    Returns:
        A dict with `width` (int) and `intensity` (float).
    """
    abi = _get_abi()
    arr = np.ascontiguousarray(y, dtype=np.float32)
    ptr = arr.ctypes.data_as(POINTER(c_float))
    width = c_size_t(0)
    intensity = c_double(0.0)
    _check("find_noise_level", abi.find_noise_level(
        ptr, len(arr), ctypes.byref(width), ctypes.byref(intensity),
    ))
    return {"width": int(width.value), "intensity": float(intensity.value)}


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

    source = getattr(file, "_remote_source", None)
    if source is not None and abi.plan_scans is not None:
        source.prefetch(_plan_scans(file, query_type, a, b, level))

    buf = _Buf()
    _check("get_scans", abi.get_scans(
        file.handle,
        c_uint8(query_type), c_double(a), c_double(b), c_uint8(level),
        ctypes.byref(buf),
    ))
    scans = buf_to_scans(abi, buf)
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
    return buf_to_records(abi, buf)


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
        raw_idx = item.get("idx", item.get("index"))
        if raw_idx is None or int(raw_idx) < 0:
            raise ValueError(f"get_peaks_from_chrom: item {i} needs a non-negative integer idx")
        indices[i] = int(raw_idx)
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
    return buf_to_records(abi, buf)


def find_features(
    file: SampleFile,
    from_rt: float,
    to_rt: float,
    *,
    eic: Optional[Dict[str, float]] = None,
    grid: Optional[Dict[str, float]] = None,
    options: Optional[Dict[str, Any]] = None,
    cores: int = 1,
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
        "min_r2":               _FindPeakDefaults.min_r2,
        "shape":                _FindPeakDefaults.shape,
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
        ctypes.byref(buf),
    ))
    return buf_to_records(abi, buf)


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
        "min_r2":               _FindPeakDefaults.min_r2,
        "shape":                _FindPeakDefaults.shape,
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
    return buf_to_records(abi, buf)


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
        "min_r2":               _FindPeakDefaults.min_r2,
        "shape":                _FindPeakDefaults.shape,
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
        ctypes.byref(buf),
    ))
    return buf_to_records(abi, buf)
