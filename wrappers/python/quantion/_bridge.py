from __future__ import annotations

import ctypes
import struct
import os
import platform
import sys
from ctypes import (
    POINTER,
    c_double,
    c_float,
    c_int32,
    c_size_t,
    c_uint8,
    c_uint32,
    c_void_p,
)
from pathlib import Path
from typing import Optional

import numpy as np


class _Buf(ctypes.Structure):
    _fields_ = [("ptr", POINTER(c_uint8)), ("len", c_size_t)]


class PeakOptions(ctypes.Structure):
    _fields_ = [
        ("min_integral",          c_double),
        ("min_intensity",         c_double),
        ("min_peak_width_points", c_int32),
        ("shape",                 c_int32),
        ("noise",                 c_double),
        ("auto_noise",            c_int32),
        ("auto_baseline",         c_int32),
        ("lambda_",               c_int32),
        ("max_iterations",        c_int32),
        ("allow_overlap",         c_int32),
        ("_pad2",                 c_int32),
        ("min_snr",               c_double),
        ("min_r2",                c_double),
        ("kernel_size",           c_int32),
    ]


assert ctypes.sizeof(PeakOptions) == 80, (
    f"PeakOptions is {ctypes.sizeof(PeakOptions)} bytes — expected 80"
)

QUANTION_ABI_VERSION = 1

_CODE_MSG = {
    1: "invalid arguments",
    2: "panic inside Rust",
    4: "parse error",
    5: "encode error",
    6: "fast EIC path unavailable: this .ion file has no usable spectrum bounds (A3); re-encode it with the current Ionic to use the fast EIC path",
}


def _check(name: str, code: int) -> None:
    if code != 0:
        msg = _CODE_MSG.get(code, f"unknown error code {code}")
        raise RuntimeError(f"quantion/{name}: {msg} (code={code})")


READ_RANGE = ctypes.CFUNCTYPE(
    c_int32, c_void_p, ctypes.c_uint64, ctypes.c_uint64, POINTER(c_uint8)
)


class _ABI:
    _REQUIRED: dict[str, tuple[list, object]] = {
        "quantion_abi_version": ([], ctypes.c_uint32),
        "quantion_sizeof_peak_options": ([], c_size_t),
        "parse_mzml": (
            [POINTER(c_uint8), c_size_t, POINTER(c_void_p)], c_int32),
        "parse_bin": ([POINTER(c_uint8), c_size_t, c_size_t, POINTER(c_void_p)], c_int32),
        "parse_ion_path": ([ctypes.c_char_p, c_size_t, POINTER(c_void_p)], c_int32),
        "parse_ion_source": ([READ_RANGE, c_void_p, c_size_t, POINTER(c_void_p)], c_int32),
        "plan_open": ([POINTER(c_uint8), c_size_t, POINTER(_Buf)], c_int32),
        "plan_eic": ([c_void_p, c_double, c_double, c_double, c_double, c_double,
             POINTER(_Buf)], c_int32),
        "free_mzml": ([c_void_p], None),
        "free_": ([POINTER(c_uint8), c_size_t], None),
        "bin_to_mzml": ([c_void_p, POINTER(_Buf)], c_int32),
        "mzml_to_bin": ([c_void_p, POINTER(_Buf), c_uint8, c_uint8], c_int32),
        "convert_mzml_file_to_ion_file": (
            [ctypes.c_char_p, ctypes.c_char_p, c_uint8, c_uint8], c_int32),
        "get_peak": ([POINTER(c_double), POINTER(c_double), c_size_t,
             c_double, c_double, POINTER(PeakOptions), POINTER(_Buf)], c_int32),
        "calculate_eic": ([c_void_p, c_double, c_double, c_double,
             c_double, c_double, POINTER(_Buf), POINTER(_Buf)], c_int32),
        "get_scans": ([c_void_p, c_uint8, c_double, c_double, c_uint8, POINTER(_Buf)], c_int32),
        "get_peaks_from_eic": ([c_void_p,
             POINTER(c_double), POINTER(c_double), POINTER(c_double),
             POINTER(c_uint32), POINTER(c_uint32),
             POINTER(c_uint8), c_size_t,
             c_size_t, c_double, c_double,
             POINTER(PeakOptions), c_size_t, POINTER(_Buf)], c_int32),
        "get_peaks_from_chrom": ([c_void_p,
             POINTER(c_uint32), POINTER(c_double), POINTER(c_double),
             c_size_t, POINTER(PeakOptions), c_size_t, POINTER(_Buf)], c_int32),
        "find_peaks": ([POINTER(c_double), POINTER(c_double), c_size_t,
             POINTER(PeakOptions), POINTER(_Buf)], c_int32),
        "fit_peak": ([POINTER(c_double), POINTER(c_double), c_size_t,
             c_double, c_double, c_int32, POINTER(_Buf)], c_int32),
        "draw_peak": ([POINTER(c_double), c_size_t, c_int32,
             c_double, c_double, c_double, c_double, POINTER(_Buf)], c_int32),
        "calculate_baseline": ([POINTER(c_double), c_size_t, c_int32, c_int32,
             POINTER(_Buf)], c_int32),
        "find_feature": ([c_void_p,
             POINTER(c_double), POINTER(c_double), POINTER(c_double),
             POINTER(c_uint32), POINTER(c_uint32),
             POINTER(c_uint8), c_size_t,
             c_size_t, c_size_t,
             c_double, c_double, c_double, c_double,
             POINTER(PeakOptions), POINTER(_Buf)], c_int32),
        "find_features": ([c_void_p,
            c_double, c_double, c_double, c_double,
            c_double, c_double, c_double,
            POINTER(PeakOptions), c_int32, POINTER(_Buf)], c_int32),
        "get_features": ([ctypes.c_char_p,
            c_double, c_double, c_double, c_double,
            c_double, c_double, c_double,
            c_double, c_double, c_double,
            c_int32, POINTER(PeakOptions), c_int32, POINTER(_Buf)], c_int32),
        "find_noise_level": ([POINTER(c_float), c_size_t,
             POINTER(c_size_t), POINTER(c_double)], c_int32),
    }

    _OPTIONAL: dict[str, tuple[list, object]] = {
        "plan_scans": ([c_void_p, c_uint8, c_double, c_double, c_uint8,
             POINTER(_Buf)], c_int32),
    }

    def __init__(self, lib: ctypes.CDLL) -> None:
        for name, (argtypes, restype) in self._REQUIRED.items():
            try:
                fn = getattr(lib, name)
                fn.argtypes = argtypes
                fn.restype  = restype
                setattr(self, name, fn)
            except AttributeError:
                raise RuntimeError(
                    f"quantion: required symbol '{name}' not found in shared library"
                )
        for name, (argtypes, restype) in self._OPTIONAL.items():
            try:
                fn = getattr(lib, name)
                fn.argtypes = argtypes
                fn.restype  = restype
                setattr(self, name, fn)
            except AttributeError:
                setattr(self, name, None)


def _read_and_free(abi: _ABI, buf: _Buf) -> bytes:
    try:
        return ctypes.string_at(buf.ptr, buf.len)
    finally:
        if buf.ptr:
            abi.free_(buf.ptr, buf.len)
            buf.ptr = None
            buf.len = 0


def buf_to_str(abi: _ABI, buf: _Buf) -> str:
    return _read_and_free(abi, buf).decode("utf-8")


def buf_to_bytes(abi: _ABI, buf: _Buf) -> bytes:
    return _read_and_free(abi, buf)


def buf_to_f64(abi: _ABI, buf: _Buf) -> np.ndarray:
    raw = _read_and_free(abi, buf)
    return np.frombuffer(raw, dtype=np.float64).copy()


BRIDGE_MAGIC = 0x42544E51
BRIDGE_LAYOUT_VERSION = 1
BRIDGE_HEADER_BYTES = 32
SECTION_ENTRY_BYTES = 24

PAYLOAD_SCANS = 1
PAYLOAD_ION_IMAGE = 2

ELEMENT_F64 = 1
ELEMENT_U32 = 2
ELEMENT_U64 = 3
ELEMENT_U8 = 4

SECTION_POINT_STARTS = 1
SECTION_MZ = 2
SECTION_INTENSITY = 3
SECTION_RT = 4
SECTION_RT_SECONDS = 5
SECTION_BASE_PEAK_MZ = 6
SECTION_SELECTED_ION_MZ = 7
SECTION_BASE_PEAK_INT = 8
SECTION_TOTAL_ION_CURRENT = 9
SECTION_MS_LEVEL = 10
SECTION_POLARITY = 11
SECTION_POSITION_X = 12
SECTION_POSITION_Y = 13
SECTION_POSITION_Z = 14
SECTION_IMAGE_SHAPE = 15
SECTION_IMAGE_DATA = 16
SECTION_IMAGE_COUNTS = 17

_ELEMENT_SIZES = {ELEMENT_F64: 8, ELEMENT_U32: 4, ELEMENT_U64: 8, ELEMENT_U8: 1}

_SECTION_TYPES = {
    SECTION_POINT_STARTS: ELEMENT_U64,
    SECTION_MZ: ELEMENT_F64,
    SECTION_INTENSITY: ELEMENT_F64,
    SECTION_RT: ELEMENT_F64,
    SECTION_RT_SECONDS: ELEMENT_F64,
    SECTION_BASE_PEAK_MZ: ELEMENT_F64,
    SECTION_SELECTED_ION_MZ: ELEMENT_F64,
    SECTION_BASE_PEAK_INT: ELEMENT_F64,
    SECTION_TOTAL_ION_CURRENT: ELEMENT_F64,
    SECTION_MS_LEVEL: ELEMENT_F64,
    SECTION_POLARITY: ELEMENT_F64,
    SECTION_POSITION_X: ELEMENT_F64,
    SECTION_POSITION_Y: ELEMENT_F64,
    SECTION_POSITION_Z: ELEMENT_F64,
    SECTION_IMAGE_SHAPE: ELEMENT_U32,
    SECTION_IMAGE_DATA: ELEMENT_F64,
    SECTION_IMAGE_COUNTS: ELEMENT_U32,
}

_METADATA_SECTIONS = (
    ("rt_seconds", SECTION_RT_SECONDS),
    ("base_peak_mz", SECTION_BASE_PEAK_MZ),
    ("selected_ion_mz", SECTION_SELECTED_ION_MZ),
    ("base_peak_int", SECTION_BASE_PEAK_INT),
    ("total_ion_current", SECTION_TOTAL_ION_CURRENT),
    ("ms_level", SECTION_MS_LEVEL),
    ("polarity", SECTION_POLARITY),
    ("position_x", SECTION_POSITION_X),
    ("position_y", SECTION_POSITION_Y),
    ("position_z", SECTION_POSITION_Z),
)


class Bridge:
    def __init__(self, raw: bytes, payload_kind: int, record_count: int, sections: dict):
        self.raw = raw
        self.payload_kind = payload_kind
        self.record_count = record_count
        self.sections = sections


def _fail(reason: str):
    raise ValueError(f"quantion bridge: {reason}")


def read_bridge(raw: bytes) -> Bridge:
    total = len(raw)
    if total < BRIDGE_HEADER_BYTES:
        _fail("buffer is shorter than the header")

    magic, layout, kind, count, table_offset, total_bytes, records = struct.unpack_from(
        "<IHHIIQQ", raw, 0
    )
    if magic != BRIDGE_MAGIC:
        _fail("magic does not match")
    if layout != BRIDGE_LAYOUT_VERSION:
        _fail(f"layout version {layout} is not supported")
    if total_bytes != total:
        _fail("total bytes does not match the buffer")
    if table_offset != BRIDGE_HEADER_BYTES:
        _fail("section table does not start after the header")
    if count * SECTION_ENTRY_BYTES > total - BRIDGE_HEADER_BYTES:
        _fail("section table does not fit")

    sections = {}
    reach = BRIDGE_HEADER_BYTES + count * SECTION_ENTRY_BYTES
    for index in range(count):
        start = BRIDGE_HEADER_BYTES + index * SECTION_ENTRY_BYTES
        section_id, element_type, offset, length = struct.unpack_from("<IIQQ", raw, start)
        if section_id in sections:
            _fail(f"section {section_id} appears twice")
        if offset > total or length > total - offset:
            _fail(f"section {section_id} runs past the end")
        if offset % 8 != 0:
            _fail(f"section {section_id} is not aligned to eight bytes")
        if offset < reach:
            _fail(f"section {section_id} overlaps the section before it")
        size = _ELEMENT_SIZES.get(element_type)
        if size is None:
            _fail(f"section {section_id} has an unknown element type")
        if length % size != 0:
            _fail(f"section {section_id} is not a whole number of elements")
        expected = _SECTION_TYPES.get(section_id)
        if expected is not None and expected != element_type:
            _fail(f"section {section_id} has the wrong element type")
        reach = offset + length
        sections[section_id] = (element_type, offset, length)

    return Bridge(raw, kind, records, sections)


def expect_kind(bridge: Bridge, wanted: int) -> None:
    if bridge.payload_kind != wanted:
        _fail(f"expected payload kind {wanted} but found {bridge.payload_kind}")


def _take_section(bridge: Bridge, section_id: int):
    section = bridge.sections.get(section_id)
    if section is None:
        _fail(f"section {section_id} is missing")
    return section


def read_section(bridge: Bridge, section_id: int, dtype) -> np.ndarray:
    _, offset, length = _take_section(bridge, section_id)
    return np.frombuffer(bridge.raw, dtype=dtype, count=length // np.dtype(dtype).itemsize,
                         offset=offset)


def read_scans_bridge(raw: bytes) -> list:
    bridge = read_bridge(raw)
    expect_kind(bridge, PAYLOAD_SCANS)

    point_starts = read_section(bridge, SECTION_POINT_STARTS, np.uint64)
    mz = read_section(bridge, SECTION_MZ, np.float64)
    intensity = read_section(bridge, SECTION_INTENSITY, np.float64)
    rt = read_section(bridge, SECTION_RT, np.float64)

    count = bridge.record_count
    if len(point_starts) != count + 1:
        _fail("point starts do not match the scan count")
    if len(intensity) != len(mz):
        _fail("intensity length does not match m/z length")
    if int(point_starts[-1]) != len(mz):
        _fail("point starts do not span the m/z section")
    if np.any(np.diff(point_starts.astype(np.int64)) < 0):
        _fail("point starts do not increase")

    columns = {name: read_section(bridge, section_id, np.float64)
               for name, section_id in _METADATA_SECTIONS}

    scans = []
    for index in range(count):
        start = int(point_starts[index])
        stop = int(point_starts[index + 1])
        scans.append({
            "rt": float(rt[index]),
            "mz": mz[start:stop],
            "intensity": intensity[start:stop],
            "metadata": {name: float(column[index]) for name, column in columns.items()},
        })
    return scans


def buf_to_scans(abi: _ABI, buf: _Buf) -> list:
    return read_scans_bridge(_read_and_free(abi, buf))


PAYLOAD_PEAKS = 3
PAYLOAD_FEATURES = 4
PAYLOAD_CHROM_PEAKS = 5
PAYLOAD_FIT_RESULT = 6
PAYLOAD_EIC_PEAKS = 7
PAYLOAD_CONSENSUS_FEATURES = 8
PAYLOAD_FOUND_FEATURES = 9

_NUMBER, _COUNT, _TEXT = "number", "count", "text"

_RECORD_COLUMNS = {
    PAYLOAD_PEAKS: [
        ("from", 18, _NUMBER), ("to", 19, _NUMBER), ("rt", 20, _NUMBER),
        ("integral", 21, _NUMBER), ("intensity", 22, _NUMBER),
        ("n_points", 23, _COUNT), ("noise", 24, _NUMBER), ("r2", 25, _NUMBER),
    ],
    PAYLOAD_FEATURES: [
        ("mz", 26, _NUMBER), ("rt", 27, _NUMBER), ("from", 28, _NUMBER),
        ("to", 29, _NUMBER), ("intensity", 30, _NUMBER), ("integral", 31, _NUMBER),
        ("n_points", 32, _COUNT), ("noise", 33, _NUMBER),
    ],
    PAYLOAD_CHROM_PEAKS: [
        ("index", 34, _COUNT), ("id", (42, 43), _TEXT), ("ort", 35, _NUMBER),
        ("rt", 36, _NUMBER), ("from", 37, _NUMBER), ("to", 38, _NUMBER),
        ("intensity", 39, _NUMBER), ("integral", 40, _NUMBER),
        ("total_area", 41, _NUMBER), ("timestamp", (44, 45), _TEXT),
    ],
    PAYLOAD_FIT_RESULT: [
        ("shape", 46, _NUMBER), ("height", 47, _NUMBER), ("center", 48, _NUMBER),
        ("fwhm", 49, _NUMBER), ("tail", 50, _NUMBER), ("r2", 51, _NUMBER),
    ],
    PAYLOAD_EIC_PEAKS: [
        ("id", (52, 53), _TEXT), ("mz", 54, _NUMBER), ("ort", 55, _NUMBER),
        ("rt", 56, _NUMBER), ("from", 57, _NUMBER), ("to", 58, _NUMBER),
        ("intensity", 59, _NUMBER), ("integral", 60, _NUMBER), ("noise", 61, _NUMBER),
    ],
    PAYLOAD_CONSENSUS_FEATURES: [
        ("mz", 62, _NUMBER), ("rt", 63, _NUMBER), ("from", 64, _NUMBER),
        ("to", 65, _NUMBER), ("intensity", 66, _NUMBER), ("integral", 67, _NUMBER),
        ("frequency", 68, _NUMBER),
    ],
    PAYLOAD_FOUND_FEATURES: [
        ("id", (69, 70), _TEXT), ("mz", 71, _NUMBER), ("rt", 72, _NUMBER),
        ("from", 73, _NUMBER), ("to", 74, _NUMBER), ("intensity", 75, _NUMBER),
        ("integral", 76, _NUMBER), ("n_points", 77, _COUNT), ("noise", 78, _NUMBER),
    ],
}


def read_records(raw: bytes) -> list:
    bridge = read_bridge(raw)
    columns = _RECORD_COLUMNS.get(bridge.payload_kind)
    if columns is None:
        _fail(f"payload kind {bridge.payload_kind} is not a record table")

    taken = []
    for name, where, kind in columns:
        if kind == _NUMBER:
            taken.append((name, read_section(bridge, where, np.float64), None))
        elif kind == _COUNT:
            taken.append((name, read_section(bridge, where, np.uint32), None))
        else:
            starts = read_section(bridge, where[0], np.uint64)
            _, offset, length = _take_section(bridge, where[1])
            taken.append((name, starts, raw[offset:offset + length]))

    rows = []
    for index in range(bridge.record_count):
        row = {}
        for name, values, text_bytes in taken:
            if text_bytes is None:
                row[name] = values[index].item()
            else:
                row[name] = text_bytes[int(values[index]):int(values[index + 1])].decode("utf-8")
        rows.append(row)
    return rows


def buf_to_records(abi: _ABI, buf: _Buf) -> list:
    return read_records(_read_and_free(abi, buf))


def _platform_dir() -> str:
    machine = platform.machine().lower()
    arm = machine in ("aarch64", "arm64")
    if sys.platform == "linux":
        return "linux-arm64" if arm else "linux-x86_64"
    if sys.platform == "darwin":
        return "macos-arm64" if arm else "macos-x86_64"
    if sys.platform == "win32":
        return "windows-arm64" if arm else "windows-x86_64"
    raise RuntimeError(f"quantion: unsupported platform {sys.platform}/{machine}")


def _default_lib_name() -> str:
    names = {
        "linux": "libquantion.so",
        "darwin": "libquantion.dylib",
        "win32": "libquantion.dll",
    }
    return names.get(sys.platform, "libquantion.so")


def _is_package_root(directory: Path) -> bool:
    manifest = directory / "pyproject.toml"
    if not manifest.is_file():
        return False
    return 'name = "quantion"' in manifest.read_text(encoding="utf-8")


def _artifacts_root() -> Optional[Path]:
    from_env = os.environ.get("QUANTION_ARTIFACTS_ROOT")
    if from_env:
        chosen = Path(from_env)
        return chosen if chosen.is_dir() else None
    for parent in Path(__file__).resolve().parents:
        candidate = parent / "artifacts"
        if candidate.is_dir():
            return candidate
        if _is_package_root(parent):
            return None
    return None


def _version_key(path: Path) -> tuple:
    return tuple(int(part) if part.isdigit() else 0 for part in path.name.split("."))


def _own_version() -> str:
    from quantion import __version__
    return __version__


def _packaged_lib_path() -> Optional[Path]:
    candidate = Path(__file__).resolve().parent / "native" / _platform_dir() / _default_lib_name()
    return candidate if candidate.exists() else None


def _bundled_lib_path() -> Optional[Path]:
    root = _artifacts_root()
    if root is not None:
        wanted = _platform_dir()
        name = _default_lib_name()

        direct = root / wanted / name
        if direct.exists():
            return direct

        own = _own_version()
        if own:
            preferred = root / own / wanted / name
            if preferred.exists():
                return preferred

        versions = sorted(
            (p for p in root.iterdir()
             if p.is_dir() and (p / wanted / name).exists() and (p / "manifest.json").exists()),
            key=_version_key,
            reverse=True,
        )
        if versions:
            return versions[0] / wanted / name

    return _packaged_lib_path()


def load_library(path: Optional[str] = None) -> tuple[ctypes.CDLL, _ABI]:
    if path is None:
        path = os.environ.get("QUANTION_LIB")
    if path is None:
        candidate = _bundled_lib_path()
        if candidate is not None:
            path = str(candidate)
    if path is None:
        raise FileNotFoundError(
            f"quantion: no library for {_platform_dir()}. "
            f"Build one with 'make {_platform_dir()}', "
            "or set QUANTION_LIB to a library file, "
            "or QUANTION_ARTIFACTS_ROOT to the artifacts folder."
        )
    try:
        lib = ctypes.CDLL(str(path))
    except OSError as exc:
        raise OSError(f"Failed to dlopen '{path}': {exc}") from exc
    abi = _ABI(lib)
    actual = abi.quantion_abi_version()
    if actual != QUANTION_ABI_VERSION:
        raise RuntimeError(
            f"quantion ABI version mismatch: library={actual}, "
            f"wrapper={QUANTION_ABI_VERSION}. Rebuild or update the package."
        )
    # Defence-in-depth: verify the PeakOptions struct layout matches the binary.
    # If the Rust struct ever changes without bumping the ABI version, this
    # catches it at import time instead of silently corrupting peak detection.
    expected = ctypes.sizeof(PeakOptions)
    actual_size = abi.quantion_sizeof_peak_options()
    if actual_size != expected:
        raise RuntimeError(
            f"quantion PeakOptions size mismatch: library={actual_size} bytes, "
            f"wrapper={expected} bytes. The native binary and Python wrapper "
            f"are out of sync — reinstall quantion."
        )
    return lib, abi