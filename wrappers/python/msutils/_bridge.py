from __future__ import annotations

import ctypes
import json
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
import threading

import numpy as np


class _Buf(ctypes.Structure):
    _fields_ = [("ptr", POINTER(c_uint8)), ("len", c_size_t)]


class PeakOptions(ctypes.Structure):
    _fields_ = [
        ("integral_threshold",     c_double),
        ("intensity_threshold",    c_double),
        ("width_threshold",        c_int32),
        ("_pad",                   c_int32),
        ("noise",                  c_double),
        ("auto_noise",             c_int32),
        ("auto_baseline",          c_int32),
        ("baseline_window",        c_int32),
        ("baseline_window_factor", c_int32),
        ("allow_overlap",          c_int32),
        ("window_size",            c_int32),
        ("sn_ratio",               c_double),
    ]


assert ctypes.sizeof(PeakOptions) == 64, (
    f"PeakOptions is {ctypes.sizeof(PeakOptions)} bytes — expected 64"
)

_CODE_MSG = {
    1: "invalid arguments",
    2: "panic inside Rust",
    4: "parse error",
    5: "encode error",
    6: "not supported for this file type",
}


def _check(name: str, code: int) -> None:
    if code != 0:
        msg = _CODE_MSG.get(code, f"unknown error code {code}")
        raise RuntimeError(f"msutils/{name}: {msg} (code={code})")


class _ABI:
    _REQUIRED: dict[str, tuple[list, object]] = {
        "parse_mzml": (
            [POINTER(c_uint8), c_size_t, POINTER(c_void_p)], c_int32),
        "parse_bin": ([POINTER(c_uint8), c_size_t, c_size_t, POINTER(c_void_p)], c_int32),
        "free_mzml": ([c_void_p], None),
        "free_": ([POINTER(c_uint8), c_size_t], None),
        "bin_to_json": ([c_void_p, POINTER(_Buf)], c_int32),
        "bin_to_mzml": ([c_void_p, POINTER(_Buf)], c_int32),
        "mzml_to_bin": ([c_void_p, POINTER(_Buf), c_uint8, c_uint8], c_int32),
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
            POINTER(PeakOptions), c_int32, c_int32, c_int32, POINTER(_Buf)], c_int32),
        "get_features": ([ctypes.c_char_p,
            c_double, c_double, c_double, c_double,
            c_double, c_double, c_double,
            c_double, c_double, c_double,
            c_int32, POINTER(PeakOptions), c_int32, c_int32, c_int32, POINTER(_Buf)], c_int32),
        "find_noise_level": ([POINTER(c_float), c_size_t], c_float),
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
                    f"msutils: required symbol '{name}' not found in shared library"
                )


def _read_and_free(abi: _ABI, buf: _Buf) -> bytes:
    try:
        return bytes(buf.ptr[:buf.len])
    finally:
        if buf.ptr:
            abi.free_(buf.ptr, buf.len)
            buf.ptr = None
            buf.len = 0


def buf_to_json(abi: _ABI, buf: _Buf):
    return json.loads(_read_and_free(abi, buf).decode("utf-8"))


def buf_to_str(abi: _ABI, buf: _Buf) -> str:
    return _read_and_free(abi, buf).decode("utf-8")


def buf_to_bytes(abi: _ABI, buf: _Buf) -> bytes:
    return _read_and_free(abi, buf)


def buf_to_f64(abi: _ABI, buf: _Buf) -> np.ndarray:
    raw = _read_and_free(abi, buf)
    return np.frombuffer(raw, dtype=np.float64).copy()


def _default_lib_name() -> str:
    names = {"linux": "libmsutils.so", "darwin": "libmsutils.dylib", "win32": "msutils.dll"}
    return names.get(sys.platform, "libmsutils.so")


def _bundled_lib_path() -> Optional[Path]:
    here = Path(__file__).parent
    machine = platform.machine().lower()

    if sys.platform == "linux":
        if machine in ("x86_64", "amd64"):
            candidate = here / "native" / "linux-x86_64" / "libmsutils.so"
            return candidate if candidate.exists() else None
        if machine in ("aarch64", "arm64"):
            candidate = here / "native" / "linux-arm64" / "libmsutils.so"
            return candidate if candidate.exists() else None

    if sys.platform == "darwin":
        if machine in ("x86_64", "amd64"):
            candidate = here / "native" / "macos-x86_64" / "libmsutils.dylib"
            return candidate if candidate.exists() else None
        if machine in ("arm64", "aarch64"):
            candidate = here / "native" / "macos-arm64" / "libmsutils.dylib"
            return candidate if candidate.exists() else None

    if sys.platform == "win32":
        candidate = here / _default_lib_name()
        return candidate if candidate.exists() else None

    return None


def _start_stderr_capture():
    r_fd, w_fd = os.pipe()
    os.dup2(w_fd, 2)
    os.close(w_fd)
    def reader():
        with os.fdopen(r_fd, 'r', errors='replace') as f:
            for line in f:
                print(line, end='', file=sys.stdout, flush=True)
    threading.Thread(target=reader, daemon=True).start()

def load_library(path: Optional[str] = None) -> tuple[ctypes.CDLL, _ABI]:
    _start_stderr_capture()
    if path is None:
        path = os.environ.get("MSUTILS_LIB")
    if path is None:
        candidate = _bundled_lib_path()
        if candidate is not None:
            path = str(candidate)
    if path is None:
        raise FileNotFoundError(
            "Cannot find the msutils shared library. "
            "Set the MSUTILS_LIB environment variable to the path of the .so/.dylib/.dll."
        )
    try:
        lib = ctypes.CDLL(str(path))
    except OSError as exc:
        raise OSError(f"Failed to dlopen '{path}': {exc}") from exc
    return lib, _ABI(lib)