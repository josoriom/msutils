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
        "bin_to_json": ([c_void_p, POINTER(_Buf)], c_int32),
        "bin_to_mzml": ([c_void_p, POINTER(_Buf)], c_int32),
        "mzml_to_bin": ([c_void_p, POINTER(_Buf), c_uint8, c_uint8], c_int32),
        "convert_mzml_file_to_ion_file": (
            [ctypes.c_char_p, ctypes.c_char_p, c_uint8, c_uint8, c_uint8], c_int32),
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