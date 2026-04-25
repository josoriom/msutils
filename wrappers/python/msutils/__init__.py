from __future__ import annotations

from msutils._bridge import _ABI, load_library
from msutils import api as _api
from msutils.api import (
    parse_mzml,
    parse_bin,
    bin_to_json,
    convert_bin_to_mzml,
    mzml_to_bin,
    calculate_eic,
    get_scans,
    find_peaks,
    get_peak,
    find_noise_level,
    calculate_baseline,
    get_peaks_from_eic,
    get_peaks_from_chrom,
    find_features,
    find_feature,
    get_features,
)
from msutils.file import MzMlFile
from msutils._pack import pack_peak_options, encode_target_ids, unpack_targets
from msutils._shared import to_cores

__all__ = [
    "MsUtils",
    "parse_mzml",
    "parse_bin",
    "bin_to_json",
    "convert_bin_to_mzml",
    "mzml_to_bin",
    "calculate_eic",
    "get_scans",
    "find_peaks",
    "get_peak",
    "find_noise_level",
    "calculate_baseline",
    "get_peaks_from_eic",
    "get_peaks_from_chrom",
    "find_features",
    "find_feature",
    "get_features",
    "MzMlFile",
    "pack_peak_options",
    "encode_target_ids",
    "unpack_targets",
    "to_cores",
]


class MsUtils:
    def __init__(self, abi: _ABI) -> None:
        self._abi = abi
        _api._set_abi(abi)

    def parse_mzml(self, data: bytes) -> MzMlFile:
        return parse_mzml(data)

    def parse_bin(self, data: bytes) -> MzMlFile:
        return parse_bin(data)

    def bin_to_json(self, file: MzMlFile):
        return bin_to_json(file)

    def convert_bin_to_mzml(self, file: MzMlFile) -> str:
        return convert_bin_to_mzml(file)

    def mzml_to_bin(self, file: MzMlFile, level: int = 12, f32_compress: bool = False) -> bytes:
        return mzml_to_bin(file, level, f32_compress)

    def calculate_eic(self, file: MzMlFile, target_mz: float, from_rt: float, to_rt: float, ppm_tol: float = 20.0, mz_tol: float = 0.005):
        return calculate_eic(file, target_mz, from_rt, to_rt, ppm_tol, mz_tol)

    def get_scans(self, file: MzMlFile, *, rt_range=None, rt=None, mz_range=None, mz=None, level: int = 1):
        return get_scans(file, rt_range=rt_range, rt=rt, mz_range=mz_range, mz=mz, level=level)

    def find_peaks(self, x, y, options=None):
        return find_peaks(x, y, options)

    def get_peak(self, x, y, rt: float, range_: float, options=None):
        return get_peak(x, y, rt, range_, options)

    def find_noise_level(self, y) -> float:
        return find_noise_level(y)

    def calculate_baseline(self, y, baseline_window: int = 0, baseline_window_factor: int = 0):
        return calculate_baseline(y, baseline_window, baseline_window_factor)

    def get_peaks_from_eic(self, file: MzMlFile, targets, from_rt: float = 0.5, to_rt: float = 5.0, options=None, cores: int = 1):
        return get_peaks_from_eic(file, targets, from_rt, to_rt, options, cores)

    def get_peaks_from_chrom(self, file: MzMlFile, items, options=None, cores: int = 1):
        return get_peaks_from_chrom(file, items, options, cores)

    def find_features(self, file: MzMlFile, from_rt: float, to_rt: float, *, eic=None, grid=None, options=None, cores: int = 1):
        return find_features(file, from_rt, to_rt, eic=eic, grid=grid, options=options, cores=cores)

    def find_feature(self, file: MzMlFile, targets, *, scan_eic=None, eic=None, options=None, cores: int = 1):
        return find_feature(file, targets, scan_eic=scan_eic, eic=eic, options=options, cores=cores)

    def get_features(self, directory_path: str, from_rt: float, to_rt: float, *, eic=None, grid=None, grouping=None, options=None, cores: int = 1):
        return get_features(directory_path, from_rt, to_rt, eic=eic, grid=grid, grouping=grouping, options=options, cores=cores)

    def __repr__(self) -> str:
        return "<MsUtils ready>"


_, _abi = load_library()
_api._set_abi(_abi)