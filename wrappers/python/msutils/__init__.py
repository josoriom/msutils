from __future__ import annotations

from msutils._bridge import _ABI, load_library
from msutils import api as _api
from msutils.api import (
    parse_mzml,
    parse_ion,
    parse_ion_url,
    ion_to_json,
    ion_to_mzml,
    mzml_to_ion,
    mzml_to_ion_file,
    calculate_eic,
    get_scans,
    find_peaks,
    get_peak,
    fit_peak,
    draw_peak,
    find_noise_level,
    calculate_baseline,
    get_peaks_from_eic,
    get_peaks_from_chrom,
    find_features,
    find_feature,
    get_features,
)
from msutils.file import SampleFile
from msutils._pack import pack_peak_options, encode_target_ids, unpack_targets
from msutils._shared import to_cores

__all__ = [
    "MsUtils",
    "parse_mzml",
    "parse_ion",
    "parse_ion_url",
    "ion_to_json",
    "ion_to_mzml",
    "mzml_to_ion",
    "mzml_to_ion_file",
    "calculate_eic",
    "get_scans",
    "find_peaks",
    "get_peak",
    "fit_peak",
    "draw_peak",
    "find_noise_level",
    "calculate_baseline",
    "get_peaks_from_eic",
    "get_peaks_from_chrom",
    "find_features",
    "find_feature",
    "get_features",
    "SampleFile",
    "pack_peak_options",
    "encode_target_ids",
    "unpack_targets",
    "to_cores",
]


class MsUtils:
    def __init__(self, abi: _ABI) -> None:
        self._abi = abi
        _api._set_abi(abi)

    def parse_mzml(self, data: bytes) -> SampleFile:
        return parse_mzml(data)

    def parse_ion(self, data: bytes) -> SampleFile:
        return parse_ion(data)

    def parse_ion_url(self, url: str, max_cache_size: int = 0) -> SampleFile:
        return parse_ion_url(url, max_cache_size)

    def ion_to_json(self, file: SampleFile):
        return ion_to_json(file)

    def ion_to_mzml(self, file: SampleFile) -> str:
        return ion_to_mzml(file)

    def mzml_to_ion(self, file: SampleFile, level: int = 12, f32_compress: bool = False) -> bytes:
        return mzml_to_ion(file, level, f32_compress)

    def mzml_to_ion_file(self, input_path: str, output_path: str, level: int = 12, f32_compress: bool = False, section_on_disk: bool = False) -> None:
        return mzml_to_ion_file(input_path, output_path, level, f32_compress, section_on_disk)

    def calculate_eic(self, file: SampleFile, target_mz: float, from_rt: float, to_rt: float, ppm_tol: float = 20.0, mz_tol: float = 0.005):
        return calculate_eic(file, target_mz, from_rt, to_rt, ppm_tol, mz_tol)

    def get_scans(self, file: SampleFile, *, rt_range=None, rt=None, mz_range=None, mz=None, level: int = 1):
        return get_scans(file, rt_range=rt_range, rt=rt, mz_range=mz_range, mz=mz, level=level)

    def find_peaks(self, x, y, options=None):
        return find_peaks(x, y, options)

    def get_peak(self, x, y, rt: float, range_: float, options=None):
        return get_peak(x, y, rt, range_, options)

    def fit_peak(self, x, y, rt: float, intensity: float, shape: str = "emg"):
        return fit_peak(x, y, rt, intensity, shape)

    def draw_peak(self, x, params):
        return draw_peak(x, params)

    def find_noise_level(self, y) -> float:
        return find_noise_level(y)

    def calculate_baseline(self, y, lambda_: int = 0, max_iterations: int = 0):
        return calculate_baseline(y, lambda_, max_iterations)

    def get_peaks_from_eic(self, file: SampleFile, targets, from_rt: float = 0.5, to_rt: float = 5.0, options=None, cores: int = 1):
        return get_peaks_from_eic(file, targets, from_rt, to_rt, options, cores)

    def get_peaks_from_chrom(self, file: SampleFile, items, options=None, cores: int = 1):
        return get_peaks_from_chrom(file, items, options, cores)

    def find_features(self, file: SampleFile, from_rt: float, to_rt: float, *, eic=None, grid=None, options=None, cores: int = 1):
        return find_features(file, from_rt, to_rt, eic=eic, grid=grid, options=options, cores=cores)

    def find_feature(self, file: SampleFile, targets, *, scan_eic=None, eic=None, options=None, cores: int = 1):
        return find_feature(file, targets, scan_eic=scan_eic, eic=eic, options=options, cores=cores)

    def get_features(self, directory_path: str, from_rt: float, to_rt: float, *, eic=None, grid=None, grouping=None, options=None, cores: int = 1):
        return get_features(directory_path, from_rt, to_rt, eic=eic, grid=grid, grouping=grouping, options=options, cores=cores)

    def __repr__(self) -> str:
        return "<MsUtils ready>"


_, _abi = load_library()
_api._set_abi(_abi)