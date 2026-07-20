from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from quantion._bridge import PeakOptions

SHAPE_CODES = {"gaussian": 0, "emg": 1}
DEFAULT_SHAPE_CODE = 1


def pack_peak_options(opts: Optional[Dict[str, Any]]) -> Optional[PeakOptions]:
    if not opts:
        return None

    def _get_float(key: str) -> float:
        value = opts.get(key)
        if value is None:
            return math.nan
        parsed = float(value)
        return parsed if math.isfinite(parsed) else math.nan

    def _get_int(key: str) -> int:
        value = opts.get(key)
        if value is None:
            return 0
        if isinstance(value, bool):
            return int(value)
        try:
            parsed = float(value)
            return int(parsed) if math.isfinite(parsed) else 0
        except (TypeError, ValueError):
            return 0

    def _get_bool(key: str) -> int:
        value = opts.get(key)
        if value is None:
            return 0
        if isinstance(value, bool):
            return int(value)
        return 1 if value else 0

    def _get_shape(key: str) -> int:
        value = opts.get(key)
        if value is None:
            return DEFAULT_SHAPE_CODE
        if isinstance(value, str):
            return SHAPE_CODES.get(value.strip().lower(), DEFAULT_SHAPE_CODE)
        return _get_int(key)

    options_struct = PeakOptions()
    options_struct.min_integral          = _get_float("min_integral")
    options_struct.min_intensity         = _get_float("min_intensity")
    options_struct.min_peak_width_points = _get_int("min_peak_width_points")
    options_struct.shape                 = _get_shape("shape")
    options_struct.noise                 = _get_float("noise")
    options_struct.auto_noise            = _get_bool("auto_noise")
    options_struct.auto_baseline         = _get_bool("auto_baseline")
    options_struct.lambda_               = _get_int("lambda_")
    options_struct.max_iterations        = _get_int("max_iterations")
    options_struct.allow_overlap         = _get_bool("allow_overlap")
    options_struct._pad2                 = 0
    options_struct.min_snr               = _get_float("min_snr")
    options_struct.min_r2                = _get_float("min_r2")
    options_struct.kernel_size           = _get_int("kernel_size")
    return options_struct


def encode_target_ids(
    ids: Sequence[str],
) -> Tuple[np.ndarray, np.ndarray, bytes]:
    n       = len(ids)
    offsets = np.zeros(n, dtype=np.uint32)
    lengths = np.zeros(n, dtype=np.uint32)
    parts: List[bytes] = []
    cursor  = 0
    for i, s in enumerate(ids):
        b = (s or "").encode("utf-8")
        offsets[i] = cursor
        lengths[i] = len(b)
        parts.append(b)
        cursor += len(b)
    return offsets, lengths, b"".join(parts)


def unpack_targets(
    targets: Sequence[Dict[str, Any]],
    default_range: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[str]]:
    n      = len(targets)
    rts    = np.empty(n, dtype=np.float64)
    mzs    = np.empty(n, dtype=np.float64)
    ranges = np.empty(n, dtype=np.float64)
    ids: List[str] = []

    for i, t in enumerate(targets):
        rts[i]    = float(t.get("rt",    t.get("retention_time", 0.0)))
        mzs[i]    = float(t.get("mz",    t.get("mass",           0.0)))
        ranges[i] = float(t.get("range", t.get("window",         default_range)))
        ids.append(str(t.get("id", t.get("name", ""))))

    return rts, mzs, ranges, ids