from typing import Any

def to_cores(cores: Any) -> int:
    try:
        v = int(cores)
    except (TypeError, ValueError):
        return 1
    return max(1, v)