import ctypes
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from quantion._bridge import _ABI


class SampleFile:
    def __init__(self, handle: ctypes.c_void_p, abi: "_ABI") -> None:
        self._handle: ctypes.c_void_p | None = handle
        self._abi = abi

    @property
    def handle(self) -> ctypes.c_void_p:
        if self._handle is None:
            raise RuntimeError("SampleFile: file has been disposed")
        return self._handle

    def dispose(self) -> None:
        if self._handle is not None:
            self._abi.free_mzml(self._handle)
            self._handle = None

    def __del__(self) -> None:
        self.dispose()

    def __enter__(self) -> "SampleFile":
        return self

    def __exit__(self, *_: Any) -> None:
        self.dispose()

    def to_mzml(self) -> str:
        from quantion import api as _api
        return _api._ion_to_mzml_raw(self)

    def to_ion(self, level: int = 12, f32_compress: bool = False) -> bytes:
        from quantion import api as _api
        return _api._mzml_to_ion_raw(self, level, f32_compress)

    def __repr__(self) -> str:
        if self._handle is None:
            return "<SampleFile disposed>"
        v = ctypes.cast(self._handle, ctypes.c_void_p).value
        return f"<SampleFile handle=0x{v:x}>"
