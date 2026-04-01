import ctypes
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from msutils._bridge import _ABI


class MzMlFile:
    def __init__(self, handle: ctypes.c_void_p, abi: "_ABI") -> None:
        self._handle: ctypes.c_void_p | None = handle
        self._abi = abi

    @property
    def handle(self) -> ctypes.c_void_p:
        if self._handle is None:
            raise RuntimeError("MzMlFile: file has been disposed")
        return self._handle

    def dispose(self) -> None:
        if self._handle is not None:
            self._abi.free_mzml(self._handle)
            self._handle = None

    def __del__(self) -> None:
        self.dispose()

    def __enter__(self) -> "MzMlFile":
        return self

    def __exit__(self, *_: Any) -> None:
        self.dispose()

    def to_json(self) -> Any:
        from msutils import api as _api
        return _api._bin_to_json_raw(self)

    def to_mzml(self) -> str:
        from msutils import api as _api
        return _api._bin_to_mzml_raw(self)

    def to_bin(self, level: int = 12, f32_compress: bool = False) -> bytes:
        from msutils import api as _api
        return _api._mzml_to_bin_raw(self, level, f32_compress)

    def __repr__(self) -> str:
        if self._handle is None:
            return "<MzMlFile disposed>"
        v = ctypes.cast(self._handle, ctypes.c_void_p).value
        return f"<MzMlFile handle=0x{v:x}>"