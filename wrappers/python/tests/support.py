import http.server
import os
import socketserver
import sys
import threading
from pathlib import Path

WRAPPER_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = WRAPPER_ROOT.parents[1]
ARTIFACTS_ROOT = REPO_ROOT / "artifacts"
ION_FIXTURE = REPO_ROOT / "core" / "tests" / "fixtures" / "api" / "api.ion"

if str(WRAPPER_ROOT) not in sys.path:
    sys.path.insert(0, str(WRAPPER_ROOT))

os.environ.setdefault("QUANTION_ARTIFACTS_ROOT", str(ARTIFACTS_ROOT))

try:
    import quantion
    from quantion import _bridge
    from quantion import api

    library_problem = ""
except Exception as error:
    quantion = None
    _bridge = None
    api = None
    library_problem = (
        f"the quantion shared library did not load: {error}. "
        f"Point QUANTION_LIB at a freshly built library, "
        f"or put one under {ARTIFACTS_ROOT}."
    )

fixture_problem = (
    "" if ION_FIXTURE.is_file() else f"the ion fixture is missing: {ION_FIXTURE}"
)


class RangeServer:
    def __init__(self, path) -> None:
        self.data = Path(path).read_bytes()
        self.bytes_sent = 0
        self.requests = 0
        self.server = socketserver.TCPServer(("127.0.0.1", 0), self.build_handler())
        self.port = self.server.server_address[1]
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()

    @property
    def url(self) -> str:
        return f"http://127.0.0.1:{self.port}/api.ion"

    def reset(self) -> None:
        self.bytes_sent = 0
        self.requests = 0

    def stop(self) -> None:
        self.server.shutdown()
        self.server.server_close()

    def build_handler(self):
        counter = self

        class Handler(http.server.BaseHTTPRequestHandler):
            def do_GET(self):
                counter.requests += 1
                wanted = self.headers.get("Range")
                if wanted:
                    first, last = wanted.split("=")[1].split("-")
                    body = counter.data[int(first):int(last) + 1]
                    self.send_response(206)
                    self.send_header(
                        "Content-Range",
                        f"bytes {first}-{last}/{len(counter.data)}",
                    )
                else:
                    body = counter.data
                    self.send_response(200)
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
                counter.bytes_sent += len(body)

            def log_message(self, *args):
                pass

        return Handler


def recording_source(total: int):
    class Source(api._RemoteSource):
        def __init__(self) -> None:
            self.url = ""
            self.cache = {}
            self.bytes_fetched = 0
            self.total = total
            self.fetched = []

        def fetch(self, offset: int, length: int) -> bytes:
            self.fetched.append((offset, length))
            self.bytes_fetched += length
            return bytes(length)

    return Source()


def place_library(root: Path, *parts: str) -> Path:
    folder = root.joinpath(*parts, _bridge._platform_dir())
    folder.mkdir(parents=True, exist_ok=True)
    library = folder / _bridge._default_lib_name()
    library.write_bytes(b"")
    return library


def place_manifest(root: Path, version: str) -> Path:
    manifest = root / version / "manifest.json"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text("{}", encoding="utf-8")
    return manifest
