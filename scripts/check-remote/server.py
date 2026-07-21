"""A local HTTP server that supports Range requests and counts what it sends.

Used by check-remote to prove a wrapper fetches only the bytes it needs.
"""

import http.server
import socketserver
import threading

HEADER_BYTES = 1024
TINY_BYTES = 512


def parse_range(header: str, total: int):
    if not header.startswith("bytes=") or "," in header:
        return None
    first_text, _, last_text = header[len("bytes="):].partition("-")
    if not first_text and not last_text:
        return None
    if not first_text:
        length = int(last_text)
        if length <= 0:
            return None
        return max(0, total - length), total - 1
    first = int(first_text)
    last = int(last_text) if last_text else total - 1
    if first >= total or last < first:
        return None
    return first, min(last, total - 1)


class RangeServer:
    def __init__(self, path: str) -> None:
        self.data = open(path, "rb").read()
        self.bytes_sent = 0
        self.requests = 0
        self._server = socketserver.TCPServer(("127.0.0.1", 0), self._handler())
        self.port = self._server.server_address[1]
        self._thread = threading.Thread(target=self._server.serve_forever, daemon=True)
        self._thread.start()

    @property
    def url(self) -> str:
        return f"http://127.0.0.1:{self.port}/api.ion"

    def reset(self) -> None:
        self.bytes_sent = 0
        self.requests = 0

    def stop(self) -> None:
        self._server.shutdown()
        self._server.server_close()

    def _handler(self):
        outer = self

        class Handler(http.server.BaseHTTPRequestHandler):
            def do_GET(self):
                if self.path == "/stats":
                    self.reply_stats()
                    return
                if self.path == "/reset":
                    outer.reset()
                    self.reply_stats()
                    return
                if self.path == "/star":
                    self.reply_without_total()
                    return
                if self.path == "/tiny":
                    self.reply_range(outer.data[:TINY_BYTES])
                    return
                if self.path == "/nolength":
                    self.reply_ignoring_range()
                    return
                self.reply_range(outer.data)

            def do_HEAD(self):
                self.send_response(200)
                self.send_header("Accept-Ranges", "bytes")
                self.send_header("Content-Length", str(len(outer.data)))
                self.end_headers()

            def reply_range(self, data):
                outer.requests += 1
                total = len(data)
                wanted = self.headers.get("Range")
                if wanted:
                    span = parse_range(wanted, total)
                    if span is None:
                        self.send_response(416)
                        self.send_header("Content-Range", f"bytes */{total}")
                        self.send_header("Content-Length", "0")
                        self.end_headers()
                        return
                    first, last = span
                    body = data[first:last + 1]
                    self.send_response(206)
                    self.send_header("Accept-Ranges", "bytes")
                    self.send_header("Content-Range", f"bytes {first}-{last}/{total}")
                else:
                    body = data
                    self.send_response(200)
                    self.send_header("Accept-Ranges", "bytes")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
                outer.bytes_sent += len(body)

            def reply_without_total(self):
                outer.requests += 1
                body = outer.data[:HEADER_BYTES]
                self.send_response(206)
                self.send_header("Accept-Ranges", "bytes")
                self.send_header("Content-Range", f"bytes 0-{HEADER_BYTES - 1}/*")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
                outer.bytes_sent += len(body)

            def reply_ignoring_range(self):
                outer.requests += 1
                body = outer.data
                self.send_response(200)
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
                outer.bytes_sent += len(body)

            def reply_stats(self):
                body = (
                    f'{{"bytes_sent": {outer.bytes_sent},'
                    f' "requests": {outer.requests}}}'
                ).encode()
                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)

            def log_message(self, *args):
                pass

        return Handler
