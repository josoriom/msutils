import { test_range_read_internals } from "../utilities/backendWasm";

const { get_range_once, read_range_into_memory, allowWorkerOnly, MAX_READ_TRIES, set_wasm_memory } =
  test_range_read_internals;

describe("Range Read Integration", () => {
  let xhr_sends: number;
  let xhr_requests: Array<{ method: string; url: string; range: string }>;

  beforeEach(() => {
    xhr_sends = 0;
    xhr_requests = [];

    (global as any).XMLHttpRequest = jest.fn(() => ({
      open: jest.fn((method: string, url: string) => {
        xhr_requests.push({ method, url, range: "" });
      }),
      send: jest.fn(() => {
        xhr_sends++;
      }),
      setRequestHeader: jest.fn((header: string, value: string) => {
        if (header === "Range") {
          xhr_requests[xhr_requests.length - 1].range = value;
        }
      }),
      status: 206,
      response: null as any,
      getResponseHeader: jest.fn(() => null),
    }));
  });

  describe("get_range_once: outcome classification", () => {
    test("206 + correct body length + Content-Range → ok", () => {
      const correct_body = new Uint8Array([1, 2, 3, 4, 5]);

      (global as any).XMLHttpRequest = jest.fn(() => ({
        open: jest.fn(),
        send: jest.fn(),
        setRequestHeader: jest.fn(),
        status: 206,
        response: correct_body.buffer,
        getResponseHeader: jest.fn(() => "bytes 0-4/50000"),
      }));

      const outcome = get_range_once("http://example.com/file", 0, 4, 5);
      expect(outcome.kind).toBe("ok");
      if (outcome.kind === "ok") {
        expect(outcome.bytes.length).toBe(5);
      }
    });

    test("206 + short body → fail", () => {
      const short_body = new Uint8Array([1, 2, 3]);

      (global as any).XMLHttpRequest = jest.fn(() => ({
        open: jest.fn(),
        send: jest.fn(),
        setRequestHeader: jest.fn(),
        status: 206,
        response: short_body.buffer,
        getResponseHeader: jest.fn(() => null),
      }));

      const outcome = get_range_once("http://example.com/file", 0, 4, 5);
      expect(outcome.kind).toBe("fail");
    });

    test("206 + overlong body → fail", () => {
      const long_body = new Uint8Array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);

      (global as any).XMLHttpRequest = jest.fn(() => ({
        open: jest.fn(),
        send: jest.fn(),
        setRequestHeader: jest.fn(),
        status: 206,
        response: long_body.buffer,
        getResponseHeader: jest.fn(() => null),
      }));

      const outcome = get_range_once("http://example.com/file", 0, 4, 5);
      expect(outcome.kind).toBe("fail");
    });

    test("200 (range not honored) → fail", () => {
      const body = new Uint8Array(5);

      (global as any).XMLHttpRequest = jest.fn(() => ({
        open: jest.fn(),
        send: jest.fn(),
        setRequestHeader: jest.fn(),
        status: 200,
        response: body.buffer,
        getResponseHeader: jest.fn(() => null),
      }));

      const outcome = get_range_once("http://example.com/file", 0, 4, 5);
      expect(outcome.kind).toBe("fail");
    });

    test("404 → fail", () => {
      (global as any).XMLHttpRequest = jest.fn(() => ({
        open: jest.fn(),
        send: jest.fn(),
        setRequestHeader: jest.fn(),
        status: 404,
        response: new ArrayBuffer(0),
        getResponseHeader: jest.fn(() => null),
      }));

      const outcome = get_range_once("http://example.com/file", 0, 4, 5);
      expect(outcome.kind).toBe("fail");
    });

    test("416 → fail", () => {
      (global as any).XMLHttpRequest = jest.fn(() => ({
        open: jest.fn(),
        send: jest.fn(),
        setRequestHeader: jest.fn(),
        status: 416,
        response: new ArrayBuffer(0),
        getResponseHeader: jest.fn(() => null),
      }));

      const outcome = get_range_once("http://example.com/file", 0, 4, 5);
      expect(outcome.kind).toBe("fail");
    });

    test("503 (server error) → retry", () => {
      (global as any).XMLHttpRequest = jest.fn(() => ({
        open: jest.fn(),
        send: jest.fn(),
        setRequestHeader: jest.fn(),
        status: 503,
        response: new ArrayBuffer(0),
        getResponseHeader: jest.fn(() => null),
      }));

      const outcome = get_range_once("http://example.com/file", 0, 4, 5);
      expect(outcome.kind).toBe("retry");
    });

    test("thrown send() → retry", () => {
      (global as any).XMLHttpRequest = jest.fn(() => ({
        open: jest.fn(),
        send: jest.fn(() => {
          throw new Error("Network error");
        }),
        setRequestHeader: jest.fn(),
        status: 206,
        response: new ArrayBuffer(0),
        getResponseHeader: jest.fn(() => null),
      }));

      const outcome = get_range_once("http://example.com/file", 0, 4, 5);
      expect(outcome.kind).toBe("retry");
    });
  });

  describe("Content-Range validation", () => {
    test("correct Content-Range header passes", () => {
      const body = new Uint8Array(5);

      (global as any).XMLHttpRequest = jest.fn(() => ({
        open: jest.fn(),
        send: jest.fn(),
        setRequestHeader: jest.fn(),
        status: 206,
        response: body.buffer,
        getResponseHeader: jest.fn((header: string) =>
          header === "Content-Range" ? "bytes 0-4/50000" : null,
        ),
      }));

      const outcome = get_range_once("http://example.com/file", 0, 4, 5);
      expect(outcome.kind).toBe("ok");
    });

    test("mismatched Content-Range header fails", () => {
      const body = new Uint8Array(5);

      (global as any).XMLHttpRequest = jest.fn(() => ({
        open: jest.fn(),
        send: jest.fn(),
        setRequestHeader: jest.fn(),
        status: 206,
        response: body.buffer,
        getResponseHeader: jest.fn((header: string) =>
          header === "Content-Range" ? "bytes 100-104/50000" : null,
        ),
      }));

      const outcome = get_range_once("http://example.com/file", 0, 4, 5);
      expect(outcome.kind).toBe("fail");
    });

    test("missing Content-Range header fails for 206", () => {
      const body = new Uint8Array(5);

      (global as any).XMLHttpRequest = jest.fn(() => ({
        open: jest.fn(),
        send: jest.fn(),
        setRequestHeader: jest.fn(),
        status: 206,
        response: body.buffer,
        getResponseHeader: jest.fn(() => null),
      }));

      const outcome = get_range_once("http://example.com/file", 0, 4, 5);
      expect(outcome.kind).toBe("fail");
    });
  });

  describe("unsigned offset reconstruction", () => {
    test("large offset: offset_lo=0xC0000000, offset_hi=1 → Range: bytes=7516192768-...", () => {
      const body = new Uint8Array(5);
      const range_headers: string[] = [];

      (global as any).XMLHttpRequest = jest.fn(() => ({
        open: jest.fn(),
        send: jest.fn(),
        setRequestHeader: jest.fn((header: string, value: string) => {
          if (header === "Range") {
            range_headers.push(value);
          }
        }),
        status: 206,
        response: body.buffer,
        getResponseHeader: jest.fn(() => null),
      }));

      const offset_lo = 0xc0000000;
      const offset_hi = 1;
      const len = 5;

      read_range_into_memory("http://example.com/file", offset_lo, offset_hi, len, 0);

      expect(range_headers[0]).toContain("bytes=7516192768-");
    });
  });

  describe("retry loop behavior", () => {
    test("throw on first attempt, succeed on second → returns ok", () => {
      let attempts = 0;
      const body = new Uint8Array(5);

      (global as any).XMLHttpRequest = jest.fn(() => ({
        open: jest.fn(),
        send: jest.fn(() => {
          attempts++;
          if (attempts === 1) throw new Error("First try fails");
        }),
        setRequestHeader: jest.fn(),
        status: 206,
        response: body.buffer,
        getResponseHeader: jest.fn(() => "bytes 0-4/50000"),
      }));

      const memory = new WebAssembly.Memory({ initial: 2 });
      set_wasm_memory(memory);

      const result = read_range_into_memory("http://example.com/file", 0, 0, 5, 0);
      expect(result).toBe(0);
      expect(attempts).toBe(2);

      set_wasm_memory(null);
    });

    test("fail-fast (404) after exactly 1 request (no retry)", () => {
      let attempts = 0;

      (global as any).XMLHttpRequest = jest.fn(() => ({
        open: jest.fn(),
        send: jest.fn(() => {
          attempts++;
        }),
        setRequestHeader: jest.fn(),
        status: 404,
        response: new ArrayBuffer(0),
        getResponseHeader: jest.fn(() => null),
      }));

      const memory = new WebAssembly.Memory({ initial: 2 });
      set_wasm_memory(memory);

      const result = read_range_into_memory("http://example.com/file", 0, 0, 5, 0);
      expect(result).toBe(-1);
      expect(attempts).toBe(1);

      set_wasm_memory(null);
    });

    test("retry on 5xx up to MAX_READ_TRIES, then fail", () => {
      let attempts = 0;

      (global as any).XMLHttpRequest = jest.fn(() => ({
        open: jest.fn(),
        send: jest.fn(() => {
          attempts++;
        }),
        setRequestHeader: jest.fn(),
        status: 503,
        response: new ArrayBuffer(0),
        getResponseHeader: jest.fn(() => null),
      }));

      const memory = new WebAssembly.Memory({ initial: 2 });
      set_wasm_memory(memory);

      const result = read_range_into_memory("http://example.com/file", 0, 0, 5, 0);
      expect(result).toBe(-1);
      expect(attempts).toBe(MAX_READ_TRIES);

      set_wasm_memory(null);
    });
  });

  describe("memory writing", () => {
    test("len === 0 returns 0 without sending request", () => {
      const memory = new WebAssembly.Memory({ initial: 2 });
      set_wasm_memory(memory);

      const send_fn = jest.fn();
      (global as any).XMLHttpRequest = jest.fn(() => ({
        open: jest.fn(),
        send: send_fn,
        setRequestHeader: jest.fn(),
        status: 206,
        response: new ArrayBuffer(0),
        getResponseHeader: jest.fn(() => null),
      }));

      const result = read_range_into_memory("http://example.com/file", 0, 0, 0, 0);
      expect(result).toBe(0);
      expect(send_fn).not.toHaveBeenCalled();

      set_wasm_memory(null);
    });

    test("successful read writes bytes to memory at dest_ptr", () => {
      const memory = new WebAssembly.Memory({ initial: 2 });
      set_wasm_memory(memory);

      const payload = new Uint8Array([10, 20, 30, 40, 50]);

      (global as any).XMLHttpRequest = jest.fn(() => ({
        open: jest.fn(),
        send: jest.fn(),
        setRequestHeader: jest.fn(),
        status: 206,
        response: payload.buffer,
        getResponseHeader: jest.fn(() => "bytes 0-4/50000"),
      }));

      const dest_ptr = 100;
      const result = read_range_into_memory("http://example.com/file", 0, 0, 5, dest_ptr);

      expect(result).toBe(0);
      const written = new Uint8Array(memory.buffer, dest_ptr, 5);
      expect(written).toEqual(payload);

      set_wasm_memory(null);
    });
  });

  describe("worker guard", () => {
    test("allow_worker_only throws in main thread environment", () => {
      (global as any).window = {};
      (global as any).document = {};

      expect(() => {
        allowWorkerOnly();
      }).toThrow("parseIon(URL) must run inside a Web Worker");

      delete (global as any).window;
      delete (global as any).document;
    });

    test("allowWorkerOnly passes in worker environment", () => {
      delete (global as any).window;
      delete (global as any).document;

      expect(() => {
        allowWorkerOnly();
      }).not.toThrow();
    });
  });
});
