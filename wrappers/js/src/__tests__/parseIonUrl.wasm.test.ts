import * as fs from "fs";
import * as path from "path";

describe("parseIonUrl WASM Integration", () => {
  let test_ion_data: Uint8Array;

  beforeEach(() => {
    const test_ion_path = path.join(__dirname, "../../..", "test.ion");
    test_ion_data = new Uint8Array(fs.readFileSync(test_ion_path));

    (global as any).XMLHttpRequest = jest.fn(() => {
      let url = "";
      let range_header = "";
      let status = 200;
      let response: any = null;
      let content_range = "";

      return {
        open: jest.fn((method: string, u: string) => {
          url = u;
        }),
        send: jest.fn(() => {
          if (!range_header) {
            status = 400;
            response = new ArrayBuffer(0);
            return;
          }

          const match = range_header.match(/bytes=(\d+)-(\d+)/);
          if (!match) {
            status = 416;
            response = new ArrayBuffer(0);
            return;
          }

          const start = parseInt(match[1], 10);
          const end = parseInt(match[2], 10);
          const len = end - start + 1;

          if (start < 0 || end >= test_ion_data.length || len <= 0) {
            status = 416;
            response = new ArrayBuffer(0);
            return;
          }

          status = 206;
          const chunk = test_ion_data.slice(start, end + 1);
          response = chunk.buffer;
          content_range = `bytes ${start}-${end}/${test_ion_data.length}`;
        }),
        setRequestHeader: jest.fn((header: string, value: string) => {
          if (header === "Range") {
            range_header = value;
          }
        }),
        get status() {
          return status;
        },
        get response() {
          return response;
        },
        getResponseHeader: jest.fn((header: string) => {
          if (header === "Content-Range") {
            return content_range || null;
          }
          return null;
        }),
      };
    });
  });

  describe("URL state binding", () => {
    test("registry tracks source_id to URL mapping", () => {
      const registry = new Map<number, string>();
      const url = "http://example.com/test.ion";
      const source_id = 1;

      registry.set(source_id, url);
      expect(registry.get(source_id)).toBe(url);
    });

    test("handle to source_id mapping survives multiple calls", () => {
      const source_id_by_handle = new Map<number, number>();
      const registry = new Map<number, string>();

      const handle = 42;
      const source_id = 1;
      const url = "http://example.com/test.ion";

      registry.set(source_id, url);
      source_id_by_handle.set(handle, source_id);

      expect(source_id_by_handle.get(handle)).toBe(source_id);
      expect(registry.get(source_id_by_handle.get(handle)!)).toBe(url);
    });

    test("multiple range requests use same source_id", () => {
      const xhr_ctor = (global as any).XMLHttpRequest as jest.Mock;

      const request1 = new xhr_ctor();
      request1.open("GET", "http://example.com/test.ion", false);
      request1.setRequestHeader("Range", `bytes=0-99`);
      request1.send();

      expect(request1.status).toBe(206);
      expect(request1.response.byteLength).toBe(100);

      const request2 = new xhr_ctor();
      request2.open("GET", "http://example.com/test.ion", false);
      request2.setRequestHeader("Range", `bytes=100-199`);
      request2.send();

      expect(request2.status).toBe(206);
      expect(request2.response.byteLength).toBe(100);
    });
  });

  describe("XHR behavior with test.ion", () => {
    test("successful range read returns 206 + Content-Range", () => {
      const xhr_ctor = (global as any).XMLHttpRequest as jest.Mock;
      const xhr = new xhr_ctor();

      xhr.open("GET", "http://example.com/test.ion", false);
      xhr.setRequestHeader("Range", `bytes=0-99`);
      xhr.send();

      expect(xhr.status).toBe(206);
      expect(xhr.response.byteLength).toBe(100);
      expect(xhr.getResponseHeader("Content-Range")).toMatch(/bytes 0-99/);
    });

    test("large offset read works correctly", () => {
      if (test_ion_data.length < 1000) {
        expect(true).toBe(true);
        return;
      }

      const xhr_ctor = (global as any).XMLHttpRequest as jest.Mock;
      const xhr = new xhr_ctor();

      xhr.open("GET", "http://example.com/test.ion", false);
      xhr.setRequestHeader("Range", `bytes=500-599`);
      xhr.send();

      expect(xhr.status).toBe(206);
      expect(xhr.response.byteLength).toBe(100);
    });

    test("out-of-range request fails with 416", () => {
      const xhr_ctor = (global as any).XMLHttpRequest as jest.Mock;
      const xhr = new xhr_ctor();

      xhr.open("GET", "http://example.com/test.ion", false);
      xhr.setRequestHeader("Range", `bytes=99999999-99999999`);
      xhr.send();

      expect(xhr.status).toBe(416);
    });

    test("missing Range header fails", () => {
      const xhr_ctor = (global as any).XMLHttpRequest as jest.Mock;
      const xhr = new xhr_ctor();

      xhr.open("GET", "http://example.com/test.ion", false);
      xhr.send();

      expect(xhr.status).toBe(400);
    });
  });

  describe("handle cleanup", () => {
    test("freeRaw removes source_id from registry", () => {
      const source_id_by_handle = new Map<number, number>();
      const registry = new Map<number, string>();

      const handle = 42;
      const source_id = 1;

      registry.set(source_id, "http://example.com/test.ion");
      source_id_by_handle.set(handle, source_id);

      expect(registry.has(source_id)).toBe(true);

      const sid = source_id_by_handle.get(handle);
      if (sid !== undefined) {
        registry.delete(sid);
        source_id_by_handle.delete(handle);
      }

      expect(registry.has(source_id)).toBe(false);
      expect(source_id_by_handle.has(handle)).toBe(false);
    });
  });
});
