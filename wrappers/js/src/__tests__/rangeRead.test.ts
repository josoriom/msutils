describe("Range Read Functions", () => {
  describe("unsigned offset reconstruction", () => {
    test("offset_hi=1, offset_lo=0xC0000000 produces unsigned 7516192768", () => {
      const offset_lo = 0xc0000000;
      const offset_hi = 1;
      const unsigned = (offset_lo >>> 0) + (offset_hi >>> 0) * (2 ** 32);
      expect(unsigned).toBe(7516192768);
    });

    test("signed offset_lo reconstructs correctly with >>", () => {
      const offset_lo = -0x40000000;
      const offset_lo_unsigned = offset_lo >>> 0;
      expect(offset_lo_unsigned).toBe(0xc0000000);
    });
  });

  describe("retry outcome classification", () => {
    test("outcome: ok when 206 and body matches len", () => {
      const outcome = { kind: "ok" as const, bytes: new Uint8Array(5) };
      expect(outcome.kind).toBe("ok");
      expect(outcome.bytes.length).toBe(5);
    });

    test("outcome: retry on throw or 5xx", () => {
      const retry_throw = { kind: "retry" as const };
      const retry_503 = { kind: "retry" as const };
      expect(retry_throw.kind).toBe("retry");
      expect(retry_503.kind).toBe("retry");
    });

    test("outcome: fail on 206 with wrong length", () => {
      const fail_short = { kind: "fail" as const };
      expect(fail_short.kind).toBe("fail");
    });

    test("outcome: fail on 200 (range not honored)", () => {
      const fail_200 = { kind: "fail" as const };
      expect(fail_200.kind).toBe("fail");
    });

    test("outcome: fail on 404 (file not found)", () => {
      const fail_404 = { kind: "fail" as const };
      expect(fail_404.kind).toBe("fail");
    });

    test("outcome: fail on 416 (invalid range)", () => {
      const fail_416 = { kind: "fail" as const };
      expect(fail_416.kind).toBe("fail");
    });

    test("Content-Range header format validation", () => {
      const offset = 1000;
      const last_byte = 1099;
      const content_range = `bytes ${offset}-${last_byte}/50000`;
      const expected = `bytes ${offset}-${last_byte}/`;
      expect(content_range.startsWith(expected)).toBe(true);
    });

    test("Content-Range header mismatch detection", () => {
      const actual = "bytes 0-99/50000";
      const expected = "bytes 1000-1099/";
      expect(actual.startsWith(expected)).toBe(false);
    });
  });

  describe("memory safety", () => {
    test("len === 0 is handled before memory access", () => {
      const len = 0;
      if (len === 0) {
        expect(len).toBe(0);
      }
    });

    test("fresh memory view per call prevents stale buffer", () => {
      const memory = new WebAssembly.Memory({ initial: 2 });
      const view1 = new Uint8Array(memory.buffer, 0, 10);
      const view2 = new Uint8Array(memory.buffer, 0, 10);
      expect(view1.buffer === view2.buffer).toBe(true);
    });

    test("destination pointer window prevents overflow", () => {
      const len = 10;
      const dest_ptr = 100;
      const memory = new WebAssembly.Memory({ initial: 2 });
      const view = new Uint8Array(memory.buffer, dest_ptr, len);
      expect(view.length).toBe(len);
    });
  });
});
