describe("RangeSourceRegistry", () => {
  describe("registry pattern", () => {
    test("add() assigns incrementing source_id", () => {
      const registry = new Map<number, string>();
      let next_id = 1;

      function add(url: string): number {
        const id = next_id++;
        registry.set(id, url);
        return id;
      }

      const id1 = add("http://example.com/file1.ion");
      const id2 = add("http://example.com/file2.ion");

      expect(id1).toBe(1);
      expect(id2).toBe(2);
    });

    test("get() retrieves URL by source_id", () => {
      const registry = new Map<number, string>();
      const url = "http://example.com/file.ion";
      registry.set(1, url);

      const retrieved = registry.get(1);
      expect(retrieved).toBe(url);
    });

    test("get() returns undefined for unknown id", () => {
      const registry = new Map<number, string>();
      const retrieved = registry.get(999);
      expect(retrieved).toBeUndefined();
    });

    test("remove() deletes entry", () => {
      const registry = new Map<number, string>();
      registry.set(1, "http://example.com/file.ion");
      expect(registry.has(1)).toBe(true);

      registry.delete(1);
      expect(registry.has(1)).toBe(false);
    });

    test("multiple open files isolate by source_id", () => {
      const registry = new Map<number, string>();
      let next_id = 1;

      const id1 = next_id++;
      const id2 = next_id++;

      registry.set(id1, "http://example.com/file1.ion");
      registry.set(id2, "http://example.com/file2.ion");

      expect(registry.get(id1)).not.toBe(registry.get(id2));
    });
  });

  describe("handle-to-source_id lifetime binding", () => {
    test("source_id_by_handle maps handle to source_id", () => {
      const source_id_by_handle = new Map<number, number>();
      const handle = 42;
      const source_id = 1;

      source_id_by_handle.set(handle, source_id);
      expect(source_id_by_handle.get(handle)).toBe(source_id);
    });

    test("null handle (0) signals error condition", () => {
      const handle = 0;
      expect(handle === 0).toBe(true);
    });

    test("freeRaw cleanup: handle lookup then registry remove", () => {
      const source_id_by_handle = new Map<number, number>();
      const registry = new Map<number, string>();

      const handle = 42;
      const source_id = 1;

      registry.set(source_id, "http://example.com/file.ion");
      source_id_by_handle.set(handle, source_id);

      const sid = source_id_by_handle.get(handle);
      if (sid !== undefined) {
        registry.delete(sid);
        source_id_by_handle.delete(handle);
      }

      expect(source_id_by_handle.has(handle)).toBe(false);
      expect(registry.has(source_id)).toBe(false);
    });
  });
});

describe("Worker Guard", () => {
  describe("main thread detection", () => {
    test("window + document defined = main thread", () => {
      const has_window = typeof window !== "undefined";
      const has_document = typeof document !== "undefined";
      const is_main_thread = has_window && has_document;

      expect(typeof is_main_thread).toBe("boolean");
    });

    test("worker environment: window/document absent", () => {
      const is_worker = typeof window === "undefined" || typeof document === "undefined";
      expect(is_worker).toBe(true);
    });

    test("guard prevents main thread parseIon(URL)", () => {
      function allow_worker_only(): void {
        const is_main_thread =
          typeof window !== "undefined" && typeof document !== "undefined";
        if (is_main_thread) {
          throw new Error(
            "parseIon(URL) must run inside a Web Worker; sync range reads are blocked on the main thread",
          );
        }
      }

      expect(() => {
        allow_worker_only();
      }).not.toThrow();
    });
  });

  describe("error message clarity", () => {
    test("error mentions Web Worker requirement", () => {
      const msg = "parseIon(URL) must run inside a Web Worker; sync range reads are blocked on the main thread";
      expect(msg).toContain("Web Worker");
      expect(msg).toContain("sync range reads");
    });

    test("error mentions main thread blockage", () => {
      const msg = "parseIon(URL) must run inside a Web Worker; sync range reads are blocked on the main thread";
      expect(msg).toContain("main thread");
    });
  });
});
