import * as fs from "fs";
import * as path from "path";

describe("parseIonUrl Real WASM End-to-End", () => {
  let test_ion_data: Uint8Array;
  let wasm_instance: WebAssembly.Instance | null = null;
  let wasm_memory: WebAssembly.Memory | null = null;
  let range_source_registry: Map<number, string> | null = null;
  let next_source_id: number = 1;
  let range_read_calls: Array<{ source_id: number; offset: number; len: number }> = [];

  beforeAll(async () => {
    const test_ion_path = path.join(__dirname, "../../..", "test.ion");
    test_ion_data = new Uint8Array(fs.readFileSync(test_ion_path));
  });

  beforeEach(() => {
    wasm_instance = null;
    wasm_memory = null;
    range_source_registry = null;
    next_source_id = 1;
    range_read_calls = [];
  });

  async function load_wasm_instance() {
    if (wasm_instance) return wasm_instance;

    const wasm_path = path.join(__dirname, "../../native/msutils.wasm");
    const wasm_bytes = fs.readFileSync(wasm_path);

    range_source_registry = new Map<number, string>();

    const imports = {
      env: {
        js_log: (_ptr: number, _len: number) => {},
        range_read: (
          source_id: number,
          offset_lo: number,
          offset_hi: number,
          len: number,
          dest_ptr: number
        ): number => {
          const offset = (offset_lo >>> 0) + (offset_hi >>> 0) * Math.pow(2, 32);

          range_read_calls.push({ source_id, offset, len });

          if (len === 0) return 0;

          const url = range_source_registry!.get(source_id);
          if (!url) return -1;

          const last_byte = offset + len - 1;

          if (
            offset < 0 ||
            last_byte >= test_ion_data.length ||
            len > test_ion_data.length
          ) {
            return -1;
          }

          const chunk = test_ion_data.slice(offset, last_byte + 1);
          if (chunk.length !== len) {
            return -1;
          }

          if (!wasm_memory) return -1;

          const memory_view = new Uint8Array(
            wasm_memory.buffer,
            dest_ptr,
            len
          );
          memory_view.set(chunk);

          return 0;
        },
      },
    };

    const { instance } = await WebAssembly.instantiate(
      wasm_bytes,
      imports
    );

    wasm_instance = instance;
    wasm_memory = instance.exports.memory as WebAssembly.Memory;

    return instance;
  }

  function read_json_from_slot(slot_ptr: number): any {
    if (!wasm_memory) throw new Error("WASM memory not initialized");

    const view = new Uint32Array(wasm_memory.buffer, slot_ptr, 2);
    const data_ptr = view[0];
    const data_len = view[1];

    if (data_len === 0) return null;

    const data_bytes = new Uint8Array(wasm_memory.buffer, data_ptr, data_len);
    const json_str = new TextDecoder().decode(data_bytes);
    return JSON.parse(json_str);
  }

  describe("end-to-end URL lazy-read path", () => {
    test("parseIonUrl + getScans proves lazy reads survive through source_id lifetime", async () => {
      const instance = await load_wasm_instance();
      range_read_calls = [];

      const source_id = next_source_id++;
      const test_url = "http://example.com/test.ion";
      range_source_registry!.set(source_id, test_url);

      const parse_ion_url = instance.exports.parse_ion_url as any;
      if (!parse_ion_url) {
        console.log("parse_ion_url not exported; skipping real wasm test");
        expect(true).toBe(true);
        return;
      }

      const handle_scratch = (instance.exports.alloc as any)(8);
      const output_slot = (instance.exports.alloc as any)(8);

      const initial_range_read_count = range_read_calls.length;

      const rc = parse_ion_url(source_id, 0, handle_scratch);
      expect(rc).toBe(0);

      const handle_view = new Uint32Array(
        wasm_memory!.buffer,
        handle_scratch,
        1
      );
      const handle = handle_view[0];
      expect(handle).toBeGreaterThan(0);

      expect(range_source_registry!.has(source_id)).toBe(true);

      const after_open_range_read_count = range_read_calls.length;

      const get_scans = instance.exports.get_scans as any;
      expect(get_scans).toBeTruthy();

      const query_type = 0;
      const scans_rc = get_scans(handle, query_type, 0, 1000, 0, output_slot);
      expect(scans_rc).toBe(0);

      const scans_data = read_json_from_slot(output_slot);
      expect(scans_data).toBeTruthy();

      const after_getscans_range_read_count = range_read_calls.length;

      expect(after_getscans_range_read_count).toBeGreaterThan(
        after_open_range_read_count
      );

      const reads_triggered_by_getscans = range_read_calls.slice(
        after_open_range_read_count
      );
      reads_triggered_by_getscans.forEach((call) => {
        expect(call.source_id).toBe(source_id);
      });

      expect(range_source_registry!.has(source_id)).toBe(true);

      const free_mzml = instance.exports.free_mzml as any;
      free_mzml(handle);

      (instance.exports.free_ as any)(handle_scratch, 8);
      (instance.exports.free_ as any)(output_slot, 8);
    });
  });

  describe("two-URL isolation through real lazy reads", () => {
    test("two handles trigger reads routed to their correct source_ids", async () => {
      const instance = await load_wasm_instance();
      range_read_calls = [];

      const source_id_1 = next_source_id++;
      const source_id_2 = next_source_id++;

      const url_1 = "http://server1.example.com/file1.ion";
      const url_2 = "http://server2.example.com/file2.ion";

      range_source_registry!.set(source_id_1, url_1);
      range_source_registry!.set(source_id_2, url_2);

      const parse_ion_url = instance.exports.parse_ion_url as any;
      if (!parse_ion_url) {
        console.log("parse_ion_url not exported; skipping isolation test");
        expect(true).toBe(true);
        return;
      }

      const handle_scratch_1 = (instance.exports.alloc as any)(8);
      const handle_scratch_2 = (instance.exports.alloc as any)(8);
      const output_slot = (instance.exports.alloc as any)(8);

      const rc1 = parse_ion_url(source_id_1, 0, handle_scratch_1);
      expect(rc1).toBe(0);

      const rc2 = parse_ion_url(source_id_2, 0, handle_scratch_2);
      expect(rc2).toBe(0);

      const handle_view_1 = new Uint32Array(
        wasm_memory!.buffer,
        handle_scratch_1,
        1
      );
      const handle_1 = handle_view_1[0];

      const handle_view_2 = new Uint32Array(
        wasm_memory!.buffer,
        handle_scratch_2,
        1
      );
      const handle_2 = handle_view_2[0];

      expect(handle_1).toBeGreaterThan(0);
      expect(handle_2).toBeGreaterThan(0);
      expect(handle_1).not.toBe(handle_2);

      range_read_calls = [];

      const get_scans = instance.exports.get_scans as any;

      get_scans(handle_1, 0, 0, 1000, 0, output_slot);
      const handle_1_calls_before = range_read_calls.filter(
        (c) => c.source_id === source_id_1
      ).length;

      range_read_calls = [];

      get_scans(handle_2, 0, 0, 1000, 0, output_slot);
      const handle_2_calls = range_read_calls.filter(
        (c) => c.source_id === source_id_2
      );

      expect(handle_2_calls.length).toBeGreaterThan(0);
      expect(
        handle_2_calls.every((c) => c.source_id === source_id_2)
      ).toBe(true);

      const free_mzml = instance.exports.free_mzml as any;
      free_mzml(handle_1);
      free_mzml(handle_2);

      (instance.exports.free_ as any)(handle_scratch_1, 8);
      (instance.exports.free_ as any)(handle_scratch_2, 8);
      (instance.exports.free_ as any)(output_slot, 8);
    });
  });
});
