/// <reference types="node" />

declare const __INLINE__: boolean;
declare const __WASM_DATA_URL__: string | undefined;
declare var require: any;

import type { Backend, FileHandle } from "./backend";

const IS_INLINE_BUILD = typeof __INLINE__ !== "undefined" && __INLINE__;
const OUTPUT_BUFFER_PAIR_BYTES = 8;
const PEAK_OPTIONS_STRUCT_BYTES = 64;

function resolveTextDecoder(): TextDecoder {
  return typeof TextDecoder !== "undefined"
    ? new TextDecoder("utf-8")
    : new (require("node:util").TextDecoder)("utf-8");
}

function writeF64(view: DataView, offset: number, value: unknown): void {
  view.setFloat64(
    offset,
    Number.isFinite(value as number) ? Number(value) : NaN,
    true,
  );
}

function writeI32(view: DataView, offset: number, value: unknown): void {
  view.setInt32(
    offset,
    typeof value === "boolean"
      ? value
        ? 1
        : 0
      : Number.isFinite(value as number)
        ? (value as number) | 0
        : 0,
    true,
  );
}

function toUint8View(typed: {
  buffer: ArrayBufferLike;
  byteOffset: number;
  byteLength: number;
}): Uint8Array {
  const buffer =
    typed.buffer instanceof ArrayBuffer
      ? typed.buffer
      : new Uint8Array(typed.buffer).buffer;
  return new Uint8Array(buffer, typed.byteOffset, typed.byteLength);
}

class WasmHeap {
  private memory: WebAssembly.Memory;
  private bytes: Uint8Array;
  private view: DataView;
  private readonly alloc: (size: number) => number;
  private readonly free: (ptr: number, size: number) => void;
  readonly decoder: TextDecoder;

  constructor(
    memory: WebAssembly.Memory,
    alloc: (size: number) => number,
    free: (ptr: number, size: number) => void,
  ) {
    this.memory = memory;
    this.bytes = new Uint8Array(memory.buffer);
    this.view = new DataView(memory.buffer);
    this.alloc = alloc;
    this.free = free;
    this.decoder = resolveTextDecoder();
  }

  private refresh(): void {
    if (this.bytes.buffer !== this.memory.buffer) {
      this.bytes = new Uint8Array(this.memory.buffer);
      this.view = new DataView(this.memory.buffer);
    }
  }

  allocAndWrite(data: Uint8Array): [ptr: number, len: number] {
    const ptr = this.alloc(data.length);
    this.refresh();
    this.bytes.set(data, ptr);
    return [ptr, data.length];
  }

  freeMany(allocations: [ptr: number, len: number][]): void {
    for (const [ptr, len] of allocations) this.free(ptr, len);
  }

  readU32(ptr: number): number {
    this.refresh();
    return this.view.getUint32(ptr, true);
  }

  readOutputSlot(slotPtr: number): { ptr: number; len: number } {
    this.refresh();
    return {
      ptr: this.view.getUint32(slotPtr, true),
      len: this.view.getUint32(slotPtr + 4, true),
    };
  }

  copyOut(ptr: number, len: number): Uint8Array {
    this.refresh();
    const copy = new Uint8Array(len);
    copy.set(this.bytes.subarray(ptr, ptr + len));
    return copy;
  }

  copyOutAndFree(ptr: number, len: number): Uint8Array {
    const data = this.copyOut(ptr, len);
    this.free(ptr, len);
    return data;
  }

  readJsonFromSlot<T>(slotPtr: number): T {
    const { ptr, len } = this.readOutputSlot(slotPtr);
    const bytes = this.copyOutAndFree(ptr, len);
    return JSON.parse(this.decoder.decode(bytes)) as T;
  }

  writeI32At(ptr: number, offset: number, value: unknown): void {
    this.refresh();
    writeI32(this.view, ptr + offset, value);
  }

  writeF64At(ptr: number, offset: number, value: unknown): void {
    this.refresh();
    writeF64(this.view, ptr + offset, value);
  }

  writeBytesAt(ptr: number, data: Uint8Array): void {
    this.refresh();
    this.bytes.set(data, ptr);
  }

  allocPermanent(size: number): number {
    return this.alloc(size);
  }
}

class WasmExports {
  readonly memory: WebAssembly.Memory;
  readonly alloc: (size: number) => number;
  readonly free: (ptr: number, size: number) => void;
  readonly parseMzml: (
    dataPtr: number,
    dataLen: number,
    outPtr: number,
  ) => number;
  readonly parseBin: (
    dataPtr: number,
    dataLen: number,
    maxCacheSize: number,
    outPtr: number,
  ) => number;
  readonly freeMzml: (handle: number) => void;
  readonly binToJson: (handle: number, outBuf: number) => number;
  readonly binToMzml: (handle: number, outBuf: number) => number;
  readonly mzmlToBin: (
    handle: number,
    outBuf: number,
    level: number,
    compress: number,
  ) => number;
  readonly calculateEic: (
    handle: number,
    target: number,
    from: number,
    to: number,
    ppm: number,
    mz: number,
    outX: number,
    outY: number,
  ) => number;
  readonly findPeaks: (
    xPtr: number,
    yPtr: number,
    len: number,
    optsPtr: number,
    outBuf: number,
  ) => number;
  readonly getPeak: (
    xPtr: number,
    yPtr: number,
    len: number,
    rt: number,
    range: number,
    optsPtr: number,
    outBuf: number,
  ) => number;
  readonly getPeaksFromChrom: (
    handle: number,
    idxPtr: number,
    rtPtr: number,
    winPtr: number,
    count: number,
    optsPtr: number,
    cores: number,
    outBuf: number,
  ) => number;
  readonly getPeaksFromEic: (
    handle: number,
    rtPtr: number,
    mzPtr: number,
    rangePtr: number,
    offPtr: number,
    lenPtr: number,
    idBytesPtr: number,
    idBytesLen: number,
    count: number,
    from: number,
    to: number,
    optsPtr: number,
    cores: number,
    outBuf: number,
  ) => number;
  readonly findNoiseLevel: (yPtr: number, len: number) => number;
  readonly calculateBaseline: (
    yPtr: number,
    len: number,
    window: number,
    windowFactor: number,
    outBuf: number,
  ) => number;
  readonly findFeature: (
    handle: number,
    rtPtr: number,
    mzPtr: number,
    rangePtr: number,
    offPtr: number,
    lenPtr: number,
    idBytesPtr: number,
    idBytesLen: number,
    count: number,
    cores: number,
    scanPpm: number,
    scanMz: number,
    eicPpm: number,
    eicMz: number,
    optsPtr: number,
    outBuf: number,
  ) => number;
  readonly getScans: (
    handle: number,
    queryType: number,
    a: number,
    b: number,
    level: number,
    outBuf: number,
  ) => number;

  constructor(instance: WebAssembly.Instance) {
    const ex = instance.exports;
    this.memory = ex.memory as WebAssembly.Memory;
    this.alloc = this.resolve(ex, ["alloc"]);
    this.free = this.resolve(ex, ["free_", "free"]);
    this.parseMzml = this.resolve(ex, ["parse_mzml"]);
    this.parseBin = this.resolve(ex, ["parse_bin"]);
    this.freeMzml = this.resolve(ex, ["free_mzml"]);
    this.binToJson = this.resolve(ex, ["bin_to_json"]);
    this.binToMzml = this.resolve(ex, ["bin_to_mzml"]);
    this.mzmlToBin = this.resolve(ex, ["mzml_to_bin"]);
    this.calculateEic = this.resolve(ex, ["calculate_eic"]);
    this.findPeaks = this.resolve(ex, ["find_peaks"]);
    this.getPeak = this.resolve(ex, ["get_peak"]);
    this.getPeaksFromChrom = this.resolve(ex, ["get_peaks_from_chrom"]);
    this.getPeaksFromEic = this.resolve(ex, ["get_peaks_from_eic"]);
    this.findNoiseLevel = this.resolve(ex, ["find_noise_level"]);
    this.calculateBaseline = this.resolve(ex, ["calculate_baseline"]);
    this.findFeature = this.resolve(ex, ["find_feature"]);
    this.getScans = this.resolve(ex, ["get_scans"]);
  }

  private resolve<T extends Function>(
    exports: WebAssembly.Exports,
    names: string[],
  ): T {
    for (const name of names)
      if (typeof exports[name] === "function")
        return exports[name] as unknown as T;
    throw new Error(
      `WasmExports: none of [${names.join(", ")}] found in WASM exports`,
    );
  }
}

class WasmApi {
  private readonly heap: WasmHeap;
  readonly fn: WasmExports;

  private readonly jsonOutputSlot: number;
  private readonly jsonScratchSlot: number;
  private readonly blobScratchSlot: number;
  private readonly baselineScratchSlot: number;
  private readonly peakOptionsSlot: number;
  private readonly handleScratchSlot: number;

  constructor(instance: WebAssembly.Instance) {
    this.fn = new WasmExports(instance);
    this.heap = new WasmHeap(this.fn.memory, this.fn.alloc, this.fn.free);

    this.jsonOutputSlot = this.heap.allocPermanent(OUTPUT_BUFFER_PAIR_BYTES);
    this.jsonScratchSlot = this.heap.allocPermanent(OUTPUT_BUFFER_PAIR_BYTES);
    this.blobScratchSlot = this.heap.allocPermanent(OUTPUT_BUFFER_PAIR_BYTES);
    this.baselineScratchSlot = this.heap.allocPermanent(
      OUTPUT_BUFFER_PAIR_BYTES,
    );
    this.peakOptionsSlot = this.heap.allocPermanent(PEAK_OPTIONS_STRUCT_BYTES);
    this.handleScratchSlot = this.heap.allocPermanent(4);
  }

  parseMzMLRaw(data: Uint8Array): number {
    const [ptr, len] = this.heap.allocAndWrite(data);
    try {
      const rc = this.fn.parseMzml(ptr, len, this.handleScratchSlot);
      if (rc !== 0) throw new Error("parse_mzml failed with code " + rc);
      return this.heap.readU32(this.handleScratchSlot);
    } finally {
      this.fn.free(ptr, len);
    }
  }

  parseBinRaw(data: Uint8Array, maxCacheSize = 0): number {
    const [ptr, len] = this.heap.allocAndWrite(data);
    try {
      const rc = this.fn.parseBin(
        ptr,
        len,
        maxCacheSize,
        this.handleScratchSlot,
      );
      if (rc !== 0) throw new Error("parse_bin failed with code " + rc);
      return this.heap.readU32(this.handleScratchSlot);
    } finally {
      this.fn.free(ptr, len);
    }
  }

  freeRaw(handle: number): void {
    this.fn.freeMzml(handle);
  }

  fileToJson(handle: number): any {
    const rc = this.fn.binToJson(handle, this.jsonOutputSlot);
    if (rc !== 0) throw new Error("bin_to_json failed with code " + rc);
    return this.heap.readJsonFromSlot<any>(this.jsonOutputSlot);
  }

  fileToMzml(handle: number): string {
    const rc = this.fn.binToMzml(handle, this.jsonOutputSlot);
    if (rc !== 0) throw new Error("bin_to_mzml failed with code " + rc);
    const { ptr, len } = this.heap.readOutputSlot(this.jsonOutputSlot);
    return this.heap.decoder.decode(this.heap.copyOutAndFree(ptr, len));
  }

  fileToBin(
    handle: number,
    compressionLevel: number,
    f32Compress: number,
  ): Uint8Array {
    const rc = this.fn.mzmlToBin(
      handle,
      this.blobScratchSlot,
      compressionLevel,
      f32Compress,
    );
    if (rc !== 0) throw new Error("mzml_to_bin failed with code " + rc);
    const { ptr, len } = this.heap.readOutputSlot(this.blobScratchSlot);
    return this.heap.copyOutAndFree(ptr, len);
  }

  calculateEic(
    handle: number,
    targetMz: number,
    from: number,
    to: number,
    ppmTolerance: number,
    mzTolerance: number,
  ): { x: Float64Array; y: Float64Array } {
    const rc = this.fn.calculateEic(
      handle,
      targetMz,
      from,
      to,
      ppmTolerance,
      mzTolerance,
      this.jsonScratchSlot,
      this.blobScratchSlot,
    );
    if (rc !== 0) throw new Error("calculate_eic failed with code " + rc);
    const xSlot = this.heap.readOutputSlot(this.jsonScratchSlot);
    const ySlot = this.heap.readOutputSlot(this.blobScratchSlot);
    const x = new Float64Array(
      this.heap.copyOutAndFree(xSlot.ptr, xSlot.len).buffer,
    );
    const y = new Float64Array(
      this.heap.copyOutAndFree(ySlot.ptr, ySlot.len).buffer,
    );
    return { x, y };
  }

  getScans(handle: number, queryType: number, a: number, b: number, level: number): any {
    const rc = this.fn.getScans(handle, queryType, a, b, level, this.jsonOutputSlot);
    if (rc !== 0) throw new Error("get_scans failed with code " + rc);
    return this.heap.readJsonFromSlot<any>(this.jsonOutputSlot);
  }

  findPeaks(
    x: Float64Array,
    y: Float64Array,
    packedOpts: Uint8Array | undefined,
  ): any {
    const [xPtr, xLen] = this.heap.allocAndWrite(toUint8View(x));
    const [yPtr, yLen] = this.heap.allocAndWrite(toUint8View(y));
    const rc = this.fn.findPeaks(
      xPtr,
      yPtr,
      x.length,
      this.peakOptsPtr(packedOpts),
      this.jsonOutputSlot,
    );
    this.heap.freeMany([
      [xPtr, xLen],
      [yPtr, yLen],
    ]);
    if (rc !== 0) throw new Error("find_peaks failed with code " + rc);
    return this.heap.readJsonFromSlot<any>(this.jsonOutputSlot);
  }

  getPeak(
    x: Float64Array,
    y: Float64Array,
    rt: number,
    range: number,
    packedOpts: Uint8Array | undefined,
  ): any {
    const [xPtr, xLen] = this.heap.allocAndWrite(toUint8View(x));
    const [yPtr, yLen] = this.heap.allocAndWrite(toUint8View(y));
    const rc = this.fn.getPeak(
      xPtr,
      yPtr,
      x.length,
      rt,
      range,
      this.peakOptsPtr(packedOpts),
      this.jsonOutputSlot,
    );
    this.heap.freeMany([
      [xPtr, xLen],
      [yPtr, yLen],
    ]);
    if (rc !== 0) throw new Error("get_peak failed with code " + rc);
    return this.heap.readJsonFromSlot<any>(this.jsonOutputSlot);
  }

  findNoiseLevel(y: Float32Array): number {
    if (y.length === 0) return Infinity;
    const [yPtr, yLen] = this.heap.allocAndWrite(toUint8View(y));
    const noise = this.fn.findNoiseLevel(yPtr, y.length);
    this.fn.free(yPtr, yLen);
    return noise;
  }

  calculateBaseline(
    y: Float64Array,
    baselineWindow: number,
    baselineWindowFactor: number,
  ): Float64Array {
    const [yPtr, yLen] = this.heap.allocAndWrite(toUint8View(y));
    const rc = this.fn.calculateBaseline(
      yPtr,
      y.length,
      baselineWindow,
      baselineWindowFactor,
      this.baselineScratchSlot,
    );
    this.fn.free(yPtr, yLen);
    if (rc !== 0) throw new Error("calculate_baseline failed with code " + rc);
    const { ptr, len } = this.heap.readOutputSlot(this.baselineScratchSlot);
    return new Float64Array(this.heap.copyOutAndFree(ptr, len).buffer.slice(0));
  }

  getPeaksFromEic(
    handle: number,
    rts: Float64Array,
    mzs: Float64Array,
    ranges: Float64Array,
    offsets: Uint32Array | null,
    lengths: Uint32Array | null,
    idBytes: Uint8Array | null,
    idBytesLen: number,
    count: number,
    from: number,
    to: number,
    packedOpts: Uint8Array | undefined,
    cores: number,
  ): any {
    const emptyU32 = new Uint32Array(count);
    const emptyBytes = new Uint8Array(0);
    const [rtPtr, rtLen] = this.heap.allocAndWrite(toUint8View(rts));
    const [mzPtr, mzLen] = this.heap.allocAndWrite(toUint8View(mzs));
    const [rangePtr, rangeLen] = this.heap.allocAndWrite(toUint8View(ranges));
    const [offPtr, offLen] = this.heap.allocAndWrite(
      toUint8View(offsets ?? emptyU32),
    );
    const [lenPtr, lenLen] = this.heap.allocAndWrite(
      toUint8View(lengths ?? emptyU32),
    );
    const [idPtr, idLen] = this.heap.allocAndWrite(idBytes ?? emptyBytes);
    const rc = this.fn.getPeaksFromEic(
      handle,
      rtPtr,
      mzPtr,
      rangePtr,
      offPtr,
      lenPtr,
      idPtr,
      idBytesLen,
      count,
      from,
      to,
      this.peakOptsPtr(packedOpts),
      cores,
      this.jsonOutputSlot,
    );
    this.heap.freeMany([
      [rtPtr, rtLen],
      [mzPtr, mzLen],
      [rangePtr, rangeLen],
      [offPtr, offLen],
      [lenPtr, lenLen],
      [idPtr, idLen],
    ]);
    if (rc !== 0) throw new Error("get_peaks_from_eic failed with code " + rc);
    return this.heap.readJsonFromSlot<any>(this.jsonOutputSlot);
  }

  getPeaksFromChrom(
    handle: number,
    indices: Uint32Array,
    rts: Float64Array,
    windows: Float64Array,
    count: number,
    packedOpts: Uint8Array | undefined,
    cores: number,
  ): any {
    const [idxPtr, idxLen] = this.heap.allocAndWrite(toUint8View(indices));
    const [rtPtr, rtLen] = this.heap.allocAndWrite(toUint8View(rts));
    const [winPtr, winLen] = this.heap.allocAndWrite(toUint8View(windows));
    const rc = this.fn.getPeaksFromChrom(
      handle,
      idxPtr,
      rtPtr,
      winPtr,
      count,
      this.peakOptsPtr(packedOpts),
      cores,
      this.jsonOutputSlot,
    );
    this.heap.freeMany([
      [idxPtr, idxLen],
      [rtPtr, rtLen],
      [winPtr, winLen],
    ]);
    if (rc !== 0)
      throw new Error("get_peaks_from_chrom failed with code " + rc);
    return this.heap.readJsonFromSlot<any>(this.jsonOutputSlot);
  }

  findFeature(
    handle: number,
    rts: Float64Array,
    mzs: Float64Array,
    ranges: Float64Array,
    offsets: Uint32Array | null,
    lengths: Uint32Array | null,
    idBytes: Uint8Array | null,
    idBytesLen: number,
    count: number,
    cores: number,
    scanPpm: number,
    scanMz: number,
    eicPpm: number,
    eicMz: number,
    packedOpts: Uint8Array | undefined,
  ): any {
    const emptyU32 = new Uint32Array(count);
    const emptyBytes = new Uint8Array(0);
    const [rtPtr, rtLen] = this.heap.allocAndWrite(toUint8View(rts));
    const [mzPtr, mzLen] = this.heap.allocAndWrite(toUint8View(mzs));
    const [rangePtr, rangeLen] = this.heap.allocAndWrite(toUint8View(ranges));
    const [offPtr, offLen] = this.heap.allocAndWrite(
      toUint8View(offsets ?? emptyU32),
    );
    const [lenPtr, lenLen] = this.heap.allocAndWrite(
      toUint8View(lengths ?? emptyU32),
    );
    const [idPtr, idLen] = this.heap.allocAndWrite(idBytes ?? emptyBytes);
    const rc = this.fn.findFeature(
      handle,
      rtPtr,
      mzPtr,
      rangePtr,
      offPtr,
      lenPtr,
      idPtr,
      idBytesLen,
      count,
      cores,
      scanPpm,
      scanMz,
      eicPpm,
      eicMz,
      this.peakOptsPtr(packedOpts),
      this.jsonOutputSlot,
    );
    this.heap.freeMany([
      [rtPtr, rtLen],
      [mzPtr, mzLen],
      [rangePtr, rangeLen],
      [offPtr, offLen],
      [lenPtr, lenLen],
      [idPtr, idLen],
    ]);
    if (rc !== 0) throw new Error("find_feature failed with code " + rc);
    return this.heap.readJsonFromSlot<any>(this.jsonOutputSlot);
  }

  private peakOptsPtr(packed: Uint8Array | undefined): number {
    if (!packed) return 0;
    if (packed.length !== PEAK_OPTIONS_STRUCT_BYTES)
      throw new RangeError(
        `peakOptsPtr: expected ${PEAK_OPTIONS_STRUCT_BYTES} bytes, got ${packed.length}`,
      );
    this.heap.writeBytesAt(this.peakOptionsSlot, packed);
    return this.peakOptionsSlot;
  }
}

async function loadWasm(): Promise<WasmApi> {
  const imports = { env: { js_log: (_ptr: number, _len: number) => {} } };

  if (IS_INLINE_BUILD) {
    if (!__WASM_DATA_URL__)
      throw new Error("Inline build is missing __WASM_DATA_URL__");
    const bytes = await fetch(__WASM_DATA_URL__).then((r) => r.arrayBuffer());
    const { instance } = await WebAssembly.instantiate(bytes, imports);
    return new WasmApi(instance);
  }

  if (typeof process !== "undefined" && (process as any).versions?.node) {
    const nfs = require("node:fs");
    const npath = require("node:path");
    const nurl = require("node:url");
    const dir =
      typeof __dirname !== "undefined"
        ? __dirname
        : npath.dirname(nurl.fileURLToPath((0, eval)("import.meta").url));
    const bytes = nfs.readFileSync(npath.join(dir, "msutils.wasm"));
    const { instance } = await WebAssembly.instantiate(bytes, imports);
    return new WasmApi(instance);
  }

  throw new Error(
    "Browser ESM build requires the UMD inline bundle (msutils.js)",
  );
}

export class WasmBackend implements Backend {
  private api: WasmApi | null = null;
  ready = false;
  readonly handlesAreGcFinalized = false;

  async init(): Promise<void> {
    if (this.api) return;
    this.api = await loadWasm();
    this.ready = true;
  }

  private getApi(): WasmApi {
    if (!this.api)
      throw new Error("WasmBackend not initialized — call init() first");
    return this.api;
  }

  parseMzML(data: Uint8Array): FileHandle {
    return this.getApi().parseMzMLRaw(data);
  }

  parseBin(data: Uint8Array, maxCacheSize = 0): FileHandle {
    return this.getApi().parseBinRaw(data, maxCacheSize);
  }

  freeFile(handle: FileHandle): void {
    this.getApi().freeRaw(handle as number);
  }

  fileToJson(handle: FileHandle): any {
    return this.getApi().fileToJson(handle as number);
  }

  fileToMzml(handle: FileHandle): string {
    return this.getApi().fileToMzml(handle as number);
  }

  fileToBin(
    handle: FileHandle,
    level: number,
    f32Compress: boolean,
  ): Uint8Array {
    return this.getApi().fileToBin(
      handle as number,
      level,
      f32Compress ? 1 : 0,
    );
  }

  calculateEic(
    handle: FileHandle,
    targetMz: number,
    from: number,
    to: number,
    ppmTol: number,
    mzTol: number,
  ): { x: Float64Array; y: Float64Array } {
    return this.getApi().calculateEic(
      handle as number,
      targetMz,
      from,
      to,
      ppmTol,
      mzTol,
    );
  }

  getScans(
    handle: FileHandle,
    queryType: number,
    a: number,
    b: number,
    level: number,
  ): any {
    return this.getApi().getScans(handle as number, queryType, a, b, level);
  }

  findPeaks(
    x: Float64Array,
    y: Float64Array,
    packedOpts: Uint8Array | undefined,
  ): any {
    return this.getApi().findPeaks(x, y, packedOpts);
  }

  getPeak(
    x: Float64Array,
    y: Float64Array,
    rt: number,
    range: number,
    packedOpts: Uint8Array | undefined,
  ): any {
    return this.getApi().getPeak(x, y, rt, range, packedOpts);
  }

  findNoiseLevel(y: Float64Array | Float32Array): number {
    const y32 = y instanceof Float32Array ? y : new Float32Array(y);
    return this.getApi().findNoiseLevel(y32);
  }

  calculateBaseline(
    y: Float64Array,
    baselineWindow: number,
    baselineWindowFactor: number,
  ): Float64Array {
    return this.getApi().calculateBaseline(
      y,
      baselineWindow,
      baselineWindowFactor,
    );
  }

  getPeaksFromEic(
    handle: FileHandle,
    rts: Float64Array,
    mzs: Float64Array,
    ranges: Float64Array,
    offsets: Uint32Array | null,
    lengths: Uint32Array | null,
    idBytes: Uint8Array | null,
    idBytesLen: number,
    count: number,
    from: number,
    to: number,
    packedOpts: Uint8Array | undefined,
    cores: number,
  ): any {
    return this.getApi().getPeaksFromEic(
      handle as number,
      rts,
      mzs,
      ranges,
      offsets,
      lengths,
      idBytes,
      idBytesLen,
      count,
      from,
      to,
      packedOpts,
      cores,
    );
  }

  getPeaksFromChrom(
    handle: FileHandle,
    indices: Uint32Array,
    rts: Float64Array,
    windows: Float64Array,
    count: number,
    packedOpts: Uint8Array | undefined,
    cores: number,
  ): any {
    return this.getApi().getPeaksFromChrom(
      handle as number,
      indices,
      rts,
      windows,
      count,
      packedOpts,
      cores,
    );
  }

  findFeatures(
    _handle: FileHandle,
    _from: number,
    _to: number,
    _eicPpm: number,
    _eicMz: number,
    _gridStart: number,
    _gridEnd: number,
    _gridStep: number,
    _packedOpts: Uint8Array | undefined,
    _cores: number,
  ): any {
    throw new Error("findFeatures is not supported in the WASM build");
  }

  findFeature(
    handle: FileHandle,
    rts: Float64Array,
    mzs: Float64Array,
    ranges: Float64Array,
    offsets: Uint32Array | null,
    lengths: Uint32Array | null,
    idBytes: Uint8Array | null,
    idBytesLen: number,
    count: number,
    cores: number,
    scanPpm: number,
    scanMz: number,
    eicPpm: number,
    eicMz: number,
    packedOpts: Uint8Array | undefined,
  ): any {
    return this.getApi().findFeature(
      handle as number,
      rts,
      mzs,
      ranges,
      offsets,
      lengths,
      idBytes,
      idBytesLen,
      count,
      cores,
      scanPpm,
      scanMz,
      eicPpm,
      eicMz,
      packedOpts,
    );
  }
}
