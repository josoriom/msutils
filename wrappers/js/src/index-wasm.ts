/// <reference types="node" />

declare const __INLINE__: boolean;
declare const __WASM_DATA_URL__: string | undefined;
declare var require: any;

const IS_INLINE_BUILD = typeof __INLINE__ !== "undefined" && __INLINE__;
const OUTPUT_BUFFER_PAIR_BYTES = 8;
const PEAK_OPTIONS_STRUCT_BYTES = 64;

export type PeakOptions = Partial<{
  integralThreshold: number;
  intensityThreshold: number;
  widthThreshold: number;
  noise: number;
  autoNoise: boolean | number;
  autoBaseline: boolean | number;
  baselineWindow: number;
  baselineWindowFactor: number;
  allowOverlap: boolean | number;
  windowSize: number;
  snRatio: number;
}>;

export type BaselineOptions = Partial<{
  baselineWindow: number;
  baselineWindowFactor: number;
}>;

export type Peak = {
  from: number;
  to: number;
  rt: number;
  integral: number;
  intensity: number;
  ratio: number;
  np: number;
};

export type Target = {
  id?: string;
  rt: number;
  mz: number;
  range: number;
};

export type ChromItem = {
  id?: string;
  idx?: number;
  index?: number;
  rt: number;
  window?: number;
  range?: number;
};

export type ChromPeak = {
  index?: number;
  id?: string;
  ort: number;
  rt: number;
  from: number;
  to: number;
  intensity: number;
  integral: number;
};

export type FindFeaturesOptions = {
  eic?: { ppmTolerance?: number; mzTolerance?: number };
  grid?: { start?: number; end?: number; stepSize?: number };
  findPeak?: PeakOptions;
  cores?: number;
};

export type Feature = {
  mz: number;
  rt: number;
  from: number;
  to: number;
  intensity: number;
  integral: number;
  ratio: number;
  np: number;
};

export type FoundFeature = {
  id: string;
  mz: number;
  rt: number;
  from: number;
  to: number;
  intensity: number;
  integral: number;
  ratio: number;
  np: number;
  noise: number;
};

export interface FromTo {
  from: number;
  to: number;
}

export type ScanMetadataEntry = {
  section: string;
  accession?: string;
  name: string;
  value: string;
  unitAccession?: string;
  unitName?: string;
};

export type CentroidScan = {
  rt: number;
  mz: number[];
  intensity: number[];
  metadata: ScanMetadataEntry[];
};

type BinaryInput = Uint8Array | ArrayBuffer;

const snakeToCamel = (s: string): string =>
  s.replace(/_([a-z])/g, (_, char) => char.toUpperCase());

function camelizeKeys(value: any): any {
  if (Array.isArray(value)) return value.map(camelizeKeys);
  if (value && typeof value === "object" && !ArrayBuffer.isView(value)) {
    const result: any = {};
    for (const [key, val] of Object.entries(value))
      result[snakeToCamel(key)] = camelizeKeys(val);
    return result;
  }
  return value;
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

function toUint8(input: BinaryInput): Uint8Array {
  return input instanceof Uint8Array ? input : new Uint8Array(input);
}

function resolveTextDecoder(): TextDecoder {
  return typeof TextDecoder !== "undefined"
    ? new TextDecoder("utf-8")
    : new (require("node:util").TextDecoder)("utf-8");
}

function resolveTextEncoder(): TextEncoder {
  return typeof TextEncoder !== "undefined"
    ? new TextEncoder()
    : new (require("node:util").TextEncoder)();
}

function encodeTargetIds(ids: string[]): {
  packedBytes: Uint8Array;
  offsets: Uint32Array;
  lengths: Uint32Array;
} {
  const encoder = resolveTextEncoder();
  const encoded = ids.map((id) => encoder.encode(id));
  const totalBytes = encoded.reduce((sum, e) => sum + e.length, 0);
  const packedBytes = new Uint8Array(totalBytes);
  const offsets = new Uint32Array(ids.length);
  const lengths = new Uint32Array(ids.length);
  let cursor = 0;
  for (let i = 0; i < ids.length; i++) {
    offsets[i] = cursor >>> 0;
    lengths[i] = encoded[i].length >>> 0;
    packedBytes.set(encoded[i], cursor);
    cursor += encoded[i].length;
  }
  return { packedBytes, offsets, lengths };
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

export function packPeakOptions(opts?: PeakOptions): Uint8Array | undefined {
  if (!opts) return undefined;
  const view = new DataView(new ArrayBuffer(PEAK_OPTIONS_STRUCT_BYTES));
  writeF64(view, 0, opts.integralThreshold);
  writeF64(view, 8, opts.intensityThreshold);
  writeI32(view, 16, opts.widthThreshold);
  view.setInt32(20, 0, true);
  writeF64(view, 24, opts.noise);
  writeI32(view, 32, opts.autoNoise);
  writeI32(view, 36, opts.autoBaseline);
  writeI32(view, 40, opts.baselineWindow);
  writeI32(view, 44, opts.baselineWindowFactor);
  writeI32(view, 48, opts.allowOverlap);
  writeI32(view, 52, opts.windowSize);
  writeF64(view, 56, opts.snRatio);
  return new Uint8Array(view.buffer);
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
  readonly collectScans: (
    handle: number,
    from: number,
    to: number,
    level: number,
    includeMetadata: number,
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
    this.collectScans = this.resolve(ex, ["collect_scans"]);
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
  private readonly fn: WasmExports;

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

  parseMzMLToFile(data: Uint8Array): MzMlFile {
    const [ptr, len] = this.heap.allocAndWrite(data);
    try {
      const rc = this.fn.parseMzml(ptr, len, this.handleScratchSlot);
      if (rc !== 0) throw new Error("parse_mzml failed with code " + rc);
      return new MzMlFile(this.heap.readU32(this.handleScratchSlot));
    } finally {
      this.fn.free(ptr, len);
    }
  }

  parseBinToFile(data: Uint8Array): MzMlFile {
    const [ptr, len] = this.heap.allocAndWrite(data);
    try {
      const rc = this.fn.parseBin(ptr, len, this.handleScratchSlot);
      if (rc !== 0) throw new Error("parse_bin failed with code " + rc);
      return new MzMlFile(this.heap.readU32(this.handleScratchSlot));
    } finally {
      this.fn.free(ptr, len);
    }
  }

  freeFile(handle: number): void {
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
    fromTo: FromTo,
    ppmTolerance: number,
    mzTolerance: number,
  ): { x: Float64Array; y: Float64Array } {
    const { from, to } = fromTo;
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

  collectScans(
    handle: number,
    from: number,
    to: number,
    level = 1,
    includeMetadata = false,
  ): CentroidScan[] {
    const rc = this.fn.collectScans(
      handle,
      from,
      to,
      level,
      includeMetadata ? 1 : 0,
      this.jsonOutputSlot,
    );
    if (rc !== 0) throw new Error("collect_scans failed with code " + rc);
    return camelizeKeys(
      this.heap.readJsonFromSlot<CentroidScan[]>(this.jsonOutputSlot),
    );
  }

  getPeaksFromEic(
    handle: number,
    targets: Target[],
    fromTo: { from: number; to: number },
    peakOpts?: PeakOptions,
    cores = 1,
  ): FoundFeature[] {
    if (targets.length === 0) return [];
    const { rts, mzs, ranges, ids } = unpackTargets(targets);
    const { packedBytes, offsets, lengths } = encodeTargetIds(ids);
    const allocations = this.writeTargetArrays(
      rts,
      mzs,
      ranges,
      offsets,
      lengths,
      packedBytes,
    );
    const [
      rtPtr,
      rtLen,
      mzPtr,
      mzLen,
      rangePtr,
      rangeLen,
      offPtr,
      offLen,
      lenPtr,
      lenLen,
      idPtr,
      idLen,
    ] = allocations.flat();
    const rc = this.fn.getPeaksFromEic(
      handle,
      rtPtr,
      mzPtr,
      rangePtr,
      offPtr,
      lenPtr,
      idPtr,
      idLen,
      targets.length,
      fromTo.from,
      fromTo.to,
      this.peakOptsPtr(peakOpts),
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
    return this.heap.readJsonFromSlot<FoundFeature[]>(this.jsonOutputSlot);
  }

  getPeaksFromChrom(
    handle: number,
    items: ChromItem[],
    peakOpts?: PeakOptions,
    cores = 1,
  ): ChromPeak[] {
    if (items.length === 0) return [];
    const count = items.length;
    const indices = new Uint32Array(count);
    const rts = new Float64Array(count);
    const windows = new Float64Array(count);
    for (let i = 0; i < count; i++) {
      const item = items[i];
      const idx = Number.isFinite(item.idx)
        ? item.idx!
        : Number.isFinite(item.index)
          ? item.index!
          : -1;
      indices[i] = idx >= 0 ? idx >>> 0 : 0xffffffff;
      rts[i] = +item.rt;
      windows[i] = +(item.window ?? item.range ?? 0);
    }
    const [idxPtr, idxLen] = this.heap.allocAndWrite(toUint8View(indices));
    const [rtPtr, rtLen] = this.heap.allocAndWrite(toUint8View(rts));
    const [winPtr, winLen] = this.heap.allocAndWrite(toUint8View(windows));
    const rc = this.fn.getPeaksFromChrom(
      handle,
      idxPtr,
      rtPtr,
      winPtr,
      count,
      this.peakOptsPtr(peakOpts),
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
    return this.heap.readJsonFromSlot<ChromPeak[]>(this.jsonOutputSlot);
  }

  findPeaks(x: Float64Array, y: Float64Array, opts?: PeakOptions): Peak[] {
    const [xPtr, xLen] = this.heap.allocAndWrite(toUint8View(x));
    const [yPtr, yLen] = this.heap.allocAndWrite(toUint8View(y));
    const rc = this.fn.findPeaks(
      xPtr,
      yPtr,
      x.length,
      this.peakOptsPtr(opts),
      this.jsonOutputSlot,
    );
    this.heap.freeMany([
      [xPtr, xLen],
      [yPtr, yLen],
    ]);
    if (rc !== 0) throw new Error("find_peaks failed with code " + rc);
    return this.heap.readJsonFromSlot<Peak[]>(this.jsonOutputSlot);
  }

  getPeak(
    x: Float64Array,
    y: Float64Array,
    rt: number,
    range: number,
    opts?: PeakOptions,
  ): Peak {
    const [xPtr, xLen] = this.heap.allocAndWrite(toUint8View(x));
    const [yPtr, yLen] = this.heap.allocAndWrite(toUint8View(y));
    const rc = this.fn.getPeak(
      xPtr,
      yPtr,
      x.length,
      rt,
      range,
      this.peakOptsPtr(opts),
      this.jsonOutputSlot,
    );
    this.heap.freeMany([
      [xPtr, xLen],
      [yPtr, yLen],
    ]);
    if (rc !== 0) throw new Error("get_peak failed with code " + rc);
    return this.heap.readJsonFromSlot<Peak>(this.jsonOutputSlot);
  }

  findFeature(
    handle: number,
    targets: Target[],
    options: FindFeaturesOptions & {
      scanEic?: { ppmTolerance?: number; mzTolerance?: number };
      cores?: number;
    } = {},
  ): FoundFeature[] {
    if (targets.length === 0) return [];
    const { rts, mzs, ranges, ids } = unpackTargets(targets, 0.5);
    const { packedBytes, offsets, lengths } = encodeTargetIds(ids);
    const allocations = this.writeTargetArrays(
      rts,
      mzs,
      ranges,
      offsets,
      lengths,
      packedBytes,
    );
    const [
      rtPtr,
      rtLen,
      mzPtr,
      mzLen,
      rangePtr,
      rangeLen,
      offPtr,
      offLen,
      lenPtr,
      lenLen,
      idPtr,
      idLen,
    ] = allocations.flat();
    const scanPpm = options.scanEic?.ppmTolerance ?? 10;
    const scanMz = options.scanEic?.mzTolerance ?? 0.003;
    const eicPpm = options.eic?.ppmTolerance ?? 20;
    const eicMz = options.eic?.mzTolerance ?? 0.005;
    const rc = this.fn.findFeature(
      handle,
      rtPtr,
      mzPtr,
      rangePtr,
      offPtr,
      lenPtr,
      idPtr,
      idLen,
      targets.length,
      options.cores ?? 1,
      scanPpm,
      scanMz,
      eicPpm,
      eicMz,
      this.peakOptsPtr(options.findPeak),
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
    return this.heap.readJsonFromSlot<FoundFeature[]>(this.jsonOutputSlot);
  }

  findNoiseLevel(y: Float32Array): number {
    if (y.length === 0) return Infinity;
    const [yPtr, yLen] = this.heap.allocAndWrite(toUint8View(y));
    const noise = this.fn.findNoiseLevel(yPtr, y.length);
    this.fn.free(yPtr, yLen);
    return noise;
  }

  calculateBaseline(y: Float64Array, opts?: BaselineOptions): Float64Array {
    const [yPtr, yLen] = this.heap.allocAndWrite(toUint8View(y));
    const rc = this.fn.calculateBaseline(
      yPtr,
      y.length,
      opts?.baselineWindow ?? 0,
      opts?.baselineWindowFactor ?? 0,
      this.baselineScratchSlot,
    );
    this.fn.free(yPtr, yLen);
    if (rc !== 0) throw new Error("calculate_baseline failed with code " + rc);
    const { ptr, len } = this.heap.readOutputSlot(this.baselineScratchSlot);
    return new Float64Array(this.heap.copyOutAndFree(ptr, len).buffer.slice(0));
  }

  private peakOptsPtr(opts?: PeakOptions): number {
    if (!opts) return 0;
    this.heap.writeF64At(this.peakOptionsSlot, 0, opts.integralThreshold);
    this.heap.writeF64At(this.peakOptionsSlot, 8, opts.intensityThreshold);
    this.heap.writeI32At(this.peakOptionsSlot, 16, opts.widthThreshold);
    this.heap.writeF64At(this.peakOptionsSlot, 24, opts.noise);
    this.heap.writeI32At(this.peakOptionsSlot, 32, opts.autoNoise);
    this.heap.writeI32At(this.peakOptionsSlot, 36, opts.autoBaseline);
    this.heap.writeI32At(this.peakOptionsSlot, 40, opts.baselineWindow);
    this.heap.writeI32At(this.peakOptionsSlot, 44, opts.baselineWindowFactor);
    this.heap.writeI32At(this.peakOptionsSlot, 48, opts.allowOverlap);
    this.heap.writeI32At(this.peakOptionsSlot, 52, opts.windowSize);
    this.heap.writeF64At(this.peakOptionsSlot, 56, opts.snRatio);
    return this.peakOptionsSlot;
  }

  private writeTargetArrays(
    rts: Float64Array,
    mzs: Float64Array,
    ranges: Float64Array,
    offsets: Uint32Array,
    lengths: Uint32Array,
    idBytes: Uint8Array,
  ): [ptr: number, len: number][] {
    return [
      this.heap.allocAndWrite(toUint8View(rts)),
      this.heap.allocAndWrite(toUint8View(mzs)),
      this.heap.allocAndWrite(toUint8View(ranges)),
      this.heap.allocAndWrite(toUint8View(offsets)),
      this.heap.allocAndWrite(toUint8View(lengths)),
      this.heap.allocAndWrite(idBytes),
    ];
  }
}

function unpackTargets(
  targets: Target[],
  defaultRange = 0,
): {
  rts: Float64Array;
  mzs: Float64Array;
  ranges: Float64Array;
  ids: string[];
} {
  const n = targets.length;
  const rts = new Float64Array(n);
  const mzs = new Float64Array(n);
  const ranges = new Float64Array(n);
  const ids = targets.map((t, i) => {
    rts[i] = +t.rt;
    mzs[i] = +t.mz;
    ranges[i] = t.range ?? defaultRange;
    return t.id ?? "";
  });
  return { rts, mzs, ranges, ids };
}

export class MzMlFile {
  _handle: number;

  constructor(handle: number) {
    this._handle = handle;
  }

  dispose(): void {
    if (this._handle) {
      getApi().freeFile(this._handle);
      this._handle = 0;
    }
  }

  toJson(): any {
    this.assertValid("toJson");
    return camelizeKeys(getApi().fileToJson(this._handle));
  }

  toMzml(): string {
    this.assertValid("toMzml");
    return getApi().fileToMzml(this._handle);
  }

  toBin(compressionLevel = 12, f32Compress = false): Uint8Array {
    this.assertValid("toBin");
    return getApi().fileToBin(
      this._handle,
      compressionLevel,
      f32Compress ? 1 : 0,
    );
  }

  private assertValid(caller: string): void {
    if (!this._handle) throw new Error(`${caller}: file has been disposed`);
  }
}

let wasmApi: WasmApi | null = null;
let initPromise: Promise<void> | null = null;

function getApi(): WasmApi {
  if (!wasmApi)
    throw new Error("WASM not ready — an await is missing before this call");
  return wasmApi;
}

async function loadWasm(): Promise<void> {
  const imports = { env: { js_log: (_ptr: number, _len: number) => {} } };

  if (IS_INLINE_BUILD) {
    if (!__WASM_DATA_URL__)
      throw new Error("Inline build is missing __WASM_DATA_URL__");
    const bytes = await fetch(__WASM_DATA_URL__).then((r) => r.arrayBuffer());
    const { instance } = await WebAssembly.instantiate(bytes, imports);
    wasmApi = new WasmApi(instance);
    return;
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
    wasmApi = new WasmApi(instance);
    return;
  }

  throw new Error(
    "Browser ESM build requires the UMD inline bundle (msutils.js)",
  );
}

function ensureReady(): Promise<void> {
  if (!initPromise) initPromise = loadWasm();
  return initPromise;
}

function autoInit<TArgs extends any[], TReturn>(
  fn: (...args: TArgs) => TReturn,
): (...args: TArgs) => Promise<TReturn> {
  return async (...args) => {
    await ensureReady();
    return fn(...args);
  };
}

export async function init(): Promise<void> {
  await ensureReady();
}

export const parseMzML = autoInit(
  (data: BinaryInput): MzMlFile => getApi().parseMzMLToFile(toUint8(data)),
);

export const parseBin = autoInit(
  (data: BinaryInput): MzMlFile => getApi().parseBinToFile(toUint8(data)),
);

export const binToJson = autoInit((file: MzMlFile): string => {
  if (!(file instanceof MzMlFile) || !file._handle)
    throw new Error("binToJson: invalid file");
  return JSON.stringify(file.toJson());
});

export const convertBinToMzml = autoInit((file: MzMlFile): string => {
  if (!(file instanceof MzMlFile) || !file._handle)
    throw new Error("convertBinToMzml: invalid file");
  return file.toMzml();
});

export const mzmlToBin = autoInit(
  (
    file: MzMlFile,
    options: { level?: number; f32Compress?: boolean } = {},
  ): Uint8Array => {
    if (!(file instanceof MzMlFile) || !file._handle)
      throw new Error("mzmlToBin: invalid file");
    return file.toBin(options.level, options.f32Compress);
  },
);

export const calculateEic = autoInit(
  (
    file: MzMlFile,
    targetMz: number,
    fromTo: FromTo,
    ppmTolerance = 20,
    mzTolerance = 0.005,
  ): { x: Float64Array; y: Float64Array } => {
    if (!(file instanceof MzMlFile) || !file._handle)
      throw new Error("calculateEic: invalid file");
    return getApi().calculateEic(
      file._handle,
      +targetMz,
      fromTo,
      ppmTolerance,
      mzTolerance,
    );
  },
);

export const collectScans = autoInit(
  (
    file: MzMlFile,
    from: number,
    to: number,
    level = 1,
    includeMetadata = false,
  ): CentroidScan[] => {
    if (!(file instanceof MzMlFile) || !file._handle)
      throw new Error("collectScans: invalid file");
    return getApi().collectScans(
      file._handle,
      from,
      to,
      level,
      includeMetadata,
    );
  },
);

export const findPeaks = autoInit(
  (x: Float64Array, y: Float64Array, opts?: PeakOptions): Peak[] =>
    getApi().findPeaks(x, y, opts),
);

export const getPeak = autoInit(
  (
    x: Float64Array,
    y: Float64Array,
    rt: number,
    range: number,
    opts?: PeakOptions,
  ): Peak => getApi().getPeak(x, y, rt, range, opts),
);

export const findNoiseLevel = autoInit(
  (y: Float32Array | ArrayLike<number>): number =>
    getApi().findNoiseLevel(
      y instanceof Float32Array ? y : new Float32Array(y as ArrayLike<number>),
    ),
);

export const calculateBaseline = autoInit(
  (
    y: Float64Array | ArrayLike<number>,
    options?: BaselineOptions,
  ): Float64Array => {
    const y64 =
      y instanceof Float64Array ? y : new Float64Array(y as ArrayLike<number>);
    return getApi().calculateBaseline(y64, options);
  },
);

export const getPeaksFromEic = autoInit(
  (
    file: MzMlFile,
    targets: Target[],
    fromTo: { from: number; to: number },
    options?: PeakOptions,
    cores = 1,
  ): FoundFeature[] => {
    if (!(file instanceof MzMlFile) || !file._handle)
      throw new Error("getPeaksFromEic: invalid file");
    return getApi().getPeaksFromEic(
      file._handle,
      targets,
      fromTo,
      options,
      cores,
    );
  },
);

export const getPeaksFromChrom = autoInit(
  (
    file: MzMlFile,
    items: ChromItem[],
    options?: PeakOptions,
    cores = 1,
  ): ChromPeak[] => {
    if (!(file instanceof MzMlFile) || !file._handle)
      throw new Error("getPeaksFromChrom: invalid file");
    return getApi().getPeaksFromChrom(file._handle, items, options, cores);
  },
);

export const findFeature = autoInit(
  (
    file: MzMlFile,
    targets: Target[],
    options: FindFeaturesOptions & {
      scanEic?: { ppmTolerance?: number; mzTolerance?: number };
      cores?: number;
    } = {},
  ): FoundFeature[] => {
    if (!(file instanceof MzMlFile) || !file._handle)
      throw new Error("findFeature: invalid file");
    return getApi().findFeature(file._handle, targets, options);
  },
);

export async function findFeatures(
  _file: MzMlFile,
  _fromTo: { from: number; to: number },
  _options: FindFeaturesOptions = {},
): Promise<Feature[]> {
  throw new Error("findFeatures is not supported in the WASM build");
}
