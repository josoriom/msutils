/// <reference types="node" />

declare const __INLINE__: boolean;
declare const __WASM_DATA_URL__: string | undefined;
declare var require: any;

import type { Backend, FileHandle, ByteRangeResult } from "./backend";
import type { NoiseLevel, PeakOptions } from "../types/types";
import { parseAndCamelize, toUint8 } from "./shared";
import {
  fetchHeader,
  newRemoteSource,
  prefetchRanges as prefetchRangesForSource,
  type RemoteSource,
} from "./remoteSource";

const IS_INLINE_BUILD = typeof __INLINE__ !== "undefined" && __INLINE__;
const OUTPUT_BUFFER_PAIR_BYTES = 8;
const NOISE_LEVEL_SLOT_BYTES = 16;
const NOISE_LEVEL_INTENSITY_OFFSET = 8;
const PEAK_OPTIONS_STRUCT_BYTES = 80;
const QUANTION_ABI_VERSION = 1;
const ERR_FAST_PATH = 6;

function rcError(name: string, rc: number): Error {
  if (rc === ERR_FAST_PATH) {
    return new Error(
      `${name}: fast EIC path unavailable: this .ion file has no usable spectrum bounds (A3); re-encode it with the current Ionic to use the fast EIC path`,
    );
  }
  return new Error(`${name} failed with code ${rc}`);
}

class RangeSourceRegistry {
  private sources = new Map<number, RemoteSource>();
  private next_id = 1;

  add(source: RemoteSource): number {
    const id = this.next_id++;
    this.sources.set(id, source);
    return id;
  }

  get(id: number): RemoteSource | undefined {
    return this.sources.get(id);
  }

  remove(id: number): void {
    this.sources.delete(id);
  }
}

const range_sources = new RangeSourceRegistry();
let wasm_memory: WebAssembly.Memory | null = null;

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

const SHAPE_CODES: Record<string, number> = { gaussian: 0, emg: 1 };
const DEFAULT_SHAPE_CODE = 1;

function shapeCode(value: unknown): number {
  if (typeof value === "number" && Number.isFinite(value)) {
    return value;
  }
  if (typeof value === "string") {
    const code = SHAPE_CODES[value.trim().toLowerCase()];
    return code === undefined ? DEFAULT_SHAPE_CODE : code;
  }
  return DEFAULT_SHAPE_CODE;
}

function packPeakOptions(opts: PeakOptions): Uint8Array {
  const view = new DataView(new ArrayBuffer(PEAK_OPTIONS_STRUCT_BYTES));
  writeF64(view, 0, opts.minIntegral);
  writeF64(view, 8, opts.minIntensity);
  writeI32(view, 16, opts.minPeakWidthPoints);
  view.setInt32(20, shapeCode(opts.shape), true);
  writeF64(view, 24, opts.noise);
  writeI32(view, 32, opts.autoNoise);
  writeI32(view, 36, opts.autoBaseline);
  writeI32(view, 40, opts.lambda);
  writeI32(view, 44, opts.maxIterations);
  writeI32(view, 48, opts.allowOverlap);
  view.setInt32(52, 0, true);
  writeF64(view, 56, opts.minSnr);
  writeF64(view, 64, opts.minR2);
  writeI32(view, 72, opts.kernelSize);
  return new Uint8Array(view.buffer);
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

  readF64(ptr: number): number {
    this.refresh();
    return this.view.getFloat64(ptr, true);
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
    return parseAndCamelize(this.decoder.decode(bytes)) as T;
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
  readonly quantionAbiVersion: () => number;
  readonly quantionSizeofPeakOptions: () => number;
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
  readonly findNoiseLevel: (
    yPtr: number,
    len: number,
    outWidth: number,
    outIntensity: number,
  ) => number;
  readonly calculateBaseline: (
    yPtr: number,
    len: number,
    lambda: number,
    maxIterations: number,
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
  readonly getIonImage: (
    handle: number,
    mz: number,
    tolerance: number,
    level: number,
    outBuf: number,
  ) => number;
  readonly parseIonUrl:
    | ((sourceId: number, cacheBytes: number, outHandle: number) => number)
    | null;
  readonly planOpen:
    | ((headerPtr: number, headerLen: number, outBuf: number) => number)
    | null;
  readonly planEic:
    | ((
        handle: number,
        target: number,
        from: number,
        to: number,
        ppm: number,
        mzTol: number,
        outBuf: number,
      ) => number)
    | null;
  readonly imageBegin:
    | ((
        handle: number,
        mz: number,
        tolerance: number,
        level: number,
        outSession: number,
      ) => number)
    | null;
  readonly imageScanCount:
    | ((session: number, outCount: number) => number)
    | null;
  readonly imageRanges:
    | ((
        handle: number,
        session: number,
        from: number,
        count: number,
        outBuf: number,
      ) => number)
    | null;
  readonly imageFold:
    | ((handle: number, session: number, from: number, count: number) => number)
    | null;
  readonly imageFinish: ((session: number, outBuf: number) => number) | null;
  readonly imageFree: ((session: number) => void) | null;

  constructor(instance: WebAssembly.Instance) {
    const ex = instance.exports;
    this.memory = ex.memory as WebAssembly.Memory;
    this.alloc = this.resolve(ex, ["alloc"]);
    this.free = this.resolve(ex, ["free_", "free"]);
    this.quantionAbiVersion = this.resolve(ex, ["quantion_abi_version"]);
    this.quantionSizeofPeakOptions = this.resolve(ex, [
      "quantion_sizeof_peak_options",
    ]);
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
    this.getIonImage = this.resolve(ex, ["get_ion_image"]);
    this.parseIonUrl =
      typeof ex["parse_ion_url"] === "function"
        ? (ex["parse_ion_url"] as unknown as (
            sourceId: number,
            cacheBytes: number,
            outHandle: number,
          ) => number)
        : null;
    this.planOpen =
      typeof ex["plan_open"] === "function"
        ? (ex["plan_open"] as unknown as (
            headerPtr: number,
            headerLen: number,
            outBuf: number,
          ) => number)
        : null;
    this.planEic =
      typeof ex["plan_eic"] === "function"
        ? (ex["plan_eic"] as unknown as (
            handle: number,
            target: number,
            from: number,
            to: number,
            ppm: number,
            mzTol: number,
            outBuf: number,
          ) => number)
        : null;
    this.imageBegin =
      typeof ex["image_begin"] === "function"
        ? (ex["image_begin"] as unknown as (
            handle: number,
            mz: number,
            tolerance: number,
            level: number,
            outSession: number,
          ) => number)
        : null;
    this.imageScanCount =
      typeof ex["image_scan_count"] === "function"
        ? (ex["image_scan_count"] as unknown as (
            session: number,
            outCount: number,
          ) => number)
        : null;
    this.imageRanges =
      typeof ex["image_ranges"] === "function"
        ? (ex["image_ranges"] as unknown as (
            handle: number,
            session: number,
            from: number,
            count: number,
            outBuf: number,
          ) => number)
        : null;
    this.imageFold =
      typeof ex["image_fold"] === "function"
        ? (ex["image_fold"] as unknown as (
            handle: number,
            session: number,
            from: number,
            count: number,
          ) => number)
        : null;
    this.imageFinish =
      typeof ex["image_finish"] === "function"
        ? (ex["image_finish"] as unknown as (
            session: number,
            outBuf: number,
          ) => number)
        : null;
    this.imageFree =
      typeof ex["image_free"] === "function"
        ? (ex["image_free"] as unknown as (session: number) => void)
        : null;
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

const IMAGE_BATCH_SCANS = 256;
const IMAGE_BATCH_BUDGET_BYTES = 192 * 1024 * 1024;

function total_range_bytes(ranges: ByteRangeResult[]): bigint {
  let total = 0n;
  for (const range of ranges) total += range.length;
  return total;
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
  private readonly noiseLevelSlot: number;
  private readonly source_id_by_handle = new Map<number, number>();

  constructor(instance: WebAssembly.Instance) {
    this.fn = new WasmExports(instance);

    const abi = this.fn.quantionAbiVersion();
    if (abi !== QUANTION_ABI_VERSION) {
      throw new Error(
        `quantion WASM ABI version mismatch: binary=${abi}, ` +
          `wrapper=${QUANTION_ABI_VERSION}. Rebuild the WASM module.`,
      );
    }
    const size = this.fn.quantionSizeofPeakOptions();
    if (size !== PEAK_OPTIONS_STRUCT_BYTES) {
      throw new Error(
        `quantion WASM PeakOptions size mismatch: binary=${size} bytes, ` +
          `wrapper=${PEAK_OPTIONS_STRUCT_BYTES} bytes. ` +
          `Native binary and JS wrapper are out of sync — rebuild.`,
      );
    }

    this.heap = new WasmHeap(this.fn.memory, this.fn.alloc, this.fn.free);

    this.jsonOutputSlot = this.heap.allocPermanent(OUTPUT_BUFFER_PAIR_BYTES);
    this.jsonScratchSlot = this.heap.allocPermanent(OUTPUT_BUFFER_PAIR_BYTES);
    this.blobScratchSlot = this.heap.allocPermanent(OUTPUT_BUFFER_PAIR_BYTES);
    this.baselineScratchSlot = this.heap.allocPermanent(
      OUTPUT_BUFFER_PAIR_BYTES,
    );
    this.peakOptionsSlot = this.heap.allocPermanent(PEAK_OPTIONS_STRUCT_BYTES);
    this.handleScratchSlot = this.heap.allocPermanent(4);
    this.noiseLevelSlot = this.heap.allocPermanent(NOISE_LEVEL_SLOT_BYTES);
  }

  private assertLocalHandle(handle: number, functionName: string): void {
    if (this.source_id_by_handle.has(handle)) {
      throw new Error(
        `${functionName}: remote browser Ion handles need plan_peaks_from_eic; use Node, Python, R, or a local/in-memory file`,
      );
    }
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

  private readRangesFromSlot(slotPtr: number): ByteRangeResult[] {
    const { ptr, len } = this.heap.readOutputSlot(slotPtr);
    const bytes = this.heap.copyOutAndFree(ptr, len);

    if (bytes.length % 16 !== 0) {
      throw new Error(
        `planned range buffer length must be divisible by 16, got ${bytes.length}`,
      );
    }

    const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
    const ranges: ByteRangeResult[] = [];

    for (let offset = 0; offset < bytes.length; offset += 16) {
      ranges.push({
        offset: view.getBigUint64(offset, true),
        length: view.getBigUint64(offset + 8, true),
      });
    }

    return ranges;
  }

  planOpen(header: Uint8Array): ByteRangeResult[] {
    if (!this.fn.planOpen) {
      throw new Error("planOpen not in WASM exports");
    }

    const [ptr, len] = this.heap.allocAndWrite(header);
    try {
      const rc = this.fn.planOpen(ptr, len, this.blobScratchSlot);
      if (rc !== 0) throw new Error(`plan_open failed with code ${rc}`);
      return this.readRangesFromSlot(this.blobScratchSlot);
    } finally {
      this.fn.free(ptr, len);
    }
  }

  planEic(
    handle: number,
    targetMz: number,
    from: number,
    to: number,
    ppmTol: number,
    mzTol: number,
  ): ByteRangeResult[] {
    if (!this.fn.planEic) {
      throw new Error("planEic not in WASM exports");
    }

    const rc = this.fn.planEic(
      handle,
      targetMz,
      from,
      to,
      ppmTol,
      mzTol,
      this.blobScratchSlot,
    );
    if (rc !== 0) throw new Error(`plan_eic failed with code ${rc}`);
    return this.readRangesFromSlot(this.blobScratchSlot);
  }

  private imageBegin(
    handle: number,
    mz: number,
    tolerance: number,
    level: number,
  ): number {
    if (!this.fn.imageBegin) throw new Error("imageBegin not in WASM exports");
    const rc = this.fn.imageBegin(
      handle,
      mz,
      tolerance,
      level,
      this.handleScratchSlot,
    );
    if (rc !== 0) throw new Error(`image_begin failed with code ${rc}`);
    return this.heap.readU32(this.handleScratchSlot);
  }

  private imageScanCount(session: number): number {
    if (!this.fn.imageScanCount)
      throw new Error("imageScanCount not in WASM exports");
    const rc = this.fn.imageScanCount(session, this.handleScratchSlot);
    if (rc !== 0) throw new Error(`image_scan_count failed with code ${rc}`);
    return this.heap.readU32(this.handleScratchSlot);
  }

  private imageRanges(
    handle: number,
    session: number,
    from: number,
    count: number,
  ): ByteRangeResult[] {
    if (!this.fn.imageRanges)
      throw new Error("imageRanges not in WASM exports");
    const rc = this.fn.imageRanges(
      handle,
      session,
      from,
      count,
      this.blobScratchSlot,
    );
    if (rc !== 0) throw new Error(`image_ranges failed with code ${rc}`);
    return this.readRangesFromSlot(this.blobScratchSlot);
  }

  private imageFold(
    handle: number,
    session: number,
    from: number,
    count: number,
  ): void {
    if (!this.fn.imageFold) throw new Error("imageFold not in WASM exports");
    const rc = this.fn.imageFold(handle, session, from, count);
    if (rc !== 0) throw new Error(`image_fold failed with code ${rc}`);
  }

  private imageFinish(session: number): any {
    if (!this.fn.imageFinish)
      throw new Error("imageFinish not in WASM exports");
    const rc = this.fn.imageFinish(session, this.jsonOutputSlot);
    if (rc !== 0) throw new Error(`image_finish failed with code ${rc}`);
    return this.heap.readJsonFromSlot<any>(this.jsonOutputSlot);
  }

  private imageFree(session: number): void {
    if (this.fn.imageFree) this.fn.imageFree(session);
  }

  private async computeImageStreaming(
    source: RemoteSource,
    handle: number,
    mz: number,
    tolerance: number,
    level: number,
    onProgress?: (fetched: number, total: number, heldBytes: number) => void,
  ): Promise<any> {
    const session = this.imageBegin(handle, mz, tolerance, level);
    try {
      const total = this.imageScanCount(session);
      onProgress?.(0, total, source.cache.bytes());
      source.cache.pin();

      let done = 0;
      while (done < total) {
        let count = Math.min(IMAGE_BATCH_SCANS, total - done);
        let ranges = this.imageRanges(handle, session, done, count);
        while (
          count > 1 &&
          total_range_bytes(ranges) > BigInt(IMAGE_BATCH_BUDGET_BYTES)
        ) {
          count = Math.floor(count / 2);
          ranges = this.imageRanges(handle, session, done, count);
        }
        await prefetchRangesForSource(source, ranges);
        this.imageFold(handle, session, done, count);
        source.cache.dropUnpinned();
        done += count;
        onProgress?.(done, total, source.cache.bytes());
      }

      return this.imageFinish(session);
    } finally {
      this.imageFree(session);
    }
  }

  async parseIonUrl(url: string, cacheSize: number): Promise<number> {
    if (!this.fn.parseIonUrl) {
      throw new Error("parse_ion_url not in WASM exports");
    }

    const source = newRemoteSource(url);
    const sourceId = range_sources.add(source);

    try {
      const headerBytes = await fetchHeader(source);
      await prefetchRangesForSource(source, this.planOpen(headerBytes));

      const rc = this.fn.parseIonUrl(
        sourceId,
        cacheSize,
        this.handleScratchSlot,
      );
      if (rc !== 0) {
        throw new Error(`parse_ion_url failed with code ${rc}`);
      }

      const handle = this.heap.readU32(this.handleScratchSlot);
      if (handle === 0) {
        throw new Error("parse_ion_url returned a null handle");
      }

      this.source_id_by_handle.set(handle, sourceId);
      return handle;
    } catch (error) {
      range_sources.remove(sourceId);
      throw error;
    }
  }

  freeRaw(handle: number): void {
    try {
      this.fn.freeMzml(handle);
    } finally {
      const source_id = this.source_id_by_handle.get(handle);
      if (source_id !== undefined) {
        range_sources.remove(source_id);
        this.source_id_by_handle.delete(handle);
      }
    }
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

  private calculateEicNow(
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
    if (rc !== 0) throw rcError("calculate_eic", rc);
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

  async calculateEic(
    handle: number,
    targetMz: number,
    from: number,
    to: number,
    ppmTolerance: number,
    mzTolerance: number,
  ): Promise<{ x: Float64Array; y: Float64Array }> {
    const sourceId = this.source_id_by_handle.get(handle);

    if (sourceId !== undefined) {
      const source = range_sources.get(sourceId);
      if (!source) {
        throw new Error("calculateEic: remote source is missing");
      }

      const ranges = this.planEic(
        handle,
        targetMz,
        from,
        to,
        ppmTolerance,
        mzTolerance,
      );
      await prefetchRangesForSource(source, ranges);
    }

    return this.calculateEicNow(
      handle,
      targetMz,
      from,
      to,
      ppmTolerance,
      mzTolerance,
    );
  }

  getScans(
    handle: number,
    queryType: number,
    a: number,
    b: number,
    level: number,
  ): any {
    this.assertLocalHandle(handle, "getScans");
    const rc = this.fn.getScans(
      handle,
      queryType,
      a,
      b,
      level,
      this.jsonOutputSlot,
    );
    if (rc !== 0) throw new Error("get_scans failed with code " + rc);
    return this.heap.readJsonFromSlot<any>(this.jsonOutputSlot);
  }

  private getIonImageNow(
    handle: number,
    mz: number,
    tolerance: number,
    level: number,
  ): any {
    const rc = this.fn.getIonImage(
      handle,
      mz,
      tolerance,
      level,
      this.jsonOutputSlot,
    );
    if (rc !== 0) throw new Error("get_ion_image failed with code " + rc);
    return this.heap.readJsonFromSlot<any>(this.jsonOutputSlot);
  }

  async getIonImage(
    handle: number,
    mz: number,
    tolerance: number,
    level: number,
    onProgress?: (fetched: number, total: number, heldBytes: number) => void,
  ): Promise<any> {
    const sourceId = this.source_id_by_handle.get(handle);

    if (sourceId !== undefined) {
      const source = range_sources.get(sourceId);
      if (!source) {
        throw new Error("getIonImage: remote source is missing");
      }
      return this.computeImageStreaming(
        source,
        handle,
        mz,
        tolerance,
        level,
        onProgress,
      );
    }

    return this.getIonImageNow(handle, mz, tolerance, level);
  }

  findPeaks(
    x: Float64Array,
    y: Float64Array,
    opts: PeakOptions | undefined,
  ): any {
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
    return this.heap.readJsonFromSlot<any>(this.jsonOutputSlot);
  }

  getPeak(
    x: Float64Array,
    y: Float64Array,
    rt: number,
    range: number,
    opts: PeakOptions | undefined,
  ): any {
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
    return this.heap.readJsonFromSlot<any>(this.jsonOutputSlot);
  }

  findNoiseLevel(y: Float32Array): NoiseLevel {
    const [yPtr, yLen] = this.heap.allocAndWrite(toUint8View(y));
    const rc = this.fn.findNoiseLevel(
      yPtr,
      y.length,
      this.noiseLevelSlot,
      this.noiseLevelSlot + NOISE_LEVEL_INTENSITY_OFFSET,
    );
    this.fn.free(yPtr, yLen);
    if (rc !== 0) throw new Error("find_noise_level failed with code " + rc);
    return {
      width: this.heap.readU32(this.noiseLevelSlot),
      intensity: this.heap.readF64(
        this.noiseLevelSlot + NOISE_LEVEL_INTENSITY_OFFSET,
      ),
    };
  }

  calculateBaseline(
    y: Float64Array,
    lambda: number,
    maxIterations: number,
  ): Float64Array {
    const [yPtr, yLen] = this.heap.allocAndWrite(toUint8View(y));
    const rc = this.fn.calculateBaseline(
      yPtr,
      y.length,
      lambda,
      maxIterations,
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
    opts: PeakOptions | undefined,
    cores: number,
  ): any {
    this.assertLocalHandle(handle, "getPeaksFromEic");
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
      this.peakOptsPtr(opts),
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
    if (rc !== 0) throw rcError("get_peaks_from_eic", rc);
    return this.heap.readJsonFromSlot<any>(this.jsonOutputSlot);
  }

  getPeaksFromChrom(
    handle: number,
    indices: Uint32Array,
    rts: Float64Array,
    windows: Float64Array,
    count: number,
    opts: PeakOptions | undefined,
    cores: number,
  ): any {
    this.assertLocalHandle(handle, "getPeaksFromChrom");
    const [idxPtr, idxLen] = this.heap.allocAndWrite(toUint8View(indices));
    const [rtPtr, rtLen] = this.heap.allocAndWrite(toUint8View(rts));
    const [winPtr, winLen] = this.heap.allocAndWrite(toUint8View(windows));
    const rc = this.fn.getPeaksFromChrom(
      handle,
      idxPtr,
      rtPtr,
      winPtr,
      count,
      this.peakOptsPtr(opts),
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
    opts: PeakOptions | undefined,
  ): any {
    this.assertLocalHandle(handle, "findFeature");
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
      this.peakOptsPtr(opts),
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

  private peakOptsPtr(opts: PeakOptions | undefined): number {
    if (!opts) return 0;
    const packed = packPeakOptions(opts);
    this.heap.writeBytesAt(this.peakOptionsSlot, packed);
    return this.peakOptionsSlot;
  }
}

function make_range_read_import(): (
  source_id: number,
  offset_lo: number,
  offset_hi: number,
  len: number,
  dest_ptr: number,
) => number {
  return (
    source_id: number,
    offset_lo: number,
    offset_hi: number,
    len: number,
    dest_ptr: number,
  ): number => {
    if (!wasm_memory) return -1;

    const source = range_sources.get(source_id >>> 0);
    if (!source) return -1;

    const offset = BigInt(offset_lo >>> 0) + (BigInt(offset_hi >>> 0) << 32n);
    const length = BigInt(len >>> 0);

    const bytes = source.cache.read(offset, length);
    if (!bytes) return -1;

    new Uint8Array(wasm_memory.buffer, dest_ptr >>> 0, len >>> 0).set(bytes);
    return 0;
  };
}

async function load_and_instantiate(
  bytes: ArrayBuffer | Uint8Array,
): Promise<WasmApi> {
  const imports = {
    env: {
      js_log: (_ptr: number, _len: number) => {},
      range_read: make_range_read_import(),
    },
  };
  const result = await WebAssembly.instantiate(bytes, imports);
  const instance = (result as any).instance as WebAssembly.Instance;
  wasm_memory = (instance.exports as any).memory as WebAssembly.Memory;
  return new WasmApi(instance);
}

async function loadWasm(): Promise<WasmApi> {
  if (IS_INLINE_BUILD) {
    if (!__WASM_DATA_URL__)
      throw new Error("Inline build is missing __WASM_DATA_URL__");
    const bytes = await fetch(__WASM_DATA_URL__).then((r) => r.arrayBuffer());
    return load_and_instantiate(bytes);
  }

  if (typeof process !== "undefined" && (process as any).versions?.node) {
    const nfs = require("node:fs");
    const npath = require("node:path");
    const nurl = require("node:url");
    const dir =
      typeof __dirname !== "undefined"
        ? __dirname
        : npath.dirname(nurl.fileURLToPath((0, eval)("import.meta").url));
    const bytes = nfs.readFileSync(npath.join(dir, "quantion.wasm"));
    return load_and_instantiate(bytes);
  }

  throw new Error(
    "Browser ESM build requires the UMD inline bundle (quantion.js)",
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

  parseIonPath(_path: string, _cacheSize = 0): FileHandle {
    throw new Error(
      "parseIon: file paths need the Node backend; pass a URL or a buffer in the browser",
    );
  }

  async parseIonRemote(url: URL, cacheSize = 0): Promise<FileHandle> {
    return this.getApi().parseIonUrl(url.href, cacheSize);
  }

  parseIonBuffer(bytes: Uint8Array, cacheSize = 0): FileHandle {
    return this.getApi().parseBinRaw(bytes, cacheSize);
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

  async calculateEic(
    handle: FileHandle,
    targetMz: number,
    from: number,
    to: number,
    ppmTol: number,
    mzTol: number,
  ): Promise<{ x: Float64Array; y: Float64Array }> {
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

  async getIonImage(
    handle: FileHandle,
    mz: number,
    tolerance: number,
    level: number,
    onProgress?: (fetched: number, total: number, heldBytes: number) => void,
  ): Promise<any> {
    return this.getApi().getIonImage(
      handle as number,
      mz,
      tolerance,
      level,
      onProgress,
    );
  }

  findPeaks(
    x: Float64Array,
    y: Float64Array,
    opts: PeakOptions | undefined,
  ): any {
    return this.getApi().findPeaks(x, y, opts);
  }

  getPeak(
    x: Float64Array,
    y: Float64Array,
    rt: number,
    range: number,
    opts: PeakOptions | undefined,
  ): any {
    return this.getApi().getPeak(x, y, rt, range, opts);
  }

  findNoiseLevel(y: Float32Array): NoiseLevel {
    const y32 = y instanceof Float32Array ? y : new Float32Array(y);
    return this.getApi().findNoiseLevel(y32);
  }

  calculateBaseline(
    y: Float64Array,
    lambda: number,
    maxIterations: number,
  ): Float64Array {
    return this.getApi().calculateBaseline(y, lambda, maxIterations);
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
    opts: PeakOptions | undefined,
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
      opts,
      cores,
    );
  }

  getPeaksFromChrom(
    handle: FileHandle,
    indices: Uint32Array,
    rts: Float64Array,
    windows: Float64Array,
    count: number,
    opts: PeakOptions | undefined,
    cores: number,
  ): any {
    return this.getApi().getPeaksFromChrom(
      handle as number,
      indices,
      rts,
      windows,
      count,
      opts,
      cores,
    );
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
    opts: PeakOptions | undefined,
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
      opts,
    );
  }
}

