import { makeParseMzML, ParseMzML } from "./utilities/parseMzML.js";
import { calculateEic } from "./utilities/calculateEic.js";

export type FindPeaksOptions = {
  integralThreshold?: number;
  intensityThreshold?: number;
  widthThreshold?: number;
  noise?: number;
  autoNoise?: boolean;
  allowOverlap?: boolean;
  windowSize?: number;
  snRatio?: number;
  autoBaseline?: boolean;
  baselineWindow?: number;
  baselineWindowFactor?: number;
};

export type BaselineOptions = {
  baselineWindow?: number;
  baselineWindowFactor?: number;
};

export type Peak = {
  from: number;
  to: number;
  rt: number;
  integral: number;
  intensity: number;
  ratio: number;
  np: number;
};

export type ChromPeakRow = {
  index: number;
  id: string;
  ort: number;
  rt: number;
  from: number;
  to: number;
  intensity: number;
  integral: number;
};

export type EicPeakRow = {
  id: string;
  mz: number;
  ort: number;
  rt: number;
  from: number;
  to: number;
  intensity: number;
  integral: number;
};

type EicItem = {
  id?: string;
  rt: number;
  mz: number;
  range?: number;
  ranges?: number;
  sd?: number;
  window?: number;
};

export type Target = {
  id?: string;
  rt: number;
  mz: number;
  range?: number;
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

export type FindFeaturesOptions = {
  scanEic?: { ppmTolerance?: number; mzTolerance?: number };
  eic?: { ppmTolerance?: number; mzTolerance?: number };
  findPeak?: FindPeaksOptions;
};

export interface Exports {
  parseMzML: ParseMzML;
  calculateEic: typeof calculateEic;
  getPeak: (
    x: Float64Array,
    y: Float64Array,
    rt: number,
    range: number,
    options?: FindPeaksOptions
  ) => Peak;
  getPeaksFromChrom: (
    bin: Uint8Array,
    items: { idx: number; rt: number; window: number }[],
    options?: FindPeaksOptions,
    cores?: number
  ) => ChromPeakRow[];
  getPeaksFromEic: (
    bin: Uint8Array,
    items: EicItem[],
    fromTo: { from: number; to: number },
    options?: FindPeaksOptions,
    cores?: number
  ) => EicPeakRow[];
  findPeaks: (
    x: Float64Array,
    y: Float64Array,
    options?: FindPeaksOptions
  ) => Peak[];
  findFeature: (
    bin: Uint8Array,
    targets: Target[],
    options?: FindFeaturesOptions
  ) => FoundFeature[];
  findNoiseLevel: (y: Float32Array | ArrayLike<number>) => number;
  calculateBaseline: (
    y: Float64Array,
    options?: BaselineOptions
  ) => Float64Array;
  collectMs1Scans: (
    bin: Uint8Array,
    fromTo: { from: number; to: number }
  ) => {
    rt: Float64Array;
    offsets: Uint32Array;
    lengths: Uint32Array;
    mz: Float64Array;
    intensity: Float64Array;
  };
  __debug: {
    memory: WebAssembly.Memory;
    exports: Record<string, any>;
    heapBytes: () => number;
  };
}

type ExportsLike = Record<string, any>;
const getExports = (obj: any): ExportsLike => {
  if (obj?.exports?.memory instanceof WebAssembly.Memory) return obj.exports;
  if (obj?.asm?.memory instanceof WebAssembly.Memory) return obj.asm;
  if (obj?.memory instanceof WebAssembly.Memory) return obj;
  throw new Error("makeApi: could not find WebAssembly exports");
};
function pickFn<T extends Function>(ex: ExportsLike, names: string[]): T {
  for (let i = 0; i < names.length; i++) {
    const f = ex[names[i]];
    if (typeof f === "function") return f as unknown as T;
  }
  throw new Error("makeApi: missing exports: " + names.join("|"));
}

const BUF_PAIR_BYTES = 8;
const SIZE_COPTS = 64;

export function makeApi(instanceOrModule: any): Exports {
  const ex = getExports(instanceOrModule);
  const memory: WebAssembly.Memory = ex.memory;

  const alloc: (n: number) => number = pickFn(ex, ["alloc"]);
  const free: (p: number, n: number) => void = pickFn(ex, ["free_", "free"]);

  const _parse_mzml_raw: (p: number, n: number, outBinBuf: number) => number =
    pickFn(ex, ["parse_mzml"]);

  const _parse_mzml_to_json_raw: (
    p: number,
    n: number,
    outJsonBuf: number,
    outBlobBuf: number
  ) => number = pickFn(ex, ["parse_mzml_to_json"]);

  const find_peaks: (
    xPtr: number,
    yPtr: number,
    len: number,
    optionsPtr: number,
    outJsonBuf: number
  ) => number = pickFn(ex, ["find_peaks"]);

  const get_peak: (
    xPtr: number,
    yPtr: number,
    len: number,
    rt: number,
    range: number,
    optionsPtr: number,
    outJsonBuf: number
  ) => number = pickFn(ex, ["get_peak"]);

  const get_peaks_from_chrom: (
    binPtr: number,
    binLen: number,
    idxsPtr: number,
    rtsPtr: number,
    rangesPtr: number,
    nItems: number,
    optionsPtr: number,
    cores: number,
    outJsonBuf: number
  ) => number = pickFn(ex, ["get_peaks_from_chrom"]);

  const get_peaks_from_eic: (
    binPtr: number,
    binLen: number,
    rtsPtr: number,
    mzsPtr: number,
    rangesPtr: number,
    idsOffPtr: number,
    idsLenPtr: number,
    idsBufPtr: number,
    idsBufLen: number,
    nItems: number,
    fromLeft: number,
    toRight: number,
    optionsPtr: number,
    cores: number,
    outJsonBuf: number
  ) => number = pickFn(ex, ["get_peaks_from_eic"]);

  const find_noise_level: (yPtr: number, len: number) => number = pickFn(ex, [
    "find_noise_level",
  ]);

  const calculate_baseline: (
    yPtr: number,
    len: number,
    baselineWindow: number,
    baselineWindowFactor: number,
    outBuf: number
  ) => number = pickFn(ex, ["calculate_baseline"]);

  const find_feature: (
    binPtr: number,
    binLen: number,
    rtsPtr: number,
    mzsPtr: number,
    rangesPtr: number,
    idsOffPtr: number,
    idsLenPtr: number,
    idsBufPtr: number,
    idsBufLen: number,
    nItems: number,
    scanPpm: number,
    scanMz: number,
    eicPpm: number,
    eicMz: number,
    optionsPtr: number,
    cores: number,
    outJsonBuf: number
  ) => number = pickFn(ex, ["find_feature"]);

  const collect_ms1_scans: (
    binPtr: number,
    binLen: number,
    fromTime: number,
    toTime: number,
    outRtBufPtr: number,
    outOffsetsBufPtr: number,
    outLengthsBufPtr: number,
    outMzBufPtr: number,
    outIntensityBufPtr: number
  ) => number = pickFn(ex, ["collect_ms1_scans"]);

  let HEAPU8 = new Uint8Array(memory.buffer);
  let HEAPDV = new DataView(memory.buffer);
  const refreshViews = () => {
    if (HEAPU8.buffer !== memory.buffer) {
      HEAPU8 = new Uint8Array(memory.buffer);
      HEAPDV = new DataView(memory.buffer);
    }
  };
  const readBuf = (bufPtr: number) => {
    refreshViews();
    const ptr = HEAPDV.getUint32(bufPtr + 0, true);
    const len = HEAPDV.getUint32(bufPtr + 4, true);
    return { ptr, len };
  };
  const heapWrite = (dstPtr: number, src: Uint8Array) => {
    refreshViews();
    HEAPU8.set(src, dstPtr);
  };
  const heapSlice = (ptr: number, len: number): Uint8Array => {
    refreshViews();
    const out = new Uint8Array(len);
    out.set(HEAPU8.subarray(ptr, ptr + len));
    return out;
  };

  const SCRATCH_A = alloc(BUF_PAIR_BYTES);
  const SCRATCH_JSON = alloc(BUF_PAIR_BYTES);
  const SCRATCH_BLOB = alloc(BUF_PAIR_BYTES);
  const SCRATCH_PEAKS = alloc(BUF_PAIR_BYTES);
  const SCRATCH_OPTS = alloc(SIZE_COPTS);
  const SCRATCH_BASE = alloc(BUF_PAIR_BYTES);

  refreshViews();

  const td: TextDecoder =
    typeof TextDecoder !== "undefined"
      ? new TextDecoder("utf-8")
      : new (require("node:util").TextDecoder)("utf-8");
  const te: TextEncoder =
    typeof TextEncoder !== "undefined"
      ? new TextEncoder()
      : new (require("node:util").TextEncoder)();

  const writeOptions = (ptr: number, o?: FindPeaksOptions) => {
    refreshViews();
    const setI32 = (off: number, v: number) =>
      HEAPDV.setInt32(ptr + off, v | 0, true);
    const setF64 = (off: number, v: number) =>
      HEAPDV.setFloat64(ptr + off, v, true);
    if (!o) {
      setF64(0, NaN);
      setF64(8, NaN);
      setI32(16, 0);
      setF64(24, NaN);
      setI32(32, 0);
      setI32(36, 0);
      setI32(40, 0);
      setI32(44, 0);
      setI32(48, 0);
      setI32(52, 0);
      setF64(56, NaN);
      return;
    }
    setF64(
      0,
      typeof o.integralThreshold === "number" ? o.integralThreshold : NaN
    );
    setF64(
      8,
      typeof o.intensityThreshold === "number" ? o.intensityThreshold : NaN
    );
    setI32(16, typeof o.widthThreshold === "number" ? o.widthThreshold : 0);
    setF64(24, typeof o.noise === "number" ? o.noise : NaN);
    setI32(32, o.autoNoise ? 1 : 0);
    setI32(36, o.autoBaseline ? 1 : 0);
    setI32(40, typeof o.baselineWindow === "number" ? o.baselineWindow : 0);
    setI32(
      44,
      typeof o.baselineWindowFactor === "number" ? o.baselineWindowFactor : 0
    );
    setI32(48, o.allowOverlap ? 1 : 0);
    setI32(52, typeof o.windowSize === "number" ? o.windowSize : 0);
    setF64(56, typeof o.snRatio === "number" ? o.snRatio : NaN);
  };

  const parse_mzml = (p: number, n: number, _slim: number, outBinBuf: number) =>
    _parse_mzml_raw(p, n, outBinBuf);

  const parse_mzml_to_json = (
    p: number,
    n: number,
    _slim: number,
    outJsonBuf: number,
    outBlobBuf: number
  ) => _parse_mzml_to_json_raw(p, n, outJsonBuf, outBlobBuf);

  const parseMzML = makeParseMzML({
    alloc,
    free,
    refreshViews,
    heapWrite,
    readBuf,
    heapSlice,
    parse_mzml,
    SCRATCH_A,
    parse_mzml_to_json,
    SCRATCH_JSON,
    SCRATCH_BLOB,
    td,
  });

  const getPeak = (
    x: Float64Array,
    y: Float64Array,
    rt: number,
    range: number,
    options?: FindPeaksOptions
  ): Peak => {
    if (!(x instanceof Float64Array))
      throw new Error("getPeak: x must be Float64Array");
    if (!(y instanceof Float64Array))
      throw new Error("getPeak: y must be Float64Array");
    if (x.length !== y.length || x.length < 3)
      throw new Error("getPeak: x,y length mismatch (>=3)");
    const xb = new Uint8Array(x.buffer, x.byteOffset, x.byteLength);
    const yb = new Uint8Array(y.buffer, y.byteOffset, y.byteLength);
    const xp = alloc(xb.length);
    const yp = alloc(yb.length);
    heapWrite(xp, xb);
    heapWrite(yp, yb);
    let pOpts = 0;
    if (options) {
      writeOptions(SCRATCH_OPTS, options);
      pOpts = SCRATCH_OPTS;
    }
    const rc = get_peak(xp, yp, x.length, rt, range, pOpts, SCRATCH_JSON);
    free(xp, xb.length);
    free(yp, yb.length);
    if (rc !== 0) throw new Error(`get_peak failed: ${rc}`);
    const { ptr, len } = readBuf(SCRATCH_JSON);
    const bytes = heapSlice(ptr, len);
    free(ptr, len);
    return JSON.parse(td.decode(bytes));
  };

  const getPeaksFromChrom = (
    bin: Uint8Array,
    items: { idx: number; rt: number; window: number }[],
    options?: FindPeaksOptions,
    cores = 1
  ): ChromPeakRow[] => {
    const n = items.length;
    if (n === 0) return [];
    const idxs = new Uint32Array(n);
    const rts = new Float64Array(n);
    const wins = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const it = items[i];
      const idx = it.idx;
      idxs[i] = idx == null || idx < 0 ? 0xffffffff : idx >>> 0;
      rts[i] = it.rt;
      wins[i] = it.window;
    }
    const binPtr = alloc(bin.length);
    heapWrite(binPtr, bin);
    const idxsU8 = new Uint8Array(
      idxs.buffer,
      idxs.byteOffset,
      idxs.byteLength
    );
    const rtsU8 = new Uint8Array(rts.buffer, rts.byteOffset, rts.byteLength);
    const winU8 = new Uint8Array(wins.buffer, wins.byteOffset, wins.byteLength);
    const idxsPtr = alloc(idxsU8.length);
    const rtsPtr = alloc(rtsU8.length);
    const winPtr = alloc(winU8.length);
    heapWrite(idxsPtr, idxsU8);
    heapWrite(rtsPtr, rtsU8);
    heapWrite(winPtr, winU8);
    let pOpts = 0;
    if (options) {
      writeOptions(SCRATCH_OPTS, options);
      pOpts = SCRATCH_OPTS;
    }
    const rc = get_peaks_from_chrom(
      binPtr,
      bin.length,
      idxsPtr,
      rtsPtr,
      winPtr,
      n,
      pOpts,
      1, // hard coded to 1 (No multicore in browser)
      SCRATCH_JSON
    );
    free(binPtr, bin.length);
    free(idxsPtr, idxsU8.length);
    free(rtsPtr, rtsU8.length);
    free(winPtr, winU8.length);
    if (rc !== 0) throw new Error(`C_get_peaks_from_chrom failed: ${rc}`);
    const { ptr, len } = readBuf(SCRATCH_JSON);
    const bytes = heapSlice(ptr, len);
    free(ptr, len);
    return JSON.parse(td.decode(bytes));
  };

  const getPeaksFromEic = (
    bin: Uint8Array,
    items: EicItem[],
    window: { from: number; to: number },
    options?: FindPeaksOptions
  ) => {
    const n = items.length;
    if (n === 0) return [];
    const rts = new Float64Array(n);
    const mzs = new Float64Array(n);
    const ranges = new Float64Array(n);
    const idStrings = new Array<string>(n);
    for (let i = 0; i < n; i++) {
      const it = items[i] ?? ({} as any);
      rts[i] = Number(it.rt) || 0;
      mzs[i] = Number(it.mz) || 0;
      const rng =
        (typeof it.range === "number" && it.range > 0 ? it.range : undefined) ??
        (typeof it.ranges === "number" && it.ranges > 0
          ? it.ranges
          : undefined) ??
        (typeof it.sd === "number" && it.sd > 0 ? it.sd : undefined) ??
        (typeof it.window === "number" && it.window > 0
          ? it.window
          : undefined) ??
        0.25;
      ranges[i] = rng;
      idStrings[i] = typeof it.id === "string" ? it.id : "";
    }
    const encodedIds = idStrings.map((s) => te.encode(s));
    let idsBufLen = 0;
    for (const b of encodedIds) idsBufLen += b.length;
    const idsBufU8 = new Uint8Array(idsBufLen);
    const idsOff = new Uint32Array(n);
    const idsLen = new Uint32Array(n);
    {
      let cur = 0;
      for (let i = 0; i < n; i++) {
        const b = encodedIds[i];
        idsOff[i] = cur >>> 0;
        idsLen[i] = b.length >>> 0;
        idsBufU8.set(b, cur);
        cur += b.length;
      }
    }
    const binPtr = alloc(bin.length);
    heapWrite(binPtr, bin);
    const rtsU8 = new Uint8Array(rts.buffer, rts.byteOffset, rts.byteLength);
    const mzsU8 = new Uint8Array(mzs.buffer, mzs.byteOffset, mzs.byteLength);
    const rngU8 = new Uint8Array(
      ranges.buffer,
      ranges.byteOffset,
      ranges.byteLength
    );
    const offU8 = new Uint8Array(
      idsOff.buffer,
      idsOff.byteOffset,
      idsOff.byteLength
    );
    const lenU8 = new Uint8Array(
      idsLen.buffer,
      idsLen.byteOffset,
      idsLen.byteLength
    );
    const rtsPtr = alloc(rtsU8.length);
    const mzsPtr = alloc(mzsU8.length);
    const rngPtr = alloc(rngU8.length);
    const offPtr = alloc(offU8.length);
    const lenPtr = alloc(lenU8.length);
    const idbPtr = alloc(idsBufU8.length);
    heapWrite(rtsPtr, rtsU8);
    heapWrite(mzsPtr, mzsU8);
    heapWrite(rngPtr, rngU8);
    heapWrite(offPtr, offU8);
    heapWrite(lenPtr, lenU8);
    heapWrite(idbPtr, idsBufU8);
    let pOpts = 0;
    if (options && Object.keys(options).length) {
      writeOptions(SCRATCH_OPTS, options);
      pOpts = SCRATCH_OPTS;
    }
    const rc = get_peaks_from_eic(
      binPtr,
      bin.length,
      rtsPtr,
      mzsPtr,
      rngPtr,
      offPtr,
      lenPtr,
      idbPtr,
      idsBufU8.length,
      n,
      window.from,
      window.to,
      pOpts,
      1, // hard coded to 1 (No multicore in browser)
      SCRATCH_JSON
    );
    free(binPtr, bin.length);
    free(rtsPtr, rtsU8.length);
    free(mzsPtr, mzsU8.length);
    free(rngPtr, rngU8.length);
    free(offPtr, offU8.length);
    free(lenPtr, lenU8.length);
    free(idbPtr, idsBufU8.length);
    if (rc !== 0) throw new Error(`get_peaks_from_eic failed: ${rc}`);
    const { ptr, len } = readBuf(SCRATCH_JSON);
    const bytes = heapSlice(ptr, len);
    free(ptr, len);
    const out = JSON.parse(td.decode(bytes));
    return out;
  };

  const findPeaks = (
    x: Float64Array,
    y: Float64Array,
    options: FindPeaksOptions = {}
  ) => {
    if (!(x instanceof Float64Array))
      throw new Error("findPeaks: x must be Float64Array");
    if (!(y instanceof Float64Array))
      throw new Error("findPeaks: y must be Float64Array");
    if (x.length !== y.length)
      throw new Error("findPeaks: x,y length mismatch");
    const xb = new Uint8Array(x.buffer, x.byteOffset, x.byteLength);
    const yb = new Uint8Array(y.buffer, y.byteOffset, y.byteLength);
    const xp = alloc(xb.length);
    const yp = alloc(yb.length);
    heapWrite(xp, xb);
    heapWrite(yp, yb);
    let pOpts = 0;
    if (options && Object.keys(options).length > 0) {
      writeOptions(SCRATCH_OPTS, options);
      pOpts = SCRATCH_OPTS;
    }
    const rc = find_peaks(xp, yp, x.length, pOpts, SCRATCH_JSON);
    free(xp, xb.length);
    free(yp, yb.length);
    if (rc !== 0) throw new Error(`find_peaks failed: ${rc}`);
    const { ptr, len } = readBuf(SCRATCH_JSON);
    const bytes = heapSlice(ptr, len);
    free(ptr, len);
    return JSON.parse(td.decode(bytes));
  };

  const findNoiseLevel = (y: Float32Array | ArrayLike<number>) => {
    const y32 =
      y instanceof Float32Array ? y : new Float32Array(y as ArrayLike<number>);
    if (y32.length === 0) return Infinity;
    const bytes = new Uint8Array(y32.buffer, y32.byteOffset, y32.byteLength);
    const p = alloc(bytes.length);
    heapWrite(p, bytes);
    const noise = find_noise_level(p, y32.length);
    free(p, bytes.length);
    return noise;
  };

  const calculateBaseline = (y: Float64Array, options?: BaselineOptions) => {
    if (!(y instanceof Float64Array))
      throw new Error("calculateBaseline: y must be Float64Array");
    const yU8 = new Uint8Array(y.buffer, y.byteOffset, y.byteLength);
    const yPtr = alloc(yU8.length);
    heapWrite(yPtr, yU8);
    const win =
      typeof options?.baselineWindow === "number"
        ? options.baselineWindow | 0
        : 0;
    const fac =
      typeof options?.baselineWindowFactor === "number"
        ? options.baselineWindowFactor | 0
        : 0;
    const rc = calculate_baseline(yPtr, y.length, win, fac, SCRATCH_BASE);
    free(yPtr, yU8.length);
    if (rc !== 0) throw new Error(`calculate_baseline failed: ${rc}`);
    const { ptr, len } = readBuf(SCRATCH_BASE);
    const bytes = heapSlice(ptr, len);
    free(ptr, len);
    return new Float64Array(bytes.buffer.slice(0));
  };

  const readAndFreeF64 = (bufPtr: number): Float64Array => {
    const { ptr, len } = readBuf(bufPtr);
    if (!ptr || !len) return new Float64Array(0);
    const bytes = heapSlice(ptr, len);
    free(ptr, len);
    return new Float64Array(bytes.buffer.slice(0));
  };
  const readAndFreeU32 = (bufPtr: number): Uint32Array => {
    const { ptr, len } = readBuf(bufPtr);
    if (!ptr || !len) return new Uint32Array(0);
    const bytes = heapSlice(ptr, len);
    free(ptr, len);
    return new Uint32Array(bytes.buffer.slice(0));
  };

  const collectMs1Scans = (
    bin: Uint8Array,
    fromTo: { from: number; to: number }
  ) => {
    const u8 = bin;
    if (!(u8 instanceof Uint8Array) || !u8.byteLength)
      throw new Error(
        "collectMs1Scans: first arg must be BIN1 bytes (Uint8Array)"
      );
    const from = +fromTo.from,
      to = +fromTo.to;
    if (!Number.isFinite(from) || !Number.isFinite(to) || to < from)
      throw new Error(`collectMs1Scans: invalid window: ${from}..${to}`);

    const binPtr = alloc(u8.length);
    heapWrite(binPtr, u8);

    const bRt = alloc(BUF_PAIR_BYTES);
    const bOff = alloc(BUF_PAIR_BYTES);
    const bLen = alloc(BUF_PAIR_BYTES);
    const bMz = alloc(BUF_PAIR_BYTES);
    const bInt = alloc(BUF_PAIR_BYTES);

    const rc = collect_ms1_scans(
      binPtr,
      u8.length,
      from,
      to,
      bRt,
      bOff,
      bLen,
      bMz,
      bInt
    );

    free(binPtr, u8.length);

    if (rc !== 0) {
      // clean up any partial outputs
      try {
        readAndFreeF64(bRt);
      } catch {}
      try {
        readAndFreeU32(bOff);
      } catch {}
      try {
        readAndFreeU32(bLen);
      } catch {}
      try {
        readAndFreeF64(bMz);
      } catch {}
      try {
        readAndFreeF64(bInt);
      } catch {}
      free(bRt, BUF_PAIR_BYTES);
      free(bOff, BUF_PAIR_BYTES);
      free(bLen, BUF_PAIR_BYTES);
      free(bMz, BUF_PAIR_BYTES);
      free(bInt, BUF_PAIR_BYTES);
      throw new Error(`collect_ms1_scans failed: ${rc}`);
    }

    const rt = readAndFreeF64(bRt);
    const offsets = readAndFreeU32(bOff);
    const lengths = readAndFreeU32(bLen);
    const mz = readAndFreeF64(bMz);
    const intensity = readAndFreeF64(bInt);

    free(bRt, BUF_PAIR_BYTES);
    free(bOff, BUF_PAIR_BYTES);
    free(bLen, BUF_PAIR_BYTES);
    free(bMz, BUF_PAIR_BYTES);
    free(bInt, BUF_PAIR_BYTES);

    return { rt, offsets, lengths, mz, intensity };
  };

  const findFeature = (
    bin: Uint8Array,
    targets: Target[],
    options: FindFeaturesOptions = {}
  ): FoundFeature[] => {
    const n = targets.length;
    if (n === 0) return [];

    const rts = new Float64Array(n);
    const mzs = new Float64Array(n);
    const ranges = new Float64Array(n);
    const ids = new Array<string>(n);
    for (let i = 0; i < n; i++) {
      const t = targets[i]!;
      rts[i] = +t.rt;
      mzs[i] = +t.mz;
      ranges[i] = t.range ?? 0.5;
      ids[i] = t.id ?? "";
    }

    const scanPpm = options.scanEic?.ppmTolerance ?? 10;
    const scanMz = options.scanEic?.mzTolerance ?? 0.003;
    const eicPpm = options.eic?.ppmTolerance ?? 20;
    const eicMz = options.eic?.mzTolerance ?? 0.005;

    const enc = ids.map((s) => te.encode(s));
    const idsBufLen = enc.reduce((a, b) => a + b.length, 0);
    const idsBufU8 = new Uint8Array(idsBufLen);
    const idsOff = new Uint32Array(n);
    const idsLen = new Uint32Array(n);

    let cur = 0;
    for (let i = 0; i < n; i++) {
      const b = enc[i];
      idsOff[i] = cur >>> 0;
      idsLen[i] = b.length >>> 0;
      idsBufU8.set(b, cur);
      cur += b.length;
    }

    const binPtr = alloc(bin.length);
    heapWrite(binPtr, bin);
    const rtsU8 = new Uint8Array(rts.buffer, rts.byteOffset, rts.byteLength);
    const mzsU8 = new Uint8Array(mzs.buffer, mzs.byteOffset, mzs.byteLength);
    const rngU8 = new Uint8Array(
      ranges.buffer,
      ranges.byteOffset,
      ranges.byteLength
    );
    const offU8 = new Uint8Array(
      idsOff.buffer,
      idsOff.byteOffset,
      idsOff.byteLength
    );
    const lenU8 = new Uint8Array(
      idsLen.buffer,
      idsLen.byteOffset,
      idsLen.byteLength
    );
    const rtsPtr = alloc(rtsU8.length),
      mzsPtr = alloc(mzsU8.length),
      rngPtr = alloc(rngU8.length);
    const offPtr = alloc(offU8.length),
      lenPtr = alloc(lenU8.length),
      idbPtr = alloc(idsBufU8.length);
    heapWrite(rtsPtr, rtsU8);
    heapWrite(mzsPtr, mzsU8);
    heapWrite(rngPtr, rngU8);
    heapWrite(offPtr, offU8);
    heapWrite(lenPtr, lenU8);
    heapWrite(idbPtr, idsBufU8);

    let pOpts = 0;
    if (options.findPeak) {
      writeOptions(SCRATCH_OPTS, options.findPeak);
      pOpts = SCRATCH_OPTS;
    }

    const rc = find_feature(
      binPtr,
      bin.length,
      rtsPtr,
      mzsPtr,
      rngPtr,
      offPtr,
      lenPtr,
      idbPtr,
      idsBufU8.length,
      n,
      scanPpm,
      scanMz,
      eicPpm,
      eicMz,
      pOpts,
      1, // cores = 1
      SCRATCH_JSON
    );

    free(binPtr, bin.length);
    free(rtsPtr, rtsU8.length);
    free(mzsPtr, mzsU8.length);
    free(rngPtr, rngU8.length);
    free(offPtr, offU8.length);
    free(lenPtr, lenU8.length);
    free(idbPtr, idsBufU8.length);
    if (rc !== 0) throw new Error(`find_feature failed: ${rc}`);

    const { ptr, len } = readBuf(SCRATCH_JSON);
    const out = JSON.parse(td.decode(heapSlice(ptr, len)));
    free(ptr, len);
    return out as FoundFeature[];
  };

  return {
    parseMzML,
    calculateEic,
    getPeak,
    getPeaksFromChrom,
    getPeaksFromEic,
    findPeaks,
    findNoiseLevel,
    calculateBaseline,
    collectMs1Scans,
    findFeature,
    __debug: { memory, exports: ex, heapBytes: () => memory.buffer.byteLength },
  };
}
