import type { Backend } from "./backend";
import { SampleFile } from "./sampleFile";
import { camelizeKeys, toCores, toUint8 } from "./shared";
import { encodeTargetIds, unpackTargets } from "./pack";
import type {
  BinaryInput,
  PeakOptions,
  BaselineOptions,
  Peak,
  Target,
  ChromItem,
  ChromPeak,
  Feature,
  FoundFeature,
  FindFeaturesOptions,
  GetFeaturesOptions,
  ConsensusFeature,
  FromTo,
  CentroidScan,
  SpectrumSummary,
  ScanQuery,
} from "../types/types";

export type {
  BinaryInput,
  PeakOptions,
  BaselineOptions,
  Peak,
  Target,
  ChromItem,
  ChromPeak,
  Feature,
  FoundFeature,
  FindFeaturesOptions,
  GetFeaturesOptions,
  ConsensusFeature,
  FromTo,
  CentroidScan,
  SpectrumSummary,
  ScanQuery,
};

const DEFAULTS = {
  eic: { mzTolerance: 0.005, ppmTolerance: 20.0 },
  grid: { start: 40, end: 1000, stepSize: 0.005 },
  grouping: {
    ppmTolerance: 10.0,
    mzTolerance: 0.0025,
    rtTolerance: 0.05,
    frequency: 1,
  },
  findPeak: {
    minIntensity: 150,
    minPeakWidthPoints: 5,
    autoNoise: true,
    autoBaseline: true,
    minSnr: 1,
  },
} as const;

let _backend: Backend | null = null;
let _initPromise: Promise<void> | null = null;

export function setBackend(b: Backend): void {
  _backend = b;
}

export function setInitPromise(p: Promise<void>): void {
  _initPromise = p;
}

function backend(): Backend {
  if (!_backend || !_backend.ready) {
    throw new Error(
      "Backend not initialized. Call init() first (WASM) or ensure the native module loaded (Node).",
    );
  }
  return _backend;
}

async function backendAsync(): Promise<Backend> {
  if (_initPromise) await _initPromise;
  return backend();
}

function assertFile(file: SampleFile, caller: string): void {
  if (!(file instanceof SampleFile) || !file._handle) {
    throw new Error(`${caller}: expects a valid SampleFile object`);
  }
}

/**
 * Parse an mzML file buffer into a {@link SampleFile}.
 *
 * @param data - Raw mzML file bytes.
 * @returns Loaded sample file.
 */
export async function parseMzML(data: BinaryInput): Promise<SampleFile> {
  const b = await backendAsync();
  const handle = b.parseMzML(toUint8(data));
  return new SampleFile(handle, b);
}

/**
 * Parse an ion binary buffer into a {@link SampleFile}.
 *
 * @param data - Raw ion file bytes.
 * @param options.maxCacheSize - Maximum scan cache size. Default 0 (unlimited).
 * @returns Loaded sample file.
 */
export async function parseIon(
  data: BinaryInput,
  options: { maxCacheSize?: number } = {},
): Promise<SampleFile> {
  const { maxCacheSize = 0 } = options;
  if (maxCacheSize < 0 || !Number.isFinite(maxCacheSize)) {
    throw new TypeError(
      "parseIon: maxCacheSize must be a non-negative integer",
    );
  }
  const b = await backendAsync();
  const handle = b.parseBin(toUint8(data), maxCacheSize);
  return new SampleFile(handle, b);
}

/**
 * Return the sample data as a JSON string.
 *
 * @param file - Loaded sample file.
 * @returns JSON string representation of the sample.
 */
export function ionToJson(file: SampleFile): string {
  assertFile(file, "ionToJson");
  return JSON.stringify(camelizeKeys(backend().fileToJson(file._handle!)));
}

/**
 * Serialize a sample back to mzML format.
 *
 * @param file - Loaded sample file.
 * @returns Full mzML file content as a string.
 */
export function ionToMzml(file: SampleFile): string {
  assertFile(file, "ionToMzml");
  return backend().fileToMzml(file._handle!);
}

/**
 * Encode a sample as compressed ion binary bytes.
 *
 * @param file - Loaded sample file.
 * @param options.level - Compression level, 0 (none) to 22 (max). Default 5.
 * @param options.f32Compress - Compress intensity values to 32-bit float. Default false.
 * @returns Raw ion binary bytes.
 */
export function mzmlToIon(
  file: SampleFile,
  options: { level?: number; f32Compress?: boolean } = {},
): Uint8Array {
  assertFile(file, "mzmlToIon");
  const { level = 5, f32Compress = false } = options;
  if (
    typeof level !== "number" ||
    !Number.isFinite(level) ||
    (level | 0) !== level ||
    level < 0 ||
    level > 22
  ) {
    throw new RangeError("mzmlToIon: level must be an integer in [0,22]");
  }
  if (typeof f32Compress !== "boolean") {
    throw new TypeError("mzmlToIon: f32Compress must be a boolean");
  }
  return backend().fileToBin(file._handle!, level, f32Compress);
}

/**
 * Extract an extracted ion chromatogram (EIC) for a target m/z.
 *
 * @param file - Loaded sample file.
 * @param mz - Target m/z value.
 * @param fromTo - Retention time range `{ from, to }` in minutes.
 * @param ppmTol - PPM tolerance. Default 20.
 * @param mzTol - Absolute m/z tolerance in Da. Default 0.005.
 * @returns Object with `x` (retention times) and `y` (intensities) arrays.
 */
export function calculateEic(
  file: SampleFile,
  mz: number,
  fromTo: FromTo,
  ppmTol = 20,
  mzTol = 0.005,
): { x: Float64Array; y: Float64Array } {
  assertFile(file, "calculateEic");
  if (typeof mz !== "number" || !Number.isFinite(mz) || mz <= 0) {
    throw new RangeError(
      "calculateEic: targetMz must be a positive finite number",
    );
  }
  const { from, to } = fromTo;
  if (!Number.isFinite(from) || !Number.isFinite(to)) {
    throw new RangeError("calculateEic: from and to must be finite numbers");
  }
  return backend().calculateEic(file._handle!, mz, from, to, ppmTol, mzTol);
}

/**
 * Detect peaks in a chromatographic trace.
 *
 * @param x - Retention time array.
 * @param y - Intensity array.
 * @param opts - Peak detection options.
 * @returns Array of detected peaks.
 */
export function findPeaks(
  x: Float64Array,
  y: Float64Array,
  opts?: PeakOptions,
): Peak[] {
  if (!(x instanceof Float64Array) || !(y instanceof Float64Array)) {
    throw new TypeError("findPeaks: x and y must be Float64Array");
  }
  if (x.length !== y.length) {
    throw new RangeError("findPeaks: x and y must have equal length");
  }
  return backend().findPeaks(x, y, opts);
}

/**
 * Find a single peak near a target retention time.
 *
 * @param x - Retention time array.
 * @param y - Intensity array.
 * @param rt - Target retention time in minutes.
 * @param range - Search window half-width in minutes.
 * @param opts - Peak detection options.
 * @returns The best matching peak.
 */
export function getPeak(
  x: Float64Array,
  y: Float64Array,
  rt: number,
  range: number,
  opts?: PeakOptions,
): Peak {
  if (!(x instanceof Float64Array) || !(y instanceof Float64Array)) {
    throw new TypeError("getPeak: x and y must be Float64Array");
  }
  if (x.length !== y.length) {
    throw new RangeError("getPeak: x and y must have equal length");
  }
  if (!Number.isFinite(rt)) {
    throw new RangeError("getPeak: rt must be a finite number");
  }
  if (!Number.isFinite(range) || range <= 0) {
    throw new RangeError("getPeak: range must be a positive finite number");
  }
  return backend().getPeak(x, y, rt, range, opts);
}

/**
 * Estimate the noise level of an intensity array.
 *
 * @param y - Intensity array.
 * @returns Estimated noise level.
 */
export function findNoiseLevel(y: Float64Array | Float32Array): number {
  return backend().findNoiseLevel(y);
}

/**
 * Compute the rolling baseline of an intensity array.
 *
 * @param y - Intensity array.
 * @param options.lambda - airPLS smoothing parameter. Default auto.
 * @param options.maxIterations - Maximum airPLS iterations. Default auto.
 * @returns Baseline array, same length as `y`.
 */
export function calculateBaseline(
  y: Float64Array | ArrayLike<number>,
  options?: BaselineOptions,
): Float64Array {
  const y64 =
    y instanceof Float64Array ? y : new Float64Array(y as ArrayLike<number>);
  return backend().calculateBaseline(
    y64,
    (options?.lambda as any) | 0,
    (options?.maxIterations as any) | 0,
  );
}

const QUERY_RT_RANGE = 0;
const QUERY_CLOSEST_RT = 1;
const QUERY_MZ_RANGE = 2;
const QUERY_CLOSEST_MZ = 3;

function assertLevel(level: number, caller: string): void {
  if (
    typeof level !== "number" ||
    !Number.isFinite(level) ||
    (level | 0) !== level ||
    level < 0 ||
    level > 255
  ) {
    throw new RangeError(`${caller}: level must be an integer in [0,255]`);
  }
}

/**
 * Retrieve centroid scans from a sample by retention time or precursor m/z.
 *
 * @param file - Loaded sample file.
 * @param query - Query by RT range, closest RT, m/z range, or closest m/z.
 * @param level - MS level (1 = MS1, 2 = MS2, etc.). Default 1.
 * @returns Array of scans for range queries, single scan (or null) for closest queries.
 */
export function getScans(
  file: SampleFile,
  query: ScanQuery,
  level = 1,
): CentroidScan[] | CentroidScan | null {
  assertFile(file, "getScans");
  assertLevel(level, "getScans");

  let queryType: number;
  let a: number;
  let b: number;

  const assertFinite = (v: number, label: string) => {
    if (!Number.isFinite(v))
      throw new RangeError(`getScans: ${label} must be finite`);
  };

  if ("rt" in query) {
    if ("from" in query.rt) {
      assertFinite(query.rt.from, "rt.from");
      assertFinite(query.rt.to, "rt.to");
      queryType = QUERY_RT_RANGE;
      a = query.rt.from;
      b = query.rt.to;
    } else {
      assertFinite(query.rt.closest, "rt.closest");
      queryType = QUERY_CLOSEST_RT;
      a = query.rt.closest;
      b = NaN;
    }
  } else {
    if ("from" in query.mz) {
      assertFinite(query.mz.from, "mz.from");
      assertFinite(query.mz.to, "mz.to");
      queryType = QUERY_MZ_RANGE;
      a = query.mz.from;
      b = query.mz.to;
    } else {
      assertFinite(query.mz.closest, "mz.closest");
      if (query.mz.closest <= 0)
        throw new RangeError("getScans: mz.closest must be positive");
      queryType = QUERY_CLOSEST_MZ;
      a = query.mz.closest;
      b = NaN;
    }
  }

  const raw = camelizeKeys(
    backend().getScans(file._handle!, queryType, a, b, level),
  );
  if (queryType === QUERY_CLOSEST_RT || queryType === QUERY_CLOSEST_MZ) {
    return raw?.length ? raw[0] : null;
  }
  return raw as CentroidScan[];
}

/**
 * Extract and pick peaks for a list of targets using their EIC traces.
 *
 * @param file - Loaded sample file.
 * @param targets - List of targets with `id`, `rt`, `mz`, and `range`.
 * @param fromTo - Retention time range `{ from, to }` in minutes.
 * @param options - Peak detection options.
 * @param cores - Number of CPU cores to use. Default 1.
 * @returns Array of found features, one per target.
 */
export function getPeaksFromEic(
  file: SampleFile,
  targets: Target[],
  fromTo: FromTo,
  options?: PeakOptions,
  cores = 1,
): FoundFeature[] {
  assertFile(file, "getPeaksFromEic");
  const { from = 0.5, to = 5 } = fromTo;
  const { rts, mzs, ranges, ids } = unpackTargets(targets);
  const { packedBytes, offsets, lengths } = encodeTargetIds(ids);
  return backend().getPeaksFromEic(
    file._handle!,
    rts,
    mzs,
    ranges,
    offsets,
    lengths,
    packedBytes,
    packedBytes.length,
    targets.length,
    from,
    to,
    options,
    toCores(cores),
  );
}

/**
 * Pick peaks from pre-computed chromatogram traces.
 *
 * @param file - Loaded sample file.
 * @param items - List of chromatogram items with `idx`, `rt`, and `window`.
 * @param options - Peak detection options.
 * @param cores - Number of CPU cores to use. Default 1.
 * @returns Array of chromatogram peaks.
 */
export function getPeaksFromChrom(
  file: SampleFile,
  items: ChromItem[],
  options?: PeakOptions,
  cores = 1,
): ChromPeak[] {
  assertFile(file, "getPeaksFromChrom");
  const n = items.length;
  const indices = new Uint32Array(n);
  const rts = new Float64Array(n);
  const windows = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const item = items[i];
    const idx = Number.isFinite(item.idx)
      ? (item.idx as number)
      : Number.isFinite(item.index)
        ? (item.index as number)
        : -1;
    indices[i] = idx >= 0 ? idx >>> 0 : 0xffffffff;
    rts[i] = +item.rt;
    windows[i] = +(item.window ?? item.range ?? 0);
  }
  return backend().getPeaksFromChrom(
    file._handle!,
    indices,
    rts,
    windows,
    n,
    options,
    toCores(cores),
  );
}

/**
 * Detect all chromatographic features across an m/z grid.
 *
 * @param file - Loaded sample file.
 * @param fromTo - Retention time range `{ from, to }` in minutes.
 * @param options - Feature detection options (eic, grid, findPeak, cores).
 * @returns Array of detected features.
 */
export function findFeatures(
  file: SampleFile,
  fromTo: FromTo,
  options: FindFeaturesOptions = {},
): Feature[] {
  assertFile(file, "findFeatures");
  const {
    eic = DEFAULTS.eic,
    grid = DEFAULTS.grid,
    findPeak = DEFAULTS.findPeak,
    cores = 1,
    useGpu = false,
    batchSize = 0,
  } = options;
  const { from, to } = fromTo;
  const eicPpm =
    typeof eic.ppmTolerance === "number" &&
    Number.isFinite(eic.ppmTolerance) &&
    eic.ppmTolerance >= 0
      ? eic.ppmTolerance
      : NaN;
  const eicMz =
    typeof eic.mzTolerance === "number" &&
    Number.isFinite(eic.mzTolerance) &&
    eic.mzTolerance >= 0
      ? eic.mzTolerance
      : NaN;
  const gridStart =
    typeof grid.start === "number" && Number.isFinite(grid.start)
      ? grid.start
      : NaN;
  const gridEnd =
    typeof grid.end === "number" && Number.isFinite(grid.end) ? grid.end : NaN;
  const gridStep =
    typeof grid.stepSize === "number" &&
    Number.isFinite(grid.stepSize) &&
    grid.stepSize > 0
      ? grid.stepSize
      : NaN;
  return backend().findFeatures(
    file._handle!,
    from,
    to,
    eicPpm,
    eicMz,
    gridStart,
    gridEnd,
    gridStep,
    findPeak,
    toCores(cores),
    useGpu ? 1 : 0,
    batchSize | 0,
  );
}

/**
 * Find a feature for each target using its EIC and scan data.
 *
 * @param file - Loaded sample file.
 * @param targets - List of targets with `id`, `rt`, `mz`, and `ranges`.
 * @param options - Detection options including `scanEic`, `eic`, `findPeak`, and `cores`.
 * @returns Array of found features, one per target.
 */
export function findFeature(
  file: SampleFile,
  targets: Target[],
  options: FindFeaturesOptions & {
    scanEic?: { ppmTolerance?: number; mzTolerance?: number };
    cores?: number;
  } = {},
): FoundFeature[] {
  assertFile(file, "findFeature");

  const {
    scanEic: { ppmTolerance: scanPpm = 10, mzTolerance: scanMz = 0.003 } = {},
    eic: { ppmTolerance: eicPpm = 20, mzTolerance: eicMz = 0.005 } = {},
  } = options;

  const check = (raw: unknown, label: string): void => {
    if (
      raw !== undefined &&
      (typeof raw !== "number" || !Number.isFinite(raw) || raw < 0)
    ) {
      throw new TypeError(`${label} must be a finite, non-negative number`);
    }
  };
  check(scanPpm, "scanEic.ppmTolerance");
  check(scanMz, "scanEic.mzTolerance");
  check(eicPpm, "eic.ppmTolerance");
  check(eicMz, "eic.mzTolerance");

  const { rts, mzs, ranges, ids } = unpackTargets(targets, 0.5);
  const { packedBytes, offsets, lengths } = encodeTargetIds(ids);

  return backend().findFeature(
    file._handle!,
    rts,
    mzs,
    ranges,
    offsets,
    lengths,
    packedBytes,
    packedBytes.length,
    targets.length,
    toCores(options.cores),
    +scanPpm,
    +scanMz,
    +eicPpm,
    +eicMz,
    options.findPeak,
  );
}

/**
 * Detect and group consensus features across multiple samples in a directory.
 * Node.js only — not available in the WASM build.
 *
 * @param directoryPath - Path to a directory containing ion files.
 * @param fromTo - Retention time range `{ from, to }` in minutes.
 * @param options - Feature detection and grouping options.
 * @returns Array of consensus features aligned across samples.
 */
export function getFeatures(
  directoryPath: string,
  fromTo: FromTo,
  options: GetFeaturesOptions = {},
): ConsensusFeature[] {
  const back = backend();
  if (!back.getFeatures) {
    throw new Error("getFeatures is not supported in the WASM build");
  }
  if (typeof directoryPath !== "string" || directoryPath.length === 0) {
    throw new Error("getFeatures: Valid directory path is required");
  }
  const {
    eic = DEFAULTS.eic,
    grid = DEFAULTS.grid,
    grouping = DEFAULTS.grouping,
    findPeak = DEFAULTS.findPeak,
    cores = 1,
    useGpu = false,
    batchSize = 0,
  } = options;
  const { from, to } = fromTo;
  const eicPpm = Number.isFinite(eic.ppmTolerance)
    ? (eic.ppmTolerance as number)
    : NaN;
  const eicMz = Number.isFinite(eic.mzTolerance)
    ? (eic.mzTolerance as number)
    : NaN;
  const gridStart = Number.isFinite(grid.start) ? (grid.start as number) : NaN;
  const gridEnd = Number.isFinite(grid.end) ? (grid.end as number) : NaN;
  const gridStep = Number.isFinite(grid.stepSize)
    ? (grid.stepSize as number)
    : NaN;
  const groupPpm = grouping.ppmTolerance ?? DEFAULTS.grouping.ppmTolerance;
  const groupMz = grouping.mzTolerance ?? DEFAULTS.grouping.mzTolerance;
  const groupRt = grouping.rtTolerance ?? DEFAULTS.grouping.rtTolerance;
  const prevalence = grouping.frequency ?? DEFAULTS.grouping.frequency;
  return back.getFeatures!(
    directoryPath,
    +from,
    +to,
    +eicPpm,
    +eicMz,
    +gridStart,
    +gridEnd,
    +gridStep,
    +groupPpm,
    +groupMz,
    +groupRt,
    +prevalence,
    findPeak,
    toCores(cores),
    useGpu ? 1 : 0,
    batchSize | 0,
  );
}
