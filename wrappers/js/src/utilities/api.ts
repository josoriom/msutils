import type { Backend } from "./backend";
import { MzMlFile } from "./mzmlFile";
import { camelizeKeys, toCores, toUint8 } from "./shared";
import { packPeakOptions, encodeTargetIds, unpackTargets } from "./pack";
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
};

const DEFAULTS = {
  eic: { mzTolerance: 0.005, ppmTolerance: 10.0 },
  grid: { start: 40, end: 1000, stepSize: 0.005 },
  grouping: {
    ppmTolerance: 5.0,
    mzTolerance: 0.0025,
    rtTolerance: 0.05,
    frequency: 1,
  },
  findPeak: {
    intensityThreshold: 150,
    widthThreshold: 5,
    autoNoise: true,
    autoBaseline: true,
    snRatio: 1,
  },
} as const;

let _backend: Backend | null = null;

export function setBackend(b: Backend): void {
  _backend = b;
}

function backend(): Backend {
  if (!_backend || !_backend.ready) {
    throw new Error(
      "Backend not initialized. Call init() first (WASM) or ensure the native module loaded (Node).",
    );
  }
  return _backend;
}

function assertFile(file: MzMlFile, caller: string): void {
  if (!(file instanceof MzMlFile) || !file._handle) {
    throw new Error(`${caller}: expects a valid MzMlFile object`);
  }
}

export function parseMzML(data: BinaryInput): MzMlFile {
  const handle = backend().parseMzML(toUint8(data));
  return new MzMlFile(handle, backend());
}

export function parseBin(data: BinaryInput): MzMlFile {
  const handle = backend().parseBin(toUint8(data));
  return new MzMlFile(handle, backend());
}

export function binToJson(file: MzMlFile): string {
  assertFile(file, "binToJson");
  return JSON.stringify(camelizeKeys(backend().fileToJson(file._handle!)));
}

export function convertBinToMzml(file: MzMlFile): string {
  assertFile(file, "convertBinToMzml");
  return backend().fileToMzml(file._handle!);
}

export function mzmlToBin(
  file: MzMlFile,
  options: { level?: number; f32Compress?: boolean } = {},
): Uint8Array {
  assertFile(file, "mzmlToBin");
  const { level = 5, f32Compress = false } = options;
  if (
    typeof level !== "number" ||
    !Number.isFinite(level) ||
    (level | 0) !== level ||
    level < 0 ||
    level > 22
  ) {
    throw new RangeError("mzmlToBin: level must be an integer in [0,22]");
  }
  if (typeof f32Compress !== "boolean") {
    throw new TypeError("mzmlToBin: f32Compress must be a boolean");
  }
  return backend().fileToBin(file._handle!, level, f32Compress);
}

export function calculateEic(
  file: MzMlFile,
  targetMz: number,
  fromTo: FromTo,
  ppmTol = 20,
  mzTol = 0.005,
): { x: Float64Array; y: Float64Array } {
  assertFile(file, "calculateEic");
  const { from, to } = fromTo;
  return backend().calculateEic(
    file._handle!,
    +targetMz,
    from,
    to,
    ppmTol,
    mzTol,
  );
}

export function findPeaks(
  x: Float64Array,
  y: Float64Array,
  opts?: PeakOptions,
): Peak[] {
  return backend().findPeaks(x, y, packPeakOptions(opts));
}

export function getPeak(
  x: Float64Array,
  y: Float64Array,
  rt: number,
  range: number,
  opts?: PeakOptions,
): Peak {
  return backend().getPeak(x, y, rt, range, packPeakOptions(opts));
}

export function findNoiseLevel(y: Float64Array | Float32Array): number {
  return backend().findNoiseLevel(y);
}

export function calculateBaseline(
  y: Float64Array | ArrayLike<number>,
  options?: BaselineOptions,
): Float64Array {
  const y64 =
    y instanceof Float64Array ? y : new Float64Array(y as ArrayLike<number>);
  return backend().calculateBaseline(
    y64,
    (options?.baselineWindow as any) | 0,
    (options?.baselineWindowFactor as any) | 0,
  );
}

export function collectScans(
  file: MzMlFile,
  fromTo: FromTo,
  level = 1,
  includeMetadata = false,
): CentroidScan[] {
  assertFile(file, "collectScans");
  if (
    typeof level !== "number" ||
    !Number.isFinite(level) ||
    (level | 0) !== level ||
    level < 0 ||
    level > 255
  ) {
    throw new RangeError("collectScans: level must be an integer in [0,255]");
  }
  const { from, to } = fromTo;
  return camelizeKeys(
    backend().collectScans(file._handle!, from, to, level, includeMetadata),
  );
}

export function getPeaksFromEic(
  file: MzMlFile,
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
    packPeakOptions(options),
    toCores(cores),
  );
}

export function getPeaksFromChrom(
  file: MzMlFile,
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
    packPeakOptions(options),
    toCores(cores),
  );
}

export function findFeatures(
  file: MzMlFile,
  fromTo: FromTo,
  options: FindFeaturesOptions = {},
): Feature[] {
  assertFile(file, "findFeatures");
  const {
    eic = DEFAULTS.eic,
    grid = DEFAULTS.grid,
    findPeak = DEFAULTS.findPeak,
    cores = 1,
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
    packPeakOptions(findPeak) ?? undefined,
    toCores(cores),
  );
}

export function findFeature(
  file: MzMlFile,
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
    packPeakOptions(options.findPeak),
  );
}

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
  return back.getFeatures(
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
    packPeakOptions(findPeak) ?? undefined,
    toCores(cores),
  );
}
