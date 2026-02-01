import * as fs from "fs";
import * as path from "path";

function firstExisting(...candidates: string[]) {
  for (const p of candidates) if (fs.existsSync(p)) return p;
  return candidates[0];
}

function platformLibPath(process: NodeJS.Process): string {
  const base = path.join(__dirname, "..", "native");
  const { platform, arch } = process;

  const file =
    platform === "win32"
      ? "msutils.dll"
      : platform === "darwin"
        ? "libmsutils.dylib"
        : "libmsutils.so";

  let dir: string;

  if (platform === "darwin") {
    dir =
      arch === "arm64"
        ? firstExisting(
            path.join(base, "darwin-arm64"),
            path.join(base, "macos-arm64"),
          )
        : firstExisting(
            path.join(base, "darwin-x64"),
            path.join(base, "macos-x86_64"),
          );
  } else if (platform === "linux") {
    dir =
      arch === "arm64"
        ? firstExisting(
            path.join(base, "linux-arm64-gnu"),
            path.join(base, "linux-arm64"),
          )
        : firstExisting(
            path.join(base, "linux-x64-gnu"),
            path.join(base, "linux-x86_64"),
          );
  } else if (platform === "win32") {
    dir = firstExisting(
      path.join(base, "win32-x64"),
      path.join(base, "windows-x86_64"),
    );
  } else {
    throw new Error(`Unsupported ${platform}/${arch}`);
  }

  return path.join(dir, file);
}

const addonPath = path.join(
  __dirname,
  "..",
  "build",
  "Release",
  "msutils.node",
);
const native = require(addonPath);
if (typeof native.bind === "function") native.bind(platformLibPath(process));

type BinaryInput = Uint8Array | ArrayBuffer;

function toBuffer(v: BinaryInput): Buffer {
  if (Buffer.isBuffer(v)) return v;
  if (v instanceof ArrayBuffer) return Buffer.from(v);
  return Buffer.from(v.buffer, v.byteOffset, v.byteLength);
}

function toCores(cores: unknown): number {
  const v = typeof cores === "number" && Number.isFinite(cores) ? cores | 0 : 1;
  return v > 0 ? v : 1;
}

function parseJson<T>(s: string): T {
  return JSON.parse(s) as T;
}

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

export function packPeakOptions(opts?: PeakOptions): Buffer | undefined {
  if (!opts) return undefined;

  const b = Buffer.allocUnsafe(64);

  const f64 = (v: unknown, off: number) =>
    b.writeDoubleLE(Number.isFinite(v as number) ? Number(v) : NaN, off);

  const i32 = (v: unknown, off: number) =>
    b.writeInt32LE(
      typeof v === "boolean"
        ? v
          ? 1
          : 0
        : Number.isFinite(v as number)
          ? (v as number) | 0
          : 0,
      off,
    );

  f64(opts.integralThreshold, 0);
  f64(opts.intensityThreshold, 8);
  i32(opts.widthThreshold, 16);
  b.writeInt32LE(0, 20);
  f64(opts.noise, 24);
  i32(opts.autoNoise, 32);
  i32(opts.autoBaseline, 36);
  i32(opts.baselineWindow, 40);
  i32(opts.baselineWindowFactor, 44);
  i32(opts.allowOverlap, 48);
  i32(opts.windowSize, 52);
  f64(opts.snRatio, 56);

  return b;
}

const parseMzMLFn = native.parseMzml || native.parseMzML;

export function parseMzML(data: BinaryInput): Buffer {
  return parseMzMLFn(toBuffer(data)) as Buffer;
}

export function binToJson(bin: BinaryInput): string {
  return native.binToJson(toBuffer(bin)) as string;
}

export function convertBinToMzml(bin: BinaryInput): string {
  return native.binToMzML(toBuffer(bin)) as string;
}

export function calculateEic(
  bin: BinaryInput,
  targets: number,
  from: number,
  to: number,
  ppmTol = 20,
  mzTol = 0.005,
) {
  const b = toBuffer(bin);
  return native.calculateEic(b, +targets, from, to, ppmTol, mzTol) as {
    x: Float64Array;
    y: Float64Array;
  };
}

export function findPeaks(
  x: Float64Array,
  y: Float64Array,
  opts?: PeakOptions,
) {
  return parseJson<Peak[]>(
    native.findPeaks(x, y, packPeakOptions(opts)) as string,
  );
}

export function getPeak(
  x: Float64Array,
  y: Float64Array,
  rt: number,
  range: number,
  opts?: PeakOptions,
) {
  return parseJson<Peak>(
    native.getPeak(x, y, rt, range, packPeakOptions(opts)) as string,
  );
}

export const findNoiseLevel = native.findNoiseLevel as (
  y: Float64Array,
) => number;

export function getPeaksFromEic(
  bin: BinaryInput,
  targets: Target[],
  fromLeft = 0.5,
  toRight = 0.5,
  options?: PeakOptions,
  cores = 1,
) {
  const n = targets.length;
  const rts = new Float64Array(n);
  const mzs = new Float64Array(n);
  const rng = new Float64Array(n);
  let ids: string[] | undefined;

  for (let i = 0; i < n; i++) {
    const t = targets[i];
    rts[i] = +t.rt;
    mzs[i] = +t.mz;
    rng[i] = +t.range;

    const id = t.id;
    if (id && id.length) {
      if (!ids) ids = new Array(n);
      ids[i] = id;
    }
  }

  const json = native.getPeaksFromEic(
    toBuffer(bin),
    rts,
    mzs,
    rng,
    ids,
    +fromLeft,
    +toRight,
    packPeakOptions(options),
    toCores(cores),
  ) as string;

  return parseJson<
    Array<{
      id?: string;
      mz: number;
      ort: number;
      rt: number;
      from: number;
      to: number;
      intensity: number;
      integral: number;
      noise: number;
    }>
  >(json);
}

export function getPeaksFromChrom(
  bin: BinaryInput,
  items: ChromItem[],
  options?: PeakOptions,
  cores = 1,
) {
  const n = items.length;
  const idxs = new Uint32Array(n);
  const rts = new Float64Array(n);
  const rng = new Float64Array(n);

  for (let i = 0; i < n; i++) {
    const it = items[i];
    const idx = Number.isFinite(it.idx)
      ? (it.idx as number)
      : Number.isFinite(it.index)
        ? (it.index as number)
        : -1;

    idxs[i] = idx >= 0 ? idx >>> 0 : 0xffffffff;
    rts[i] = +it.rt;
    rng[i] = +(it.window ?? it.range ?? 0);
  }

  const json = native.getPeaksFromChrom(
    toBuffer(bin),
    idxs,
    rts,
    rng,
    packPeakOptions(options),
    toCores(cores),
  ) as string;

  return parseJson<ChromPeak[]>(json);
}

export function calculateBaseline(
  y: Float64Array | ArrayLike<number>,
  options?: BaselineOptions,
): Float64Array {
  const y64 =
    y instanceof Float64Array ? y : new Float64Array(y as ArrayLike<number>);
  const win = (options?.baselineWindow as any) | 0;
  const fac = (options?.baselineWindowFactor as any) | 0;
  return native.calculateBaseline(y64, win, fac) as Float64Array;
}

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

export function findFeatures(
  data: BinaryInput,
  fromTo: { from: number; to: number },
  options: FindFeaturesOptions = {},
): Feature[] {
  const {
    eic = { mzTolerance: 0.0025, ppmTolerance: 5.0 },
    grid = { start: 20, end: 700, stepSize: 0.005 },
    findPeak = {},
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

  const s = native.findFeatures(
    toBuffer(data),
    from,
    to,
    eicPpm,
    eicMz,
    gridStart,
    gridEnd,
    gridStep,
    packPeakOptions(findPeak) ?? null,
    toCores(cores),
  ) as string;

  return parseJson<Feature[]>(s);
}

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

export function findFeature(
  data: BinaryInput,
  targets: Target[],
  options: FindFeaturesOptions & {
    scanEic?: { ppmTolerance?: number; mzTolerance?: number };
    cores?: number;
  } = {},
): FoundFeature[] {
  let {
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

  const peakBuf = packPeakOptions(options.findPeak);

  const n = targets.length;
  const rts = new Float64Array(n);
  const mzs = new Float64Array(n);
  const wins = new Float64Array(n);
  let ids: string[] | undefined;

  for (let i = 0; i < n; i++) {
    const t = targets[i];
    rts[i] = +t.rt;
    mzs[i] = +t.mz;
    wins[i] =
      typeof t.range === "number" && Number.isFinite(t.range) ? t.range : 0.5;

    const id = t.id;
    if (id && id.length) {
      if (!ids) ids = new Array(n);
      ids[i] = id;
    }
  }

  const s = native.findFeature(
    toBuffer(data),
    rts,
    mzs,
    wins,
    ids,
    +scanPpm,
    +scanMz,
    +eicPpm,
    +eicMz,
    peakBuf ?? null,
    toCores(options.cores),
  ) as string;

  return parseJson<FoundFeature[]>(s);
}

export function convertMzmlToBin(
  data: BinaryInput,
  options: { level?: number; f32Compress?: boolean } = {},
): Buffer {
  const { level = 5, f32Compress = false } = options;

  if (
    typeof level !== "number" ||
    !Number.isFinite(level) ||
    (level | 0) !== level ||
    level < 0 ||
    level > 22
  ) {
    throw new RangeError(
      "convertMzmlToBin: level must be an integer in [0,22]",
    );
  }

  if (typeof f32Compress !== "boolean") {
    throw new TypeError("convertMzmlToBin: f32Compress must be a boolean");
  }

  return native.convertMzmlToBin(
    toBuffer(data),
    level,
    f32Compress ? 1 : 0,
  ) as Buffer;
}

export function parseBin(bin: BinaryInput): Buffer {
  return native.parseBin(toBuffer(bin)) as Buffer;
}

module.exports = {
  parseMzML,
  binToJson,
  convertBinToMzml,
  calculateEic,
  getPeak,
  findPeaks,
  findNoiseLevel,
  getPeaksFromEic,
  getPeaksFromChrom,
  calculateBaseline,
  findFeatures,
  findFeature,
  convertMzmlToBin,
  parseBin,
};
