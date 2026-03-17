import type { PeakOptions, Target } from "../types/types";

const PEAK_OPTIONS_STRUCT_BYTES = 64;

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

export function encodeTargetIds(ids: string[]): {
  packedBytes: Uint8Array;
  offsets: Uint32Array;
  lengths: Uint32Array;
} {
  const encoder = new TextEncoder();
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

export function unpackTargets(
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
