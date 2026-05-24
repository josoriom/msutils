import type { Target } from "../types/types";

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
