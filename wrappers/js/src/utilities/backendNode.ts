/// <reference types="node" />

import * as fs from "fs";
import * as path from "path";
import type { Backend, FileHandle } from "./backend";

function firstExisting(...candidates: string[]): string {
  for (const p of candidates) if (fs.existsSync(p)) return p;
  return candidates[0];
}

function platformLibPath(proc: NodeJS.Process): string {
  const base = path.join(__dirname, "..", "..", "native");
  const { platform, arch } = proc;

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

function toNodeBuffer(data: Uint8Array): Buffer {
  if (Buffer.isBuffer(data)) return data;
  return Buffer.from(data.buffer, data.byteOffset, data.byteLength);
}

function unpackIds(
  count: number,
  offsets: Uint32Array | null,
  lengths: Uint32Array | null,
  idBytes: Uint8Array | null,
): string[] | undefined {
  if (!offsets || !lengths || !idBytes) return undefined;
  const decoder = new TextDecoder();
  const ids: string[] = new Array(count);
  for (let i = 0; i < count; i++) {
    const off = offsets[i];
    const len = lengths[i];
    ids[i] = len > 0 ? decoder.decode(idBytes.subarray(off, off + len)) : "";
  }
  return ids;
}

export class NodeBackend implements Backend {
  private native: any;
  ready = false;

  constructor(proc: NodeJS.Process) {
    const addonPath = path.join(
      __dirname,
      "..",
      "..",
      "build",
      "Release",
      "msutils.node",
    );
    this.native = require(addonPath);
    if (typeof this.native.bind === "function") {
      this.native.bind(platformLibPath(proc));
    }
    this.ready = true;
  }

  async init(): Promise<void> {}

  parseMzML(data: Uint8Array): FileHandle {
    return this.native.parseMzML(toNodeBuffer(data));
  }

  parseBin(data: Uint8Array, maxCacheSize = 0): FileHandle {
    return this.native.parseBin(toNodeBuffer(data), maxCacheSize);
  }

  freeFile(handle: FileHandle): void {
    this.native.dispose(handle);
  }

  fileToJson(handle: FileHandle): any {
    return JSON.parse(this.native.binToJson(handle) as string);
  }

  fileToMzml(handle: FileHandle): string {
    return this.native.binToMzML(handle) as string;
  }

  fileToBin(
    handle: FileHandle,
    level: number,
    f32Compress: boolean,
  ): Uint8Array {
    return this.native.mzmlToBin(handle, level, f32Compress ? 1 : 0) as Buffer;
  }

  calculateEic(
    handle: FileHandle,
    targetMz: number,
    from: number,
    to: number,
    ppmTol: number,
    mzTol: number,
  ): { x: Float64Array; y: Float64Array } {
    return this.native.calculateEic(
      handle,
      targetMz,
      from,
      to,
      ppmTol,
      mzTol,
    ) as { x: Float64Array; y: Float64Array };
  }

  getScans(
    handle: FileHandle,
    queryType: number,
    a: number,
    b: number,
    level: number,
  ): any {
    return JSON.parse(
      this.native.getScans(handle, queryType, a, b, level) as string,
    );
  }

  findPeaks(
    x: Float64Array,
    y: Float64Array,
    packedOpts: Uint8Array | undefined,
  ): any {
    return JSON.parse(
      this.native.findPeaks(
        x,
        y,
        packedOpts ? toNodeBuffer(packedOpts) : undefined,
      ) as string,
    );
  }

  getPeak(
    x: Float64Array,
    y: Float64Array,
    rt: number,
    range: number,
    packedOpts: Uint8Array | undefined,
  ): any {
    return JSON.parse(
      this.native.getPeak(
        x,
        y,
        rt,
        range,
        packedOpts ? toNodeBuffer(packedOpts) : undefined,
      ) as string,
    );
  }

  findNoiseLevel(y: Float64Array | Float32Array): number {
    const y64 = y instanceof Float64Array ? y : new Float64Array(y);
    return this.native.findNoiseLevel(y64) as number;
  }

  calculateBaseline(
    y: Float64Array,
    baselineWindow: number,
    baselineWindowFactor: number,
  ): Float64Array {
    return this.native.calculateBaseline(
      y,
      baselineWindow,
      baselineWindowFactor,
    ) as Float64Array;
  }

  getPeaksFromEic(
    handle: FileHandle,
    rts: Float64Array,
    mzs: Float64Array,
    ranges: Float64Array,
    offsets: Uint32Array | null,
    lengths: Uint32Array | null,
    idBytes: Uint8Array | null,
    _idBytesLen: number,
    count: number,
    from: number,
    to: number,
    packedOpts: Uint8Array | undefined,
    cores: number,
  ): any {
    const ids = unpackIds(count, offsets, lengths, idBytes);
    return JSON.parse(
      this.native.getPeaksFromEic(
        handle,
        rts,
        mzs,
        ranges,
        ids,
        from,
        to,
        packedOpts ? toNodeBuffer(packedOpts) : undefined,
        cores,
      ) as string,
    );
  }

  getPeaksFromChrom(
    handle: FileHandle,
    indices: Uint32Array,
    rts: Float64Array,
    windows: Float64Array,
    _count: number,
    packedOpts: Uint8Array | undefined,
    cores: number,
  ): any {
    return JSON.parse(
      this.native.getPeaksFromChrom(
        handle,
        indices,
        rts,
        windows,
        packedOpts ? toNodeBuffer(packedOpts) : undefined,
        cores,
      ) as string,
    );
  }

  findFeatures(
    handle: FileHandle,
    from: number,
    to: number,
    eicPpm: number,
    eicMz: number,
    gridStart: number,
    gridEnd: number,
    gridStep: number,
    packedOpts: Uint8Array | undefined,
    cores: number,
    useGpu: number,
    batchSize: number,
  ): any {
    return JSON.parse(
      this.native.findFeatures(
        handle,
        from,
        to,
        eicPpm,
        eicMz,
        gridStart,
        gridEnd,
        gridStep,
        packedOpts ? toNodeBuffer(packedOpts) : null,
        cores,
        useGpu,
        batchSize,
      ) as string,
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
    _idBytesLen: number,
    count: number,
    cores: number,
    scanPpm: number,
    scanMz: number,
    eicPpm: number,
    eicMz: number,
    packedOpts: Uint8Array | undefined,
  ): any {
    const ids = unpackIds(count, offsets, lengths, idBytes);
    return JSON.parse(
      this.native.findFeature(
        handle,
        rts,
        mzs,
        ranges,
        ids,
        cores,
        scanPpm,
        scanMz,
        eicPpm,
        eicMz,
        packedOpts ? toNodeBuffer(packedOpts) : null,
      ) as string,
    );
  }

  getFeatures(
    dirPath: string,
    from: number,
    to: number,
    eicPpm: number,
    eicMz: number,
    gStart: number,
    gEnd: number,
    gStep: number,
    groupPpm: number,
    groupMz: number,
    groupRt: number,
    prevalence: number,
    packedOpts: Uint8Array | undefined,
    cores: number,
    useGpu: number,
    batchSize: number,
  ): any {
    return JSON.parse(
      this.native.getFeatures(
        dirPath,
        from,
        to,
        eicPpm,
        eicMz,
        gStart,
        gEnd,
        gStep,
        groupPpm,
        groupMz,
        groupRt,
        prevalence,
        packedOpts ? toNodeBuffer(packedOpts) : null,
        cores,
        useGpu,
        batchSize,
      ) as string,
    );
  }
}
