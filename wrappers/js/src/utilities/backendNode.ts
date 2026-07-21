/// <reference types="node" />

import * as fs from "fs";
import * as path from "path";
import type { Backend, FileHandle } from "./backend";
import type { NoiseLevel, PeakOptions } from "../types/types";
import { parseAndCamelize } from "./shared";
import {
  ByteRangeResult,
  RemoteSource,
  fetchHeader,
  newRemoteSource,
  prefetchRanges,
} from "./remoteSource";

function libraryFileName(platform: string): string {
  if (platform === "win32") return "libquantion.dll";
  if (platform === "darwin") return "libquantion.dylib";
  return "libquantion.so";
}

function platformDirName(platform: string, arch: string): string {
  const arm = arch === "arm64";
  if (platform === "darwin") return arm ? "macos-arm64" : "macos-x86_64";
  if (platform === "linux") return arm ? "linux-arm64" : "linux-x86_64";
  if (platform === "win32") return arm ? "windows-arm64" : "windows-x86_64";
  throw new Error(`quantion: unsupported ${platform}/${arch}`);
}

function isPackageRoot(directory: string): boolean {
  const manifest = path.join(directory, "package.json");
  if (!fs.existsSync(manifest)) return false;
  try {
    return JSON.parse(fs.readFileSync(manifest, "utf8")).name === "quantion";
  } catch {
    return false;
  }
}

function artifactsRoot(): string | null {
  const fromEnv = process.env.QUANTION_ARTIFACTS_ROOT;
  if (fromEnv) return fs.existsSync(fromEnv) ? fromEnv : null;
  let here = __dirname;
  for (;;) {
    const candidate = path.join(here, "artifacts");
    if (fs.existsSync(candidate)) return candidate;
    if (isPackageRoot(here)) return null;
    const parent = path.dirname(here);
    if (parent === here) return null;
    here = parent;
  }
}

function platformLibPath(proc: NodeJS.Process): string {
  const fromEnv = proc.env.QUANTION_LIB;
  if (fromEnv) return fromEnv;

  const root = artifactsRoot();
  const wanted = platformDirName(proc.platform, proc.arch);
  const file = libraryFileName(proc.platform);

  if (root) {
    const direct = path.join(root, wanted, file);
    if (fs.existsSync(direct)) return direct;

    const own = require("../../package.json").version as string;
    const preferred = path.join(root, own, wanted, file);
    if (fs.existsSync(preferred)) return preferred;

    const parts = (name: string) => name.split(".").map((p) => Number(p) || 0);
    const newest = fs
      .readdirSync(root)
      .filter((entry) => fs.existsSync(path.join(root, entry, wanted, file)))
      .filter((entry) => fs.existsSync(path.join(root, entry, "manifest.json")))
      .sort((left, right) => {
        const a = parts(left);
        const b = parts(right);
        for (let i = 0; i < Math.max(a.length, b.length); i++) {
          if ((a[i] ?? 0) !== (b[i] ?? 0)) return (b[i] ?? 0) - (a[i] ?? 0);
        }
        return 0;
      })[0];
    if (newest) return path.join(root, newest, wanted, file);
  }

  const bundled = path.join(__dirname, "..", "..", "native", wanted, file);
  if (fs.existsSync(bundled)) return bundled;

  throw new Error(
    `quantion: no library for ${wanted}. Build one with 'make ${wanted}', ` +
      `or set QUANTION_LIB to a library file, ` +
      `or QUANTION_ARTIFACTS_ROOT to the artifacts folder.`,
  );
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
  private remote_by_handle = new WeakMap<object, RemoteSource>();
  ready = false;
  readonly handlesAreGcFinalized = true;

  constructor(proc: NodeJS.Process) {
    const addonPath = path.join(
      __dirname,
      "..",
      "..",
      "build",
      "Release",
      "quantion.node",
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

  parseIonPath(path: string, cacheSize = 0): FileHandle {
    return this.native.parseIonPath(path, cacheSize);
  }


  planOpen(header: Uint8Array): ByteRangeResult[] {
    return this.native.planOpen(toNodeBuffer(header));
  }

  async parseIonRemote(url: URL, cacheSize = 0): Promise<FileHandle> {
    const source = newRemoteSource(url.href);
    const header = await fetchHeader(source);
    await prefetchRanges(source, this.planOpen(header));

    const serve = (offset: number, length: number): Buffer | null => {
      const bytes = source.cache.read(BigInt(offset), BigInt(length));
      return bytes ? toNodeBuffer(bytes) : null;
    };

    const handle = this.native.parseIonSource(serve, cacheSize) as FileHandle;
    this.remote_by_handle.set(handle as object, source);
    return handle;
  }

  parseIonBuffer(bytes: Uint8Array, cacheSize = 0): FileHandle {
    return this.native.parseBin(toNodeBuffer(bytes), cacheSize);
  }

  freeFile(handle: FileHandle): void {
    this.remote_by_handle.delete(handle as object);
    this.native.dispose(handle);
  }

  fileToJson(handle: FileHandle): any {
    return parseAndCamelize(this.native.binToJson(handle) as string);
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

  async calculateEic(
    handle: FileHandle,
    targetMz: number,
    from: number,
    to: number,
    ppmTol: number,
    mzTol: number,
  ): Promise<{ x: Float64Array; y: Float64Array }> {
    const source = this.remote_by_handle.get(handle as object);
    if (source) {
      const ranges = this.native.planEic(
        handle,
        targetMz,
        from,
        to,
        ppmTol,
        mzTol,
      ) as ByteRangeResult[];
      await prefetchRanges(source, ranges);
    }
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
    return parseAndCamelize(
      this.native.getScans(handle, queryType, a, b, level) as string,
    );
  }

  getIonImage(
    handle: FileHandle,
    mz: number,
    tolerance: number,
    level: number,
    _onProgress?: (fetched: number, total: number, heldBytes: number) => void,
  ): any {
    return parseAndCamelize(
      this.native.getIonImage(handle, mz, tolerance, level) as string,
    );
  }

  findPeaks(
    x: Float64Array,
    y: Float64Array,
    opts: PeakOptions | undefined,
  ): any {
    return parseAndCamelize(
      this.native.findPeaks(x, y, opts ?? null) as string,
    );
  }

  getPeak(
    x: Float64Array,
    y: Float64Array,
    rt: number,
    range: number,
    opts: PeakOptions | undefined,
  ): any {
    return parseAndCamelize(
      this.native.getPeak(x, y, rt, range, opts ?? null) as string,
    );
  }

  findNoiseLevel(y: Float64Array | Float32Array): NoiseLevel {
    const y32 = y instanceof Float32Array ? y : new Float32Array(y);
    return this.native.findNoiseLevel(y32) as NoiseLevel;
  }

  calculateBaseline(
    y: Float64Array,
    lambda: number,
    maxIterations: number,
  ): Float64Array {
    return this.native.calculateBaseline(
      y,
      lambda,
      maxIterations,
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
    opts: PeakOptions | undefined,
    cores: number,
  ): any {
    const ids = unpackIds(count, offsets, lengths, idBytes);
    return parseAndCamelize(
      this.native.getPeaksFromEic(
        handle,
        rts,
        mzs,
        ranges,
        ids,
        from,
        to,
        opts ?? null,
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
    opts: PeakOptions | undefined,
    cores: number,
  ): any {
    return parseAndCamelize(
      this.native.getPeaksFromChrom(
        handle,
        indices,
        rts,
        windows,
        opts ?? null,
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
    opts: PeakOptions | undefined,
    cores: number,
  ): any {
    return parseAndCamelize(
      this.native.findFeatures(
        handle,
        from,
        to,
        eicPpm,
        eicMz,
        gridStart,
        gridEnd,
        gridStep,
        opts ?? null,
        cores,
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
    opts: PeakOptions | undefined,
  ): any {
    const ids = unpackIds(count, offsets, lengths, idBytes);
    return parseAndCamelize(
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
        opts ?? null,
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
    opts: PeakOptions | undefined,
    cores: number,
  ): any {
    return parseAndCamelize(
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
        opts ?? null,
        cores,
      ) as string,
    );
  }
}
