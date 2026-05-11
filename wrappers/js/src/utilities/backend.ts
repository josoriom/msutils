export type FileHandle = unknown;

export interface Backend {
  readonly ready: boolean;
  readonly handlesAreGcFinalized: boolean;

  init(): Promise<void>;

  parseMzML(data: Uint8Array): FileHandle;
  parseBin(data: Uint8Array, maxCacheSize?: number): FileHandle;
  freeFile(handle: FileHandle): void;

  fileToJson(handle: FileHandle): any;
  fileToMzml(handle: FileHandle): string;
  fileToBin(
    handle: FileHandle,
    level: number,
    f32Compress: boolean,
  ): Uint8Array;

  calculateEic(
    handle: FileHandle,
    targetMz: number,
    from: number,
    to: number,
    ppmTol: number,
    mzTol: number,
  ): { x: Float64Array; y: Float64Array };

  getScans(
    handle: FileHandle,
    queryType: number,
    a: number,
    b: number,
    level: number,
  ): any;

  findPeaks(
    x: Float64Array,
    y: Float64Array,
    packedOpts: Uint8Array | undefined,
  ): any;

  getPeak(
    x: Float64Array,
    y: Float64Array,
    rt: number,
    range: number,
    packedOpts: Uint8Array | undefined,
  ): any;

  findNoiseLevel(y: Float64Array | Float32Array): number;

  calculateBaseline(
    y: Float64Array,
    baselineWindow: number,
    baselineWindowFactor: number,
  ): Float64Array;

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
    packedOpts: Uint8Array | undefined,
    cores: number,
  ): any;

  getPeaksFromChrom(
    handle: FileHandle,
    indices: Uint32Array,
    rts: Float64Array,
    windows: Float64Array,
    count: number,
    packedOpts: Uint8Array | undefined,
    cores: number,
  ): any;

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
  ): any;

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
    packedOpts: Uint8Array | undefined,
  ): any;

  getFeatures?(
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
  ): any;
}
