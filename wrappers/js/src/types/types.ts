export type BinaryInput = Uint8Array | ArrayBuffer;
export type IonInput = BinaryInput | URL | string;

export type PeakOptions = Partial<{
  minIntegral: number;
  minIntensity: number;
  minPeakWidthPoints: number;
  noise: number;
  autoNoise: boolean | number;
  autoBaseline: boolean | number;
  lambda: number;
  maxIterations: number;
  allowOverlap: boolean | number;
  minSnr: number;
  minGaussianR2: number;
}>;

export type BaselineOptions = Partial<{
  lambda: number;
  maxIterations: number;
}>;

export type Peak = {
  from: number;
  to: number;
  rt: number;
  integral: number;
  intensity: number;
  nPoints: number;
  noise: number;
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
  totalArea: number;
  timestamp: string;
};

export type Feature = {
  mz: number;
  rt: number;
  from: number;
  to: number;
  intensity: number;
  integral: number;
  nPoints: number;
};

export type FoundFeature = {
  id: string;
  mz: number;
  rt: number;
  from: number;
  to: number;
  intensity: number;
  integral: number;
  nPoints: number;
  noise: number;
};

export type FindFeaturesOptions = {
  eic?: { ppmTolerance?: number; mzTolerance?: number };
  grid?: { start?: number; end?: number; stepSize?: number };
  findPeak?: PeakOptions;
  cores?: number;
  useGpu?: boolean;
  batchSize?: number;
};

export type ConsensusFeature = {
  mz: number;
  rt: number;
  from: number;
  to: number;
  intensity: number;
  integral: number;
  frequency: number;
};

export type GetFeaturesOptions = FindFeaturesOptions & {
  grouping?: {
    ppmTolerance?: number;
    mzTolerance?: number;
    rtTolerance?: number;
    frequency?: number;
  };
};

export type FromTo = {
  from: number;
  to: number;
};

export type SpectrumSummary = {
  rtSeconds: number;
  basePeakMz: number;
  selectedIonMz: number;
  basePeakInt: number;
  totalIonCurrent: number;
  msLevel: number;
  polarity: number;
  positionX: number;
  positionY: number;
  positionZ: number;
};

export type CentroidScan = {
  rt: number;
  mz: number[];
  intensity: number[];
  metadata: SpectrumSummary;
};

export type ScanQuery =
  | { rt: { from: number; to: number } }
  | { rt: { closest: number } }
  | { mz: { from: number; to: number } }
  | { mz: { closest: number } };

export type IonImageOptions = {
  tolerance?: number;
  level?: number;
  onProgress?: (fetched: number, total: number, heldBytes: number) => void;
};

export type IonImage = {
  width: number;
  height: number;
  minX: number;
  minY: number;
  minZ: number;
  maxZ: number;
  data: number[];
  counts: number[];
};

export type ByteRange = {
  offset: number;
  length: number;
};
