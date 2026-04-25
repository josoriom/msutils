export type BinaryInput = Uint8Array | ArrayBuffer;

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
  np: number;
  frequency: number;
  mzs: number[];
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
  | { selectedMz: { from: number; to: number } }
  | { selectedMz: { closest: number } };
