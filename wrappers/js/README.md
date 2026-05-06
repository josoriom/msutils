# msutils-js

JavaScript wrappers for the msutils library for LC-MS data processing.

## Install

```bash
npm install msutils
```

## Concepts

**Sample**: A parsed mzML or ion file held in memory. Create samples with `parseMzML()` or `parseIon()`. Samples are automatically garbage collected when no longer referenced.

**XY data**: A pair of numeric arrays where `x` is retention time and `y` is intensity.

## Quick start

Load a file and extract peaks in a few lines:

```js
import * as msutils from "msutils";

const raw = await fetch("sample.mzML").then((r) => r.arrayBuffer());

const sample = await msutils.parseMzML(raw);
const eic = msutils.calculateEic(sample, 174.112, { from: 1, to: 15 });
const peaks = msutils.findPeaks(eic.x, eic.y, {
  autoNoise: true,
  snRatio: 3.0,
});

peaks.forEach((p) => {
  console.log(`rt=${p.rt.toFixed(3)}  intensity=${p.intensity.toFixed(0)}`);
});
```

## Load files

### Load mzML

Load an mzML file into a sample.

```js
const sample = await msutils.parseMzML(buffer);
```

### Load ion

Load an ion binary file into a sample.

```js
const sample = await msutils.parseIon(buffer);
const sample = await msutils.parseIon(buffer, { maxCacheSize: 1000000 });
```

## Convert

### Sample to JSON

Get JSON from a sample.

```js
const json = msutils.ionToJson(sample);
```

### Sample to mzML

Get mzML text from a sample.

```js
const mzml = msutils.ionToMzml(sample);
```

### Sample to ion

Convert a sample to compressed ion bytes.

```js
const ionBytes = msutils.mzmlToIon(sample, { level: 12, f32Compress: false });
```

## Extract from samples

### Extract EIC

Get an extracted ion chromatogram for one m/z from a sample.

```js
const eic = msutils.calculateEic(
  sample,
  174.112,
  { from: 1, to: 15 },
  20.0, // ppm tolerance
  0.005, // m/z tolerance
);
// eic.x and eic.y are Float64Array
```

### Get scans

Get scans from a sample by retention time range, m/z range, or closest value.

```js
// Get scans in RT range
const scans = msutils.getScans(sample, { rt: { from: 1, to: 15 } }, 1);

// Get scan closest to RT
const scan = msutils.getScans(sample, { rt: { closest: 5.3 } }, 1);

// Get scans in m/z range
const scans = msutils.getScans(
  sample,
  { selectedMz: { from: 100, to: 500 } },
  1,
);

// Get scan closest to m/z
const scan = msutils.getScans(sample, { selectedMz: { closest: 174.112 } }, 1);
```

## Process XY data

### Calculate baseline

Estimate a baseline for XY data.

```js
const baseline = msutils.calculateBaseline(y, {
  baselineWindow: 101,
  baselineWindowFactor: 3,
});
```

### Find noise level

Estimate the noise level in XY data.

```js
const noise = msutils.findNoiseLevel(y);
```

### Find peaks

Find all peaks in XY data.

```js
const peaks = msutils.findPeaks(x, y, {
  autoNoise: true,
  snRatio: 3.0,
});
```

### Get single peak

Find one peak near a target retention time in XY data.

```js
const peak = msutils.getPeak(x, y, 5.3, 0.6, { autoNoise: true });
```

## Find peaks from samples

### From EICs

Find peaks for many targets using extracted ion chromatograms from a sample.

```js
const targets = [
  { rt: 5.3, mz: 174.112, range: 0.5, id: "compound_A" },
  { rt: 7.1, mz: 203.156, range: 0.5, id: "compound_B" },
];

const peaks = msutils.getPeaksFromEic(
  sample,
  targets,
  { from: 0.5, to: 10.0 },
  { autoNoise: true },
  2, // cores
);
```

### From chromatograms

Find peaks from stored chromatograms in a sample.

```js
const items = [
  { idx: 0, rt: 5.3, window: 0.5 },
  { idx: 5, rt: 7.1, window: 0.5 },
];

const peaks = msutils.getPeaksFromChrom(sample, items, { autoNoise: true }, 2);
```

## Untargeted

### Feature detection

Find targeted features by m/z and retention time from a sample

```js
const targets = [
  { rt: 5.3, mz: 174.112, range: 0.5, id: "A" },
  { rt: 7.1, mz: 203.156, range: 0.5, id: "B" },
];

const features = msutils.findFeature(sample, targets, {
  scanEic: { ppmTolerance: 10, mzTolerance: 0.003 },
  eic: { ppmTolerance: 20, mzTolerance: 0.005 },
  findPeak: { autoNoise: true },
  cores: 2,
});
```

### From single sample

Find all features in a sample (Node.js only).

```js
const features = msutils.findFeatures(
  sample,
  { from: 1, to: 15 },
  {
    eic: { ppmTolerance: 10, mzTolerance: 0.005 },
    grid: { start: 40, end: 1000, stepSize: 0.005 },
    findPeak: { autoNoise: true, snRatio: 3 },
    cores: 4,
  },
);
```

### Find and align features across samples

Find and align features across many samples (Node.js only).

```js
const features = msutils.getFeatures(
  "/path/to/ion/dir",
  { from: 0, to: 20 },
  {
    eic: { ppmTolerance: 5 },
    grid: { start: 50, end: 1000, stepSize: 0.01 },
    grouping: {
      ppmTolerance: 5,
      mzTolerance: 0.0025,
      rtTolerance: 0.05,
      frequency: 2,
    },
    cores: 8,
  },
);
```

## Peak detection options

Pass peak detection options as an object:

```js
const options = {
  integralThreshold: 0.0,
  intensityThreshold: 150.0,
  widthThreshold: 5,
  noise: NaN,
  autoNoise: true,
  autoBaseline: true,
  baselineWindow: 0,
  baselineWindowFactor: 0,
  allowOverlap: false,
  windowSize: 0,
  snRatio: 3.0,
};
```
