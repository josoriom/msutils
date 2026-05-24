# msutils-py

Python wrappers for the msutils Rust library for LC-MS data processing.

## Install

```bash
pip install git+https://github.com/josoriom/msutils#subdirectory=wrappers/python
```

## Concepts

**Sample**: A parsed mzML or ion file held in memory. Create samples with `parse_mzml()` or `parse_ion()`. Use as a context manager (`with`) to ensure automatic cleanup.

**XY data**: A pair of numeric sequences where `x` is retention time and `y` is intensity.

## Quick start

Load a file and extract peaks in a few lines:

```python
import msutils

raw = open("sample.mzML", "rb").read()

with msutils.parse_mzml(raw) as sample:
    eic   = msutils.calculate_eic(sample, target_mz=174.112, from_rt=1, to_rt=15)
    peaks = msutils.find_peaks(eic["x"], eic["y"], options={"auto_noise": True, "min_snr": 3.0})
    for p in peaks:
        print(f"rt={p['rt']:.3f}  intensity={p['intensity']:.0f}")
```

The library loads automatically on import. No initialization call needed.

## Load files

### Load mzML

Load an mzML file into a sample.

```python
sample = msutils.parse_mzml(bytes)
```

### Load ion

Load an ion binary file into a sample.

```python
sample = msutils.parse_ion(bytes)
```

## Convert

### Sample to JSON

Get JSON from a sample.

```python
json = msutils.ion_to_json(sample)
```

### Sample to mzML

Get mzML text from a sample.

```python
mzml = msutils.ion_to_mzml(sample)
```

### Sample to ion

Convert a sample to compressed ion bytes.

```python
ion_bytes = msutils.mzml_to_ion(sample, level=5, f32_compress=False)
```

### Sample instance methods

Convert using sample methods.

```python
with msutils.parse_mzml(data) as sample:
    json = sample.to_json()
    mzml = sample.to_mzml()
    ion_bytes = sample.to_ion(level=12)
    sample.dispose()
```

## Extract from samples

### Extract EIC

Get an extracted ion chromatogram for one m/z from a sample.

```python
eic = msutils.calculate_eic(sample, target_mz=174.112, from_rt=1, to_rt=15, ppm_tol=20.0, mz_tol=0.005)
```

### Get scans

Get scans from a sample by retention time range and MS level.

```python
scans = msutils.get_scans(sample, rt_range=(1, 15), level=1)
```

## Process XY data

### Calculate baseline

Estimate a baseline for XY data.

```python
baseline = msutils.calculate_baseline(y, lambda_=0, max_iterations=0)
```

### Find noise level

Estimate the noise level in XY data.

```python
noise = msutils.find_noise_level(y)
```

### Find peaks

Find all peaks in XY data.

```python
peaks = msutils.find_peaks(x, y, options={"auto_noise": True, "min_snr": 3.0})
```

### Get single peak

Find one peak near a target retention time in XY data.

```python
peak = msutils.get_peak(x, y, rt=5.3, range_=0.6, options={"auto_noise": True})
```

## Find peaks from samples

### From EICs

Find peaks for many targets using extracted ion chromatograms from a sample.

```python
targets = [
    {"rt": 5.3, "mz": 174.112, "range": 0.5, "id": "compound_A"},
    {"rt": 7.1, "mz": 203.156, "range": 0.5, "id": "compound_B"},
]

peaks = msutils.get_peaks_from_eic(
    sample, targets,
    from_rt=0.5, to_rt=10.0,
    options={"auto_noise": True},
    cores=2
)
```

### From chromatograms

Find peaks from stored chromatograms in a sample.

```python
items = [
    {"idx": 0, "rt": 5.3, "window": 0.5},
    {"idx": 5, "rt": 7.1, "window": 0.5},
]

peaks = msutils.get_peaks_from_chrom(sample, items, options=None, cores=2)
```

## Untargeted

### Targeted feature detection

Find targeted features by m/z and retention time from a sample.

```python
targets = [
    {"rt": 5.3, "mz": 174.112, "range": 0.5, "id": "A"},
    {"rt": 7.1, "mz": 203.156, "range": 0.5, "id": "B"},
]

features = msutils.find_feature(
    sample, targets,
    scan_eic={"ppm_tolerance": 10, "mz_tolerance": 0.003},
    eic={"ppm_tolerance": 20, "mz_tolerance": 0.005},
    options={"auto_noise": True},
    cores=2,
)
```

### From single sample

Find all features in a sample.

```python
features = msutils.find_features(
    sample,
    from_rt=1, to_rt=15,
    eic={"ppm_tolerance": 10, "mz_tolerance": 0.005},
    grid={"start": 40, "end": 1000, "step_size": 0.005},
    options={"auto_noise": True, "min_snr": 3},
    cores=4,
)
```

### Align features across samples

Find and align features across many samples.

```python
features = msutils.get_features(
    "/path/to/ion/dir",
    from_rt=0, to_rt=20,
    eic={"ppm_tolerance": 5},
    grid={"start": 50, "end": 1000, "step_size": 0.01},
    grouping={"ppm_tolerance": 5, "rt_tolerance": 0.05, "frequency": 2},
    cores=8,
)
```

## Peak detection options

Pass peak detection options as a dictionary:

```python
options = {
    "min_integral":        0.0,
    "min_intensity":       150.0,
    "min_peak_width_points": 5,
    "noise":               float("nan"),
    "auto_noise":          True,
    "auto_baseline":       True,
    "lambda_":             0,
    "max_iterations":      0,
    "allow_overlap":       False,
    "min_snr":             3.0,
}
```
