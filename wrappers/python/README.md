# msutils-py

## Install

```bash
pip install git+https://github.com/josoriom/msutils#subdirectory=wrappers/python
```

## Usage

```python
import msutils

raw = open("sample.mzML", "rb").read()

with msutils.parse_mzml(raw) as f:
    eic   = msutils.calculate_eic(f, target_mz=174.112, from_rt=1, to_rt=15)
    peaks = msutils.find_peaks(eic["x"], eic["y"], options={"auto_noise": True, "sn_ratio": 3.0})
    for p in peaks:
        print(f"rt={p['rt']:.3f}  intensity={p['intensity']:.0f}")
```

The library loads automatically on import. No initialisation call needed.

## Full API

### Parsing

```python
f = msutils.parse_mzml(bytes)
f = msutils.parse_bin(bytes)
```

### Conversion

```python
msutils.bin_to_json(f)
msutils.convert_bin_to_mzml(f)
msutils.mzml_to_bin(f, level=5, f32_compress=False)
```

### EIC / scans

```python
msutils.calculate_eic(f, target_mz, from_rt, to_rt, ppm_tol=20.0, mz_tol=0.005)
msutils.collect_scans(f, from_rt, to_rt, level=1)
```

### Peak detection

```python
msutils.find_peaks(x, y, options=None)
msutils.get_peak(x, y, rt, range_, options=None)
msutils.find_noise_level(y)
msutils.calculate_baseline(y, baseline_window=0, baseline_window_factor=0)
```

### Batch peak extraction

```python
targets = [{"rt": 5.3, "mz": 174.112, "range": 0.5, "id": "A"}]

msutils.get_peaks_from_eic(f, targets, from_rt=0.5, to_rt=5.0, options=None, cores=1)
msutils.get_peaks_from_chrom(f, items, options=None, cores=1)
```

### Feature detection

```python
msutils.find_features(
    f, from_rt=1, to_rt=15,
    eic     = {"ppm_tolerance": 10, "mz_tolerance": 0.005},
    grid    = {"start": 40, "end": 1000, "step_size": 0.005},
    options = {"auto_noise": True, "sn_ratio": 3},
    cores   = 4,
)

msutils.find_feature(
    f, targets,
    scan_eic = {"ppm_tolerance": 10, "mz_tolerance": 0.003},
    eic      = {"ppm_tolerance": 20, "mz_tolerance": 0.005},
    options  = {"auto_noise": True},
    cores    = 2,
)

msutils.get_features(
    "/path/to/bin/dir", from_rt=0, to_rt=20,
    eic      = {"ppm_tolerance": 5},
    grid     = {"start": 50, "end": 1000, "step_size": 0.01},
    grouping = {"ppm_tolerance": 5, "rt_tolerance": 0.05, "frequency": 2},
    cores    = 8,
)
```

### Peak options

The `options` parameter accepts a plain dict:

```python
{
    "integral_threshold":     0.0,
    "intensity_threshold":    150.0,
    "width_threshold":        5,
    "noise":                  float("nan"),
    "auto_noise":             True,
    "auto_baseline":          True,
    "baseline_window":        0,
    "baseline_window_factor": 0,
    "allow_overlap":          False,
    "window_size":            0,
    "sn_ratio":               3.0,
}
```

### MzMlFile instance methods

```python
with msutils.parse_mzml(data) as f:
    f.to_json()
    f.to_mzml()
    f.to_bin(level=12)
    f.dispose()
```
