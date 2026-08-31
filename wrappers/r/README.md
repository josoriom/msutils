# quantion

R wrappers for the quantion library for LC-MS data processing.

## Install

```r
remotes::install_github("josoriom/quantion", subdir = "wrappers/r")
```

## Concepts

**Sample**: A pointer to a parsed mzML file held in memory, or to an ion file decoded on demand. Create samples with `parse_mzml()` or `parse_ion()`. Memory is freed automatically when no longer in use, or call `dispose()` to free manually.

**XY data**: A pair of numeric vectors where `x` is retention time and `y` is intensity.

## Load files

### mzML

Load an mzML file into a sample.

```r
mzml_path <- "/path/to/file.mzML"
buf <- readBin(mzml_path, "raw", file.info(mzml_path)$size)
sample <- quantion::parse_mzml(buf)
```

### ion

Load an ion binary file into a sample.

```r
ion_path <- "/path/to/file.ion"
buf <- readBin(ion_path, "raw", file.info(ion_path)$size)
sample <- quantion::parse_ion(buf)
```

## Convert

### mzML to ion

Convert a sample to compressed ion bytes.

```r
out <- quantion::mzml_to_ion(sample, level = 12L, f32_compress = FALSE)
writeBin(out, "/path/to/file.ion")
```

### Sample to mzML

Get mzML text from a sample.

```r
mzml <- quantion::ion_to_mzml(sample)
```

## Extract from samples

### Extract EIC

Get an extracted ion chromatogram for one m/z from a sample.

```r
eic <- quantion::calculate_eic(sample, targets = 288.1326, from = 2.5, to = 4.0, ppm_tolerance = 10, mz_tolerance = 0.005)
eic_df <- data.frame(x = eic$x, y = eic$y)
```

### Get scans

Get scans from a sample by retention time or m/z query.

```r
scans <- quantion::get_scans(sample, rt_from = 2.5, rt_to = 4.0, level = 1L)
```

## Process XY data

### Calculate baseline

Estimate a baseline for XY data.

```r
xs <- seq(0, 10, length.out = 1000)
ys <- sin(xs) + 0.1 * xs + 2
baseline <- quantion::calculate_baseline(ys, lambda = 0L, max_iterations = 0L)
```

### Find peaks

Find all peaks in XY data.

```r
peaks <- quantion::find_peaks(eic_df$x, eic_df$y, min_intensity = 1000, min_peak_width_points = 10L, auto_noise = TRUE, auto_baseline = TRUE)
```

### Get single peak

Find one peak near a target retention time in XY data.

```r
peak <- quantion::get_peak(eic_df$x, eic_df$y, rt = 3.4, range = 0.6, min_intensity = 1000, min_peak_width_points = 10L, auto_noise = TRUE, auto_baseline = TRUE)
```

## Find peaks from samples

### From chromatograms

Find peaks from stored chromatograms in a sample.

```r
transitions <- data.frame(
  idx    = c(0L, 5L, 18L),
  rt     = c(1.20, 2.85, 3.50),
  window = c(0.40, 0.50, 0.60)
)
peaks <- quantion::get_peaks_from_chrom(sample, transitions, cores = 2L)
```

### From EICs

Find peaks for many targets using extracted ion chromatograms from a sample.

```r
targets <- data.frame(
  id     = c("1-methylhistidine", "Leucine", "Glutamine", "Alanine", "Arginine", "Asparagine"),
  rt     = c(2.61, 4.15, 1.91, 2.39, 2.59, 1.75),
  mz     = c(340.1404, 302.1499, 317.1244, 260.1030, 345.1670, 303.1088),
  ranges = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2)
)
peaks <- quantion::get_peaks_from_eic(sample, targets, auto_noise = TRUE, auto_baseline = TRUE, cores = 2L)
```

## Untargeted

### Feature detection

Find features by m/z and retention time from a sample.

```r
targets <- data.frame(
  id     = c("compound1", "compound2"),
  rt     = c(2.61, 4.15),
  mz     = c(340.1404, 302.1499),
  ranges = c(0.2, 0.2)
)
features <- quantion::find_feature(sample, targets, auto_noise = TRUE, auto_baseline = TRUE)
```

### From single sample

Find all features in a sample.

```r
features <- quantion::find_features(
  sample,
  from = 0, to = 10,
  cores = 5L,
  auto_noise = TRUE,
  auto_baseline = TRUE,
  min_peak_width_points = 10L,
  min_intensity = 1000
)
```

### Align features across samples

Find and align features across many samples.

```r
features <- quantion::get_features(
  "/path/to/sample/directory",
  from = 0,
  to = 10,
  cores = 5L,
  auto_noise = TRUE,
  auto_baseline = TRUE
)
```
