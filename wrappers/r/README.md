# msutils

## Install

Using **remotes**:

```r
remotes::install_github("josoriom/msutils", subdir = "wrappers/r")
```

## Load files

Loads raw bytes and builds the in-memory representation. Use `parse_mzml` for .mzML and `parse_bin`() for .b64 files.

### mzML

```r
mzml_path <- "/path/to/file/file.mzML"

mzml_buffer <- readBin(mzml_path, "raw", file.info(mzml_path)$size)

mzml_file <- msutils::parse_mzml(mzml_buffer)
```

### b64

```r
b64_path  <- "/path/to/file/file.b64"

b64_buffer  <- readBin(b64_path,  "raw", file.info(b64_path)$size)

b64_file  <- msutils::parse_bin(b64_buffer)
```

## Conversion .mzML to .b64 test

Converts `mzML` files into `b64`.

```r
bin <- msutils::convert_mzml_to_bin(mzml_buffer, 12L, f32_compress = FALSE)
writeBin(bin, "/path/to/file/test_file.b64")
```

## calculate_eic

Extracts an Extracted Ion Chromatogram (EIC) for a target m/z across a retention-time interval.

```r
target <- 288.1326
from <- 2.5
to <- 4
ppm_tol <- 10
mz_tol <- 0.005

eic <- msutils::calculate_eic(b64_file, target, from, to, ppm_tol, mz_tol)
eic <- data.frame(x = eic$x, y = eic$y)

eic_plot <- ggplot(eic, aes(x = x, y = y)) +
  geom_line() +
  labs(x = "Time", y = "Intensity", title = "EIC") +
  theme_minimal()
```

## calculate_baseline

Estimates a baseline for a 1D signal

```r
xs <- seq(0, 10, length.out = 1000)
ys <- sin(xs) + 0.1 * xs + 2
baseline <- msutils::calculate_baseline(ys, baseline_window = 101, baseline_window_factor = 3)
```

## find_peaks

Detects local chromatographic peaks in a 1D trace (x=time, y=intensity).

```r
peaks <- msutils::find_peaks(eic$x, eic$y, intensity_threshold = 1000, width_threshold = 10, auto_noise = TRUE, auto_baseline = TRUE)
```

## get_peak

Targeted extraction of one peak near a specified retention time, within a given RT window.

```r
peak <- msutils::get_peak(eic$x, eic$y, rt = 3.4, range = 0.6, intensity_threshold = 1000, width_threshold = 10, auto_noise = TRUE, auto_baseline = TRUE)
```

## get_peaks_from_chrom

Work in progress!

```r
library(parallel)
transitions <- data.frame(
  idx    = c(0L, 5L, 18L),
  rt     = c(1.20, 2.85, 3.50),
  window = c(0.40, 0.50, 0.60),
  stringsAsFactors = FALSE
)
cores <- max(1L, detectCores(logical = FALSE) - 1L)
peaks <- msutils::get_peaks_from_chrom(file, transitions, cores = cores, auto_noise = FALSE)
```

## get_peaks_from_eic

Batch targeted peak detection across many features, for each feature described by {id, mz, rt, ranges}, it performs EIC extraction + peak detection.

```r
targets <- data.frame(
  id = c("1-methylhistidine", "Leucine", "Glutamine", "Alanine", "Arginine", "Asparagine"),
  rt = c(2.61, 4.15, 1.91, 2.39, 2.59, 1.75),
  mz = c(340.1404, 302.1499, 317.1244, 260.1030, 345.1670, 303.1088),
  ranges = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2),
  stringsAsFactors = FALSE
)

peaks <- msutils::get_peaks_from_eic(
  mzml_file, # or b64_file
  targets,
  auto_noise = TRUE,
  auto_baseline = TRUE,
  cores = 2L
)
```

## Find features

Untargeted feature detection across m/z–RT space.

```r
peaks <- msutils::find_features(
  mzml_file, # or b64_file
  from = 0,
  to = 10,
  cores = 5L,
  auto_noise = TRUE,
  auto_baseline = TRUE,
  width_threshold = 10L,
  intensity_threshold = 1000
)
```
