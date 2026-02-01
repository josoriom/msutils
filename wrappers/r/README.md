# msutils

## Install

Using **remotes**:

```r
remotes::install_github("josoriom/msutils", subdir = "wrappers/r")
```

## Parse an mzML (parse_mzml)

```r
path <- "/path/to/file.mzML"
bin <- readBin(path, "raw", file.info(path)$size)
file <- msutils::parse_mzml(bin)
```

## calculate_eic

```py
target <- 288.1326
from_time <- 300
to_time <- 600
ppm_tol <- 10
mz_tol <- 0.005
xy <- msutils::calculate_eic(file, target, from_time, to_time, ppm_tol, mz_tol)
x <- xy$x
y <- xy$y
```

## calculate_baseline

```py
xs <- seq(0, 10, length.out = 1000)
ys <- sin(xs) + 0.1 * xs + 2
baseline <- msutils::calculate_baseline(ys, baseline_window = 101, baseline_window_factor = 3)
```

## find_peaks

```py
xs <- seq(0, 10, length.out = 2000)
set.seed(1)
ys <- (sin(xs * 2.5))^2 + 0.05 * rnorm(length(xs)) + 0.2
peaks <- msutils::find_peaks(xs, ys, options = list(window_size = 17, sn_ratio = 2.0, auto_noise = TRUE,   allow_overlap = FALSE))
```

## get_peak

```py
xs <- seq(0, 10, length.out = 2000)
set.seed(2)
ys <- exp(-(xs - 5)^2 / (2 * 0.2^2)) + 0.02 * rnorm(length(xs))
peak <- msutils::get_peak(xs, ys, rt = 5.0, range = 0.6, options = list(window_size = 17, sn_ratio = 1.5))
```

## get_peaks_from_chrom

```py
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

```py
path <- "path/to/your/file.mzML"
targets <- data.frame(
  id = c(
    "1-methylhistidine",
    "Leucine",
    "Glutamine",
    "Alanine",
    "Arginine",
    "Asparagine"
  ),
  rt = c(2.61, 4.15, 1.91, 2.39, 2.59, 1.75),
  mz = c(340.1404, 302.1499, 317.1244, 260.1030, 345.1670, 303.1088),
  ranges = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2),
  stringsAsFactors = FALSE
)
bin <- readBin(path, "raw", file.info(path)$size)
file <- msutils::parse_mzml(bin)
peaks <- msutils::get_peaks_from_eic(file, targets, cores = 2)
```
