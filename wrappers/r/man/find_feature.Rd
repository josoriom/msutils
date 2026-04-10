\name{find_feature}
\alias{find_feature}
\title{Targeted feature extraction}
\description{
Finds one or more target features from known m/z and RT values.
}
\usage{
find_feature(
  bin,
  rt, mz, window, id = NULL,
  scan_ppm, scan_mz, eic_ppm, eic_mz,
  cores = 1L,
  integral_threshold = NaN,
  intensity_threshold = 1000,
  width_threshold = 5L,
  noise = NaN,
  auto_noise = TRUE,
  auto_baseline = TRUE,
  baseline_window = 0L,
  baseline_window_factor = 0L,
  allow_overlap = FALSE,
  window_size = 0L,
  sn_ratio = 1
)
}
\arguments{
  \item{bin}{A file handle returned by \code{\link{parse_mzml}()} or \code{\link{parse_bin}()}.}
  \item{rt}{Target retention time or vector of RT values.}
  \item{mz}{Target m/z or vector of m/z values.}
  \item{window}{Target RT window or vector of windows.}
  \item{id}{Optional feature identifiers.}
  \item{scan_ppm}{PPM tolerance used while scanning spectra.}
  \item{scan_mz}{Absolute m/z tolerance used while scanning spectra.}
  \item{eic_ppm}{PPM tolerance used for EIC extraction.}
  \item{eic_mz}{Absolute m/z tolerance used for EIC extraction.}
  \item{cores}{Number of CPU cores to use.}
  \item{integral_threshold}{Minimum peak area threshold.}
  \item{intensity_threshold}{Minimum intensity threshold.}
  \item{width_threshold}{Minimum peak width.}
  \item{noise}{Manual noise estimate.}
  \item{auto_noise}{Logical flag to estimate noise automatically.}
  \item{auto_baseline}{Logical flag to estimate baseline automatically.}
  \item{baseline_window}{Baseline window size.}
  \item{baseline_window_factor}{Baseline window factor.}
  \item{allow_overlap}{Logical flag to allow overlapping peaks.}
  \item{window_size}{Optional smoothing or local window size.}
  \item{sn_ratio}{Signal-to-noise ratio threshold.}
}
\value{
A data frame of extracted targeted features.
}
\details{
Use this when you already know the target m/z and RT values you want to check.
}
\seealso{
\code{\link{find_features}()}, \code{\link{get_peaks_from_eic}()}
}
