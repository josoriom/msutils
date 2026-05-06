\name{find_feature}
\alias{find_feature}
\title{Find targeted features by m/z and RT}
\description{
Extracts peaks for one or more targets defined by m/z and RT.
Use this when you already know what you are looking for.
}
\usage{
find_feature(
  file,
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
  \item{file}{Sample file from \code{\link{parse_mzml}()} or \code{\link{parse_ion}()}.}
  \item{rt}{Target RT or numeric vector of RTs (minutes).}
  \item{mz}{Target m/z or numeric vector of m/z values (Da).}
  \item{window}{RT search half-window or numeric vector of half-windows.}
  \item{id}{Optional character vector of feature labels. \code{NULL} uses blank labels.}
  \item{scan_ppm}{PPM tolerance for spectrum scan extraction.}
  \item{scan_mz}{Absolute m/z tolerance (Da) for spectrum scan extraction.}
  \item{eic_ppm}{PPM tolerance for EIC extraction.}
  \item{eic_mz}{Absolute m/z tolerance (Da) for EIC extraction.}
  \item{cores}{Number of CPU cores to use. Default 1.}
  \item{integral_threshold}{Minimum peak area. \code{NaN} disables.}
  \item{intensity_threshold}{Minimum apex intensity. Default 1000.}
  \item{width_threshold}{Minimum peak width in points. Default 5.}
  \item{noise}{Fixed noise level. \code{NaN} disables.}
  \item{auto_noise}{Estimate noise automatically. Default \code{TRUE}.}
  \item{auto_baseline}{Subtract rolling baseline before detection. Default \code{TRUE}.}
  \item{baseline_window}{Baseline window size. \code{0} uses a default.}
  \item{baseline_window_factor}{Scales the baseline window. \code{0} uses a default.}
  \item{allow_overlap}{Allow overlapping peak boundaries. Default \code{FALSE}.}
  \item{window_size}{Savitzky-Golay smoothing window. \code{0} selects automatically.}
  \item{sn_ratio}{Minimum signal-to-noise ratio. Default 1.}
}
\value{
A data frame with one row per target. Columns include \code{id}, \code{mz},
\code{rt}, \code{from}, \code{to}, \code{intensity}, \code{integral}.
}
\seealso{
\code{\link{find_features}()}, \code{\link{get_peaks_from_eic}()}
}
