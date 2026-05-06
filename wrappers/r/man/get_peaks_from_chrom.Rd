\name{get_peaks_from_chrom}
\alias{get_peaks_from_chrom}
\title{Find peaks from stored chromatograms}
\description{
Finds peaks from file-stored chromatograms using index, RT, and RT window.
Use this when you already know which chromatogram index to read.
}
\usage{
get_peaks_from_chrom(
  file, items, cores = 1L,
  integral_threshold = NaN,
  intensity_threshold = 1000,
  width_threshold = 5L,
  noise = NaN,
  auto_noise = FALSE,
  auto_baseline = FALSE,
  baseline_window = 0L,
  baseline_window_factor = 0L,
  allow_overlap = FALSE,
  window_size = 0L,
  sn_ratio = 1
)
}
\arguments{
  \item{file}{Sample file from \code{\link{parse_mzml}()} or \code{\link{parse_ion}()}.}
  \item{items}{Data frame or list with columns \code{idx} (chromatogram index, integer),
    \code{rt} (target retention time), and \code{window} or \code{range} (RT half-window).}
  \item{cores}{Number of CPU cores to use. Default 1.}
  \item{integral_threshold}{Minimum peak area. \code{NaN} disables.}
  \item{intensity_threshold}{Minimum apex intensity. Default 1000.}
  \item{width_threshold}{Minimum peak width in points. Default 5.}
  \item{noise}{Fixed noise level. \code{NaN} disables.}
  \item{auto_noise}{Estimate noise automatically. Default \code{FALSE}.}
  \item{auto_baseline}{Subtract rolling baseline before detection. Default \code{FALSE}.}
  \item{baseline_window}{Baseline window size. \code{0} uses a default.}
  \item{baseline_window_factor}{Scales the baseline window. \code{0} uses a default.}
  \item{allow_overlap}{Allow overlapping peak boundaries. Default \code{FALSE}.}
  \item{window_size}{Savitzky-Golay smoothing window. \code{0} selects automatically.}
  \item{sn_ratio}{Minimum signal-to-noise ratio. Default 1.}
}
\value{
A data frame with one row per input item. Columns include \code{index}, \code{rt},
\code{from}, \code{to}, \code{intensity}, \code{integral}.
}
\examples{
\donttest{
items <- data.frame(
  idx    = c(0L, 5L),
  rt     = c(1.2, 2.85),
  window = c(0.4, 0.5)
)
get_peaks_from_chrom(sample, items, cores = 2L)
}
}
\seealso{
\code{\link{get_peaks_from_eic}()}, \code{\link{get_peak}()}
}
