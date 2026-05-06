\name{get_peaks_from_eic}
\alias{get_peaks_from_eic}
\title{Find peaks for many targets using EICs}
\description{
Builds an EIC for each target and finds the best peak near the given RT.
Use this for batch targeted peak extraction from one file.
}
\usage{
get_peaks_from_eic(
  file, df, from = 0.5, to = 5, cores = 1L,
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
  \item{df}{Data frame with columns \code{id} (character), \code{rt} (numeric),
    \code{mz} (numeric), \code{ranges} (numeric RT half-window).}
  \item{from}{EIC extraction start time (minutes). Default 0.5.}
  \item{to}{EIC extraction end time (minutes). Default 5.}
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
\examples{
\donttest{
targets <- data.frame(
  id     = c("Leucine", "Alanine"),
  rt     = c(4.15, 2.39),
  mz     = c(302.1499, 260.1030),
  ranges = c(0.2, 0.2)
)
get_peaks_from_eic(sample, targets, cores = 2L)
}
}
\seealso{
\code{\link{calculate_eic}()}, \code{\link{get_peak}()}, \code{\link{get_peaks_from_chrom}()}
}
