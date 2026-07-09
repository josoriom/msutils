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
  min_integral = NaN,
  min_intensity = 1000,
  min_peak_width_points = 5L,
  noise = NaN,
  auto_noise = TRUE,
  auto_baseline = TRUE,
  lambda = 0L,
  max_iterations = 0L,
  allow_overlap = FALSE,
  min_snr = 1,
  min_r2 = 0.3
)
}
\arguments{
  \item{file}{Sample file from \code{\link{parse_mzml}()} or \code{\link{parse_ion}()}.}
  \item{df}{Data frame with columns \code{id} (character), \code{rt} (numeric),
    \code{mz} (numeric), \code{ranges} (numeric RT half-window).}
  \item{from}{EIC extraction start time (minutes). Default 0.5.}
  \item{to}{EIC extraction end time (minutes). Default 5.}
  \item{cores}{Number of CPU cores to use. Default 1.}
  \item{min_integral}{Minimum peak area. \code{NaN} disables.}
  \item{min_intensity}{Minimum apex intensity. Default 1000.}
  \item{min_peak_width_points}{Minimum peak width in points. Default 5.}
  \item{noise}{Fixed noise level. \code{NaN} disables.}
  \item{auto_noise}{Estimate noise automatically. Default \code{TRUE}.}
  \item{auto_baseline}{Subtract baseline before detection. Default \code{TRUE}.}
  \item{lambda}{airPLS smoothing parameter for baseline estimation. \code{0} uses a default.}
  \item{max_iterations}{Maximum airPLS iterations. \code{0} uses a default.}
  \item{allow_overlap}{Allow overlapping peak boundaries. Default \code{FALSE}.}
  \item{min_snr}{Minimum signal-to-noise ratio. Default 1.}
  \item{min_r2}{Minimum fit R-squared (Gaussian or EMG) for the artifact filter; \code{0} disables it (keep peaks regardless of shape). Default 0.3.}
}
\value{
A data frame with one row per target. Columns include \code{id}, \code{mz},
\code{rt}, \code{from}, \code{to}, \code{intensity}, \code{integral},
\code{noise}.
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
