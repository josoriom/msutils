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
  min_integral = NaN,
  min_intensity = 1000,
  min_peak_width_points = 5L,
  noise = NaN,
  auto_noise = FALSE,
  auto_baseline = FALSE,
  lambda = 0L,
  max_iterations = 0L,
  allow_overlap = FALSE,
  min_snr = 1
)
}
\arguments{
  \item{file}{Sample file from \code{\link{parse_mzml}()} or \code{\link{parse_ion}()}.}
  \item{items}{Data frame or list with columns \code{idx} (chromatogram index, integer),
    \code{rt} (target retention time), and \code{window} or \code{range} (RT half-window).}
  \item{cores}{Number of CPU cores to use. Default 1.}
  \item{min_integral}{Minimum peak area. \code{NaN} disables.}
  \item{min_intensity}{Minimum apex intensity. Default 1000.}
  \item{min_peak_width_points}{Minimum peak width in points. Default 5.}
  \item{noise}{Fixed noise level. \code{NaN} disables.}
  \item{auto_noise}{Estimate noise automatically. Default \code{FALSE}.}
  \item{auto_baseline}{Subtract baseline before detection. Default \code{FALSE}.}
  \item{lambda}{airPLS smoothing parameter for baseline estimation. \code{0} uses a default.}
  \item{max_iterations}{Maximum airPLS iterations. \code{0} uses a default.}
  \item{allow_overlap}{Allow overlapping peak boundaries. Default \code{FALSE}.}
  \item{min_snr}{Minimum signal-to-noise ratio. Default 1.}
}
\value{
A data frame with one row per input item. Columns include \code{index}, \code{id},
\code{ort}, \code{rt}, \code{from}, \code{to}, \code{intensity}, \code{integral},
\code{total_area}, \code{timestamp}.
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
