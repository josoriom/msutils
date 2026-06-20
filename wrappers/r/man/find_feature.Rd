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
  min_gaussian_r2 = 0.3
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
  \item{min_gaussian_r2}{Minimum Gaussian R-squared for the artifact filter; \code{0} disables it (keep peaks regardless of shape). Default 0.3.}
}
\value{
A data frame with one row per target. Columns include \code{id}, \code{mz},
\code{rt}, \code{from}, \code{to}, \code{intensity}, \code{integral},
\code{n_points}, \code{noise}.
}
\seealso{
\code{\link{find_features}()}, \code{\link{get_peaks_from_eic}()}
}
