\name{find_features}
\alias{find_features}
\title{Find all features in one file}
\description{
Scans the full m/z grid and finds all chromatographic features in one file.
No target list is needed.
}
\usage{
find_features(
  file,
  from = 0,
  to = 10,
  ppm_tolerance = 10,
  mz_tolerance = 0.005,
  grid_start = 40,
  grid_end = 1000,
  grid_step_ppm = 0.005,
  cores = 1L,
  use_gpu = FALSE,
  batch_size = 0L,
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
  \item{from}{Start retention time (minutes). Default 0.}
  \item{to}{End retention time (minutes). Default 10.}
  \item{ppm_tolerance}{PPM tolerance for EIC extraction.}
  \item{mz_tolerance}{Absolute m/z tolerance (Da) for EIC extraction.}
  \item{grid_start}{Lower m/z bound of the search grid (Da). Default 40.}
  \item{grid_end}{Upper m/z bound of the search grid (Da). Default 1000.}
  \item{grid_step_ppm}{Grid step size (ppm). Default 0.005.}
  \item{cores}{Number of CPU cores to use. Default 1.}
  \item{use_gpu}{Use GPU acceleration if available. Default \code{FALSE}.}
  \item{batch_size}{GPU batch size. \code{0} selects automatically.}
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
A data frame with one row per detected feature. Columns include \code{mz},
\code{rt}, \code{from}, \code{to}, \code{intensity}, \code{integral},
\code{n_points}.
}
\examples{
\donttest{
features <- find_features(
  file,
  from = 0,
  to = 10,
  cores = 5L,
  auto_noise = TRUE,
  auto_baseline = TRUE,
  min_intensity = 1000
)
}
}
\seealso{
\code{\link{find_feature}()}, \code{\link{get_features}()}
}
