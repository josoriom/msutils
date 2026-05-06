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
A data frame with one row per detected feature. Columns include \code{mz},
\code{rt}, \code{from}, \code{to}, \code{intensity}, \code{integral},
\code{ratio}, \code{np} (number of points).
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
  intensity_threshold = 1000
)
}
}
\seealso{
\code{\link{find_feature}()}, \code{\link{get_features}()}
}
