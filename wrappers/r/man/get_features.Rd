\name{get_features}
\alias{get_features}
\title{Find and align features across many files}
\description{
Detects features in every file in a folder and groups them by m/z and RT.
Returns a consensus feature table aligned across all samples.
}
\usage{
get_features(
  dir_path,
  from = 0,
  to = 10,
  eic_ppm_tol = 10,
  eic_mz_tol = 0.005,
  grid_start = 40,
  grid_end = 1000,
  grid_step = 0.005,
  group_ppm_tol = 5,
  group_da_tol = 0.0025,
  group_rt_tol = 0.05,
  frequency = 1L,
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
  \item{dir_path}{Path to a folder with mzML or ion files.}
  \item{from}{Start retention time (minutes). Default 0.}
  \item{to}{End retention time (minutes). Default 10.}
  \item{eic_ppm_tol}{PPM tolerance for EIC extraction. Default 10.}
  \item{eic_mz_tol}{Absolute m/z tolerance (Da) for EIC extraction. Default 0.005.}
  \item{grid_start}{Lower m/z bound of the search grid (Da). Default 40.}
  \item{grid_end}{Upper m/z bound of the search grid (Da). Default 1000.}
  \item{grid_step}{Grid step size (ppm). Default 0.005.}
  \item{group_ppm_tol}{PPM tolerance used when grouping features across samples. Default 5.}
  \item{group_da_tol}{Absolute m/z tolerance (Da) used when grouping features. Default 0.0025.}
  \item{group_rt_tol}{RT tolerance (minutes) used when grouping features. Default 0.05.}
  \item{frequency}{Minimum number of samples in which a feature must appear to be kept. Default 1.}
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
A data frame with one row per consensus feature. Columns include \code{mz},
\code{rt}, \code{from}, \code{to}, \code{intensity}, \code{integral},
\code{np}, \code{frequency} (number of samples), \code{rmz} (representative m/z).
}
\examples{
\donttest{
features <- get_features(
  "/path/to/samples",
  from = 0, to = 10,
  cores = 4L,
  frequency = 2L
)
}
}
\seealso{
\code{\link{find_features}()}, \code{\link{parse_mzml}()}, \code{\link{parse_ion}()}
}
