\name{get_features}
\alias{get_features}
\title{Align features across a directory of samples}
\description{
Finds features in each sample file and then aligns them across one folder.
}
\usage{
get_features(
  dir_path, from = 0, to = 10,
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
  \item{dir_path}{Path to a folder with sample files.}
  \item{from}{Start retention time.}
  \item{to}{End retention time.}
  \item{eic_ppm_tol}{PPM tolerance for EIC extraction.}
  \item{eic_mz_tol}{Absolute m/z tolerance for EIC extraction.}
  \item{grid_start}{Start m/z for the search grid.}
  \item{grid_end}{End m/z for the search grid.}
  \item{grid_step}{Grid step size.}
  \item{group_ppm_tol}{PPM tolerance used during alignment grouping.}
  \item{group_da_tol}{Absolute m/z tolerance used during alignment grouping.}
  \item{group_rt_tol}{RT tolerance used during alignment grouping.}
  \item{frequency}{Minimum number of samples required for a consensus feature.}
  \item{cores}{Number of CPU cores to use.}
  \item{use_gpu}{Logical flag to enable GPU acceleration when available.}
  \item{batch_size}{Optional GPU batch size.}
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
A data frame of aligned consensus features.
}
\details{
The folder should contain supported sample files such as \code{.mzML} or \code{.ion}. GPU is only used during per-sample feature detection.
}
\examples{
\dontrun{
aligned <- get_features(
  "/path/to/samples",
  from = 0,
  to = 10,
  cores = 4L,
  use_gpu = TRUE
)
}
}
\seealso{
\code{\link{find_features}()}, \code{\link{parse_mzml}()}, \code{\link{parse_bin}()}
}
