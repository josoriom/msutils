\name{find_features}
\alias{find_features}
\title{Untargeted feature detection}
\description{
Finds untargeted features in one file across m/z and retention time.
}
\usage{
find_features(
  data, from = 0, to = 10,
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
  \item{data}{A file handle returned by \code{\link{parse_mzml}()} or \code{\link{parse_bin}()}.}
  \item{from}{Start retention time.}
  \item{to}{End retention time.}
  \item{ppm_tolerance}{PPM tolerance for EIC extraction.}
  \item{mz_tolerance}{Absolute m/z tolerance for EIC extraction.}
  \item{grid_start}{Start m/z for the search grid.}
  \item{grid_end}{End m/z for the search grid.}
  \item{grid_step_ppm}{Grid step size.}
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
A data frame of detected untargeted features.
}
\details{
\code{use_gpu = TRUE} only helps when GPU support is available.
}
\examples{
\dontrun{
features <- find_features(
  file,
  from = 0,
  to = 10,
  cores = 5L,
  auto_noise = TRUE,
  auto_baseline = TRUE,
  width_threshold = 10L,
  intensity_threshold = 1000
)
}
}
\seealso{
\code{\link{find_feature}()}, \code{\link{get_features}()}
}
