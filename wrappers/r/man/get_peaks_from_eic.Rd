\name{get_peaks_from_eic}
\alias{get_peaks_from_eic}
\title{Batch targeted peak detection from EICs}
\description{
Finds many target peaks by first building EICs and then searching near the given RT values.
}
\usage{
get_peaks_from_eic(
  bin, df, from = 0.5, to = 5, cores = 1L,
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
  \item{bin}{A file handle returned by \code{\link{parse_mzml}()} or \code{\link{parse_bin}()}.}
  \item{df}{A data frame with columns \code{id}, \code{rt}, \code{mz}, and \code{ranges}.}
  \item{from}{Start retention time for the EIC extraction window.}
  \item{to}{End retention time for the EIC extraction window.}
  \item{cores}{Number of CPU cores to use.}
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
A data frame of targeted peaks.
}
\details{
The input table must have the columns \code{id}, \code{rt}, \code{mz}, and \code{ranges}.
}
\examples{
\dontrun{
targets <- data.frame(
  id = c("Leucine", "Alanine"),
  rt = c(4.15, 2.39),
  mz = c(302.1499, 260.1030),
  ranges = c(0.2, 0.2)
)
get_peaks_from_eic(file, targets, cores = 2L)
}
}
\seealso{
\code{\link{calculate_eic}()}, \code{\link{get_peak}()}, \code{\link{get_peaks_from_chrom}()}
}
