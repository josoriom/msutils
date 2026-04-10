\name{get_peaks_from_chrom}
\alias{get_peaks_from_chrom}
\title{Targeted peak detection from chromatograms}
\description{
Finds many target peaks using chromatogram indexes and time windows.
}
\usage{
get_peaks_from_chrom(
  bin, items, cores = 1L,
  integral_threshold = NaN,
  intensity_threshold = 1000,
  width_threshold = 5L,
  noise = NaN,
  auto_noise = FALSE,
  auto_baseline = FALSE,
  baseline_window = 0L,
  baseline_window_factor = 0L,
  allow_overlap = FALSE,
  window_size = 0L,
  sn_ratio = 1
)
}
\arguments{
  \item{bin}{A file handle returned by \code{\link{parse_mzml}()} or \code{\link{parse_bin}()}.}
  \item{items}{A list or data frame with \code{idx} or \code{index}, and also \code{rt} plus \code{window} or \code{range}.}
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
A data frame of targeted chromatogram peaks.
}
\details{
Use this when you already know which chromatogram index to read.
}
\examples{
\dontrun{
items <- data.frame(
  idx = c(0L, 5L),
  rt = c(1.2, 2.85),
  window = c(0.4, 0.5)
)
get_peaks_from_chrom(file, items, cores = 2L)
}
}
\seealso{
\code{\link{get_peaks_from_eic}()}, \code{\link{get_peak}()}
}
