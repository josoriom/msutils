\name{find_peaks}
\alias{find_peaks}
\title{Find peaks in a 1D trace}
\description{
Finds peaks in one numeric trace, where \code{x} is time and \code{y} is intensity.
}
\usage{
find_peaks(
  x, y,
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
  \item{x}{Numeric x values, usually retention time.}
  \item{y}{Numeric y values, usually intensity.}
  \item{integral_threshold}{Minimum peak area threshold.}
  \item{intensity_threshold}{Minimum intensity threshold.}
  \item{width_threshold}{Minimum peak width.}
  \item{noise}{Manual noise estimate.}
  \item{auto_noise}{Logical flag to estimate noise automatically.}
  \item{auto_baseline}{Logical flag to estimate baseline automatically.}
  \item{baseline_window}{Baseline window size.}
  \item{baseline_window_factor}{Baseline window factor.}
  \item{allow_overlap}{Logical flag to allow overlapping peaks.}
  \item{window_size}{Optional local window size.}
  \item{sn_ratio}{Signal-to-noise ratio threshold.}
}
\value{
A parsed R object describing detected peaks, usually a data frame.
}
\details{
This works on numeric vectors. It does not take a file handle.
}
\examples{
x <- seq(0, 10, length.out = 200)
y <- dnorm(x, mean = 4, sd = 0.2) * 10000
find_peaks(x, y, intensity_threshold = 1000, width_threshold = 5L)
}
\seealso{
\code{\link{get_peak}()}, \code{\link{calculate_baseline}()}
}
