\name{find_peaks}
\alias{find_peaks}
\title{Find all peaks in a trace}
\description{
Detects peaks in one numeric trace (\code{x} = time, \code{y} = intensity).
Returns all peaks that pass the given filters.
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
  \item{x}{Numeric retention-time vector.}
  \item{y}{Numeric intensity vector, same length as \code{x}. Minimum 3 points.}
  \item{integral_threshold}{Minimum peak area. \code{NaN} disables this filter.}
  \item{intensity_threshold}{Minimum peak apex intensity. Default 1000.}
  \item{width_threshold}{Minimum peak width in number of data points. Default 5.}
  \item{noise}{Fixed noise level. \code{NaN} disables. Use with \code{auto_noise = FALSE}.}
  \item{auto_noise}{If \code{TRUE}, estimate noise level from the signal. Default \code{TRUE}.}
  \item{auto_baseline}{If \code{TRUE}, subtract a rolling baseline before detection. Default \code{TRUE}.}
  \item{baseline_window}{Baseline rolling window size. \code{0} uses a default.}
  \item{baseline_window_factor}{Scales the baseline window. \code{0} uses a default.}
  \item{allow_overlap}{If \code{TRUE}, keep peaks whose boundaries overlap. Default \code{FALSE}.}
  \item{window_size}{Savitzky-Golay smoothing window. \code{0} selects automatically.}
  \item{sn_ratio}{Minimum signal-to-noise ratio. Default 1.}
}
\value{
A data frame with one row per peak. Columns: \code{rt}, \code{from}, \code{to},
\code{intensity}, \code{integral}, \code{np} (number of points), \code{noise}.
}
\examples{
\donttest{
x <- seq(0, 10, length.out = 200)
y <- dnorm(x, mean = 4, sd = 0.2) * 10000
find_peaks(x, y, intensity_threshold = 1000, width_threshold = 5L)
}
}
\seealso{
\code{\link{get_peak}()}, \code{\link{calculate_baseline}()}
}
