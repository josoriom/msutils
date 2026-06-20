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
  \item{x}{Numeric retention-time vector.}
  \item{y}{Numeric intensity vector, same length as \code{x}. Minimum 3 points.}
  \item{min_integral}{Minimum peak area. \code{NaN} disables this filter.}
  \item{min_intensity}{Minimum peak apex intensity. Default 1000.}
  \item{min_peak_width_points}{Minimum peak width in number of data points. Default 5.}
  \item{noise}{Fixed noise level. \code{NaN} disables. Use with \code{auto_noise = FALSE}.}
  \item{auto_noise}{If \code{TRUE}, estimate noise level from the signal. Default \code{TRUE}.}
  \item{auto_baseline}{If \code{TRUE}, subtract a rolling baseline before detection. Default \code{TRUE}.}
  \item{lambda}{airPLS smoothing parameter for baseline estimation. \code{0} uses a default.}
  \item{max_iterations}{Maximum airPLS iterations for baseline estimation. \code{0} uses a default.}
  \item{allow_overlap}{If \code{TRUE}, keep peaks whose boundaries overlap. Default \code{FALSE}.}
  \item{min_snr}{Minimum signal-to-noise ratio. Default 1.}
  \item{min_gaussian_r2}{Minimum Gaussian R-squared for the artifact filter; \code{0} disables it (keep peaks regardless of shape). Default 0.3.}
}
\value{
A data frame with one row per peak. Columns: \code{rt}, \code{from}, \code{to},
\code{intensity}, \code{integral}, \code{n_points} (number of points), \code{noise}.
}
\examples{
\donttest{
x <- seq(0, 10, length.out = 200)
y <- dnorm(x, mean = 4, sd = 0.2) * 10000
find_peaks(x, y, min_intensity = 1000, min_peak_width_points = 5L)
}
}
\seealso{
\code{\link{get_peak}()}, \code{\link{calculate_baseline}()}
}
