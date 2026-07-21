\name{get_peak}
\alias{get_peak}
\title{Find one peak near a target RT}
\description{
Finds the best peak within a given retention-time window around one target RT.
Use this when you already know approximately where to look.
}
\usage{
get_peak(
  x,
  y,
  rt,
  range,
  min_integral = NaN,
  min_intensity = .DEFAULTS$min_intensity,
  min_peak_width_points = .DEFAULTS$min_peak_width_points,
  noise = NaN,
  auto_noise = .DEFAULTS$auto_noise,
  auto_baseline = .DEFAULTS$auto_baseline,
  lambda = 0L,
  max_iterations = 0L,
  allow_overlap = FALSE,
  min_snr = .DEFAULTS$min_snr,
  min_r2 = .DEFAULTS$min_r2,
  shape = .DEFAULTS$shape,
  kernel_size = .DEFAULTS$kernel_size
)
}
\arguments{
  \item{x}{Numeric retention-time vector.}
  \item{y}{Numeric intensity vector, same length as \code{x}. Minimum 3 points.}
  \item{rt}{Target retention time.}
  \item{range}{Half-width of the search window around \code{rt}.}
  \item{min_integral}{Minimum peak area. \code{NaN} disables this filter.}
  \item{min_intensity}{Minimum peak apex intensity. Default 1000.}
  \item{min_peak_width_points}{Minimum peak width in number of data points. Default 5.}
  \item{noise}{Fixed noise level. \code{NaN} disables. Use with \code{auto_noise = FALSE}.}
  \item{auto_noise}{If \code{TRUE}, estimate noise level from the signal. Default \code{TRUE}.}
  \item{auto_baseline}{If \code{TRUE}, subtract a baseline before detection. Default \code{TRUE}.}
  \item{lambda}{airPLS smoothing parameter for baseline estimation. \code{0} uses a default.}
  \item{max_iterations}{Maximum airPLS iterations. \code{0} uses a default.}
  \item{allow_overlap}{If \code{TRUE}, keep peaks whose boundaries overlap. Default \code{FALSE}.}
  \item{min_snr}{Minimum signal-to-noise ratio. Default 1.}
  \item{min_r2}{Minimum fit R-squared (Gaussian or EMG) for the artifact filter; \code{0} disables it (keep peaks regardless of shape). Default 0.3.}
  \item{shape}{Peak shape used when fitting.}

  \item{kernel_size}{Smoothing kernel size used before peak detection.}

}
\value{
A list describing the detected peak, or an empty result if no peak was found.
Fields: \code{rt}, \code{from}, \code{to}, \code{intensity}, \code{integral},
\code{n_points}.
}
\examples{
\donttest{
x <- seq(0, 10, length.out = 200)
y <- dnorm(x, mean = 4, sd = 0.2) * 10000
get_peak(x, y, rt = 4, range = 0.6, min_intensity = 1000)
}
}
\seealso{
\code{\link{find_peaks}()}, \code{\link{calculate_eic}()}
}
