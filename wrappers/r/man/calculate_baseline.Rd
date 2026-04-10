\name{calculate_baseline}
\alias{calculate_baseline}
\title{Calculate a baseline}
\description{
Estimates a baseline for one numeric signal.
}
\usage{
calculate_baseline(y, baseline_window = 15L, baseline_window_factor = 1L)
}
\arguments{
  \item{y}{A numeric signal vector.}
  \item{baseline_window}{Window size used to estimate the baseline.}
  \item{baseline_window_factor}{Extra scaling factor for the baseline window.}
}
\value{
A numeric vector with the estimated baseline.
}
\examples{
xs <- seq(0, 10, length.out = 1000)
ys <- sin(xs) + 0.1 * xs + 2
baseline <- calculate_baseline(ys, baseline_window = 101, baseline_window_factor = 3)
plot(xs, ys, type = "l")
lines(xs, baseline, col = 2)
}
\seealso{
\code{\link{find_peaks}()}, \code{\link{get_peak}()}
}
