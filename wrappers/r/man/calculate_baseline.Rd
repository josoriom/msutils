\name{calculate_baseline}
\alias{calculate_baseline}
\title{Estimate a signal baseline}
\description{
Fits a rolling baseline to one numeric signal vector.
}
\usage{
calculate_baseline(y, baseline_window = 15L, baseline_window_factor = 1L)
}
\arguments{
  \item{y}{Numeric signal vector.}
  \item{baseline_window}{Rolling window size used to estimate the baseline.}
  \item{baseline_window_factor}{Scales the effective window. Larger values give a smoother baseline.}
}
\value{
Numeric vector, same length as \code{y}, with the estimated baseline values.
}
\examples{
\donttest{
x <- seq(0, 10, length.out = 1000)
y <- sin(x) + 0.1 * x + 2
bl <- calculate_baseline(y, baseline_window = 101L, baseline_window_factor = 3L)
plot(x, y, type = "l")
lines(x, bl, col = 2)
}
}
\seealso{
\code{\link{find_peaks}()}, \code{\link{get_peak}()}
}
