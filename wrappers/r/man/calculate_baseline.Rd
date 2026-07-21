\name{calculate_baseline}
\alias{calculate_baseline}
\title{Estimate a signal baseline}
\description{
Fits an airPLS baseline to one numeric signal vector.
}
\usage{
calculate_baseline(
  y,
  lambda = 0L,
  max_iterations = 0L
)
}
\arguments{
  \item{y}{Numeric signal vector.}
  \item{lambda}{airPLS smoothing parameter. \code{0} uses a default.}
  \item{max_iterations}{Maximum airPLS iterations. \code{0} uses a default.}
}
\value{
Numeric vector, same length as \code{y}, with the estimated baseline values.
}
\examples{
\donttest{
x <- seq(0, 10, length.out = 1000)
y <- sin(x) + 0.1 * x + 2
bl <- calculate_baseline(y)
plot(x, y, type = "l")
lines(x, bl, col = 2)
}
}
\seealso{
\code{\link{find_peaks}()}, \code{\link{get_peak}()}
}
