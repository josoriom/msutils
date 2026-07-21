\name{find_noise_level}
\alias{find_noise_level}
\title{Estimate the noise level of a signal}
\description{
Estimates the noise level of an intensity vector. The value is the intensity
below which a point is treated as noise, and it is what \code{auto_noise = TRUE}
computes for you in the peak finders.
}
\usage{
find_noise_level(
  y
)
}
\arguments{
  \item{y}{A numeric vector of intensities.}
}
\value{
A list with \code{width}, the number of points in the window the estimate came
from, and \code{intensity}, the estimated noise level.
}
\examples{
\donttest{
eic <- calculate_eic(file, 89.04768, from = 0, to = 5)
noise <- find_noise_level(eic$y)$intensity
}
}
\seealso{
\code{\link{calculate_baseline}()}, \code{\link{find_peaks}()}, \code{\link{get_peak}()}
}
