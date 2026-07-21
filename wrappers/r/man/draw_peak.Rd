\name{draw_peak}
\alias{draw_peak}
\title{Draw a peak from its shape parameters}
\description{
Evaluates a fitted peak shape at the given retention times, so the model curve
can be plotted next to the measured signal.
}
\usage{
draw_peak(
  x,
  params
)
}
\arguments{
  \item{x}{A numeric vector of retention times to evaluate the shape at.}
  \item{params}{A list of shape parameters, as returned by \code{\link{fit_peak}()}:
    \code{shape} (\code{"gaussian"} or \code{"emg"}), \code{height}, \code{center},
    \code{fwhm} and, for \code{"emg"}, \code{tail}. A missing \code{tail} is read as \code{0}.}
}
\value{
A numeric vector of intensities, the same length as \code{x}.
}
\examples{
\donttest{
eic <- calculate_eic(file, 89.04768, from = 0, to = 5)
peak <- get_peak(eic$x, eic$y, rt = 2.385, range = 0.1)
params <- fit_peak(eic$x, eic$y, peak$rt, peak$intensity)
curve <- draw_peak(eic$x, params)
}
}
\seealso{
\code{\link{fit_peak}()}, \code{\link{get_peak}()}
}
