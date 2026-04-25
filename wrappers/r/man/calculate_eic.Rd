\name{calculate_eic}
\alias{calculate_eic}
\title{Calculate an extracted ion chromatogram}
\description{
Gets an extracted ion chromatogram for one target m/z in one time range.
}
\usage{
calculate_eic(bin, targets, from, to, ppm_tolerance = 20, mz_tolerance = 0.005)
}
\arguments{
  \item{bin}{A file handle returned by \code{\link{parse_mzml}()} or \code{\link{parse_bin}()}.}
  \item{targets}{One target m/z value.}
  \item{from}{Start retention time.}
  \item{to}{End retention time.}
  \item{ppm_tolerance}{Mass tolerance in ppm.}
  \item{mz_tolerance}{Absolute mass tolerance in Da.}
}
\value{
A list with chromatogram x and y values.
}
\details{
The result is a list with \code{x} for time and \code{y} for intensity.
}
\examples{
\dontrun{
eic <- calculate_eic(file, 288.1326, 2.5, 4.0, 10, 0.005)
plot(eic$x, eic$y, type = "l")
}
}
\seealso{
\code{\link{get_scans}()}, \code{\link{get_peak}()}, \code{\link{find_peaks}()}
}
