\name{calculate_eic}
\alias{calculate_eic}
\title{Get an extracted ion chromatogram (EIC)}
\description{
Extracts an EIC for one target m/z over a given retention-time range.
}
\usage{
calculate_eic(
  bin,
  targets,
  from,
  to,
  ppm_tolerance = 20,
  mz_tolerance = 0.005
)
}
\arguments{
  \item{bin}{Sample file from \code{\link{parse_mzml}()} or \code{\link{parse_ion}()}.}
  \item{targets}{Single numeric m/z target.}
  \item{from}{Start retention time (minutes).}
  \item{to}{End retention time (minutes).}
  \item{ppm_tolerance}{Mass tolerance in ppm.}
  \item{mz_tolerance}{Absolute mass tolerance in Da.}
}
\value{
A list with \code{x} (retention time, numeric) and \code{y} (intensity, numeric).
}
\examples{
\donttest{
eic <- calculate_eic(file, 288.1326, from = 2.5, to = 4.0)
}
}
\seealso{
\code{\link{get_scans}()}, \code{\link{find_peaks}()}, \code{\link{get_peak}()}
}
