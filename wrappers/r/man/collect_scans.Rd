\name{collect_scans}
\alias{collect_scans}
\title{Collect scans in a retention-time range}
\description{
Gets scans from a loaded file in one time range.
}
\usage{
collect_scans(bin, from, to, level = 1L, include_metadata = FALSE)
}
\arguments{
  \item{bin}{A file handle returned by \code{\link{parse_mzml}()} or \code{\link{parse_bin}()}.}
  \item{from}{Start retention time.}
  \item{to}{End retention time.}
  \item{level}{MS level to collect.}
  \item{include_metadata}{Logical flag to include scan metadata.}
}
\value{
A parsed R object describing the collected scans.
}
\details{
\code{include_metadata = TRUE} returns more information, but it also makes the result larger.
}
\seealso{
\code{\link{calculate_eic}()}, \code{\link{parse_mzml}()}, \code{\link{parse_bin}()}
}
