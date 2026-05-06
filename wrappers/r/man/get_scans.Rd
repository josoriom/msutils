\name{get_scans}
\alias{get_scans}
\title{Get scans by RT or m/z query}
\description{
Fetches scans from a loaded file. Exactly one query mode must be provided per call.
}
\usage{
get_scans(file, rt_from = NULL, rt_to = NULL, rt = NULL,
          mz_from = NULL, mz_to = NULL, mz = NULL,
          level = 1L)
}
\arguments{
  \item{file}{Sample file from \code{\link{parse_mzml}()} or \code{\link{parse_ion}()}.}
  \item{rt_from}{RT range start (minutes). Use with \code{rt_to}.}
  \item{rt_to}{RT range end (minutes). Use with \code{rt_from}.}
  \item{rt}{Target RT — returns the single closest scan.}
  \item{mz_from}{m/z range lower bound (Da). Use with \code{mz_to}.}
  \item{mz_to}{m/z range upper bound (Da). Use with \code{mz_from}.}
  \item{mz}{Target precursor m/z — returns the single closest scan.}
  \item{level}{MS level to query (integer 0–255). Default 1.}
}
\value{
Range query (\code{rt_from}/\code{rt_to} or \code{mz_from}/\code{mz_to}): a list of scans,
each with \code{rt}, \code{mz}, \code{intensity}, and metadata.

Point query (\code{rt} or \code{mz}): a single scan list, or \code{NULL} if none found.
}
\examples{
\donttest{
# All MS1 scans from 1.0 to 2.5 min
scans <- get_scans(sample, rt_from = 1.0, rt_to = 2.5)

# Closest MS1 scan to 1.8 min
scan  <- get_scans(sample, rt = 1.8)

# MS2 scans with precursor m/z in 400–450 Da
scans <- get_scans(sample, mz_from = 400.0, mz_to = 450.0, level = 2L)

# Closest MS2 scan to precursor m/z 445.12 Da
scan  <- get_scans(sample, mz = 445.12, level = 2L)
}
}
\seealso{
\code{\link{calculate_eic}()}, \code{\link{parse_mzml}()}, \code{\link{parse_ion}()}
}
