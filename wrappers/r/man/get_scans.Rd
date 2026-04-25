\name{get_scans}
\alias{get_scans}
\title{Retrieve scans by query}
\description{
Fetches scans from a loaded file using one of four query modes.
Exactly one query mode must be provided per call.
}
\usage{
get_scans(bin, rt_from = NULL, rt_to = NULL, rt = NULL,
          mz_from = NULL, mz_to = NULL, mz = NULL,
          level = 1L)
}
\arguments{
  \item{bin}{A file handle returned by \code{\link{parse_mzml}()} or \code{\link{parse_bin}()}.}
  \item{rt_from}{Start of retention-time window in minutes (use with \code{rt_to}).}
  \item{rt_to}{End of retention-time window in minutes (use with \code{rt_from}).}
  \item{rt}{Target retention time in minutes — returns the single closest scan.}
  \item{mz_from}{Lower bound of \code{selected_ion_mz} window (use with \code{mz_to}).}
  \item{mz_to}{Upper bound of \code{selected_ion_mz} window (use with \code{mz_from}).}
  \item{mz}{Target \code{selected_ion_mz} — returns the single closest scan.}
  \item{level}{MS level to query (integer in 0-255). Defaults to 1.}
}
\details{
Exactly one of the four query modes must be provided:
\itemize{
  \item \code{rt_from} + \code{rt_to}: all scans whose RT falls in \code{[rt_from, rt_to]}.
  \item \code{rt}: the single scan whose RT is closest to the target.
  \item \code{mz_from} + \code{mz_to}: all scans whose \code{selected_ion_mz} falls in the window.
  \item \code{mz}: the single scan whose \code{selected_ion_mz} is closest to the target.
}
}
\value{
For range queries (\code{rt_from}/\code{rt_to} or \code{mz_from}/\code{mz_to}): a list of scans,
each containing \code{rt}, \code{mz}, \code{intensity}, and \code{metadata}.

For point queries (\code{rt} or \code{mz}): a single scan as a list with \code{rt}, \code{mz},
\code{intensity}, and \code{metadata}, or \code{NULL} if no scan was found.
}
\examples{
\dontrun{
# All MS1 scans between 1.0 and 2.5 minutes
scans <- get_scans(bin, rt_from = 1.0, rt_to = 2.5, level = 1L)

# MS1 scan closest to 1.8 minutes
scan  <- get_scans(bin, rt = 1.8, level = 1L)

# MS2 scans with selected_ion_mz between 400 and 450 Da
scans <- get_scans(bin, mz_from = 400.0, mz_to = 450.0, level = 2L)

# MS2 scan whose precursor is closest to 445.12 Da
scan  <- get_scans(bin, mz = 445.12, level = 2L)
}
}
\seealso{
\code{\link{calculate_eic}()}, \code{\link{parse_mzml}()}, \code{\link{parse_bin}()}
}
