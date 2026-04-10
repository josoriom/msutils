\name{dispose}
\alias{dispose}
\title{Dispose a native handle}
\description{
Releases a file handle created by \code{msutils}.
}
\usage{
dispose(bin)
}
\arguments{
  \item{bin}{A file handle returned by \code{\link{parse_mzml}()} or \code{\link{parse_bin}()}.}
}
\value{
Called for its side effect of freeing native memory.
}
\details{
You usually do not need this. R will free the handle later. Use this when you want to free memory now.
}
\seealso{
\code{\link{parse_mzml}()}, \code{\link{parse_bin}()}
}
