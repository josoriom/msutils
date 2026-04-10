\name{parse_bin}
\alias{parse_bin}
\title{Parse binary msutils bytes}
\description{
Reads raw binary msutils bytes and returns a file handle.
}
\usage{
parse_bin(bin, max_cache_size = 0)
}
\arguments{
  \item{bin}{A raw vector with the bytes from one binary msutils file.}
  \item{max_cache_size}{Lazy-read cache size in bytes. Bigger values may use more RAM but can speed up repeated reads. \code{0} means no cache limit.}
}
\value{
A file handle for use in other \code{msutils} functions.
}
\details{
Use this for binary files such as \code{.b64}. The result is a handle, not a data frame.
}
\examples{
\dontrun{
bin_path <- "/path/to/file.b64"
buf <- readBin(bin_path, "raw", file.info(bin_path)$size)
handle <- parse_bin(buf)
}
}
\seealso{
\code{\link{parse_mzml}()}, \code{\link{bin_to_mzml}()}, \code{\link{dispose}()}
}
