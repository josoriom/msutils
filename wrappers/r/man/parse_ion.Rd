\name{parse_ion}
\alias{parse_ion}
\title{Load an ion file}
\description{
Reads raw ion binary bytes and returns a sample file.
}
\usage{
parse_ion(file, max_cache_size = 0)
}
\arguments{
  \item{file}{Raw vector — bytes from one ion file.}
  \item{max_cache_size}{Cache size in bytes. \code{0} sets no limit. Larger values speed up repeated reads at the cost of RAM.}
}
\value{
A Sample file for use in other \code{msutils} functions.
}
\examples{
\donttest{
buf <- readBin("/path/to/file.ion", "raw", file.info("/path/to/file.ion")$size)
sample <- parse_ion(buf)
}
}
\seealso{
\code{\link{parse_mzml}()}, \code{\link{mzml_to_ion}()}, \code{\link{dispose}()}
}
