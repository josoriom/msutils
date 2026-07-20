\name{parse_ion}
\alias{parse_ion}
\title{Load an ion file}
\description{
Loads an ion file and returns a sample file. The input can be a file path,
which is read lazily from disk, or a raw vector already held in memory.
A string starting with \code{http://} or \code{https://} is read over the network
with \code{parse_ion_remote()}.
}
\usage{
parse_ion(input, max_cache_size = 0)
}
\arguments{
  \item{input}{A file path to an ion file, or a raw vector of ion bytes.}
  \item{max_cache_size}{Cache size in bytes. \code{0} sets no limit. Larger values speed up repeated reads at the cost of RAM.}
}
\value{
A Sample file for use in other \code{quantion} functions.
}
\examples{
\donttest{
sample <- parse_ion("/path/to/file.ion")

buf <- readBin("/path/to/file.ion", "raw", file.info("/path/to/file.ion")$size)
sample <- parse_ion(buf)
}
}
\seealso{
\code{\link{parse_mzml}()}, \code{\link{mzml_to_ion}()}, \code{\link{parse_ion_remote}()}, \code{\link{parse_ion_path}()}, \code{\link{dispose}()}
}
