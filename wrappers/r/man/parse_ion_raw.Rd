\name{parse_ion_raw}
\alias{parse_ion_raw}
\title{Load an ion file already held in memory}
\description{
Loads a raw vector of ion bytes and returns a sample file. Use it when the file
is already in memory and there is nothing to read from disk.
\code{\link{parse_ion}()} calls this for you when its input is a raw vector.
}
\usage{
parse_ion_raw(
  input,
  max_cache_size = 0
)
}
\arguments{
  \item{input}{A raw vector of ion bytes.}
  \item{max_cache_size}{Cache size in bytes. \code{0} sets no limit. Larger values speed up repeated reads at the cost of RAM.}
}
\value{
A Sample file for use in other \code{quantion} functions.
}
\examples{
\donttest{
buf <- readBin("/path/to/file.ion", "raw", file.info("/path/to/file.ion")$size)
sample <- parse_ion_raw(buf)
}
}
\seealso{
\code{\link{parse_ion}()}, \code{\link{parse_ion_path}()}, \code{\link{parse_ion_remote}()}, \code{\link{dispose}()}
}
