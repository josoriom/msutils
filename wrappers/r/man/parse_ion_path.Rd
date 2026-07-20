\name{parse_ion_path}
\alias{parse_ion_path}
\title{Load an ion file from disk}
\description{
Loads an ion file from a path on disk and returns a sample file. The bytes are
read lazily, so only the parts a call needs are touched.
\code{\link{parse_ion}()} calls this for you when its input is a path.
}
\usage{
parse_ion_path(path, max_cache_size = 0)
}
\arguments{
  \item{path}{A single non-empty path to an existing ion file. \code{~} is expanded.}
  \item{max_cache_size}{Cache size in bytes. \code{0} sets no limit. Larger values speed up repeated reads at the cost of RAM.}
}
\value{
A Sample file for use in other \code{quantion} functions.
}
\examples{
\donttest{
sample <- parse_ion_path("/path/to/file.ion")
}
}
\seealso{
\code{\link{parse_ion}()}, \code{\link{parse_ion_raw}()}, \code{\link{parse_ion_remote}()}, \code{\link{dispose}()}
}
