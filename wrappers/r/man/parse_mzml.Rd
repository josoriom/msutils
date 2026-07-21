\name{parse_mzml}
\alias{parse_mzml}
\title{Load an mzML file}
\description{
Reads raw mzML bytes and returns a file handle.
}
\usage{
parse_mzml(
  data
)
}
\arguments{
  \item{data}{Raw vector — bytes from one mzML file.}
}
\value{
A file handle for use in other \code{quantion} functions.
}
\examples{
\donttest{
buf <- readBin("/path/to/file.mzML", "raw", file.info("/path/to/file.mzML")$size)
handle <- parse_mzml(buf)
}
}
\seealso{
\code{\link{parse_ion}()}, \code{\link{mzml_to_ion}()}, \code{\link{dispose}()}
}
