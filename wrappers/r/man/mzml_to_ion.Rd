\name{mzml_to_ion}
\alias{mzml_to_ion}
\title{Convert a sample to ion binary bytes}
\description{
Encodes a loaded file as compressed ion binary bytes.
}
\usage{
mzml_to_ion(
  bin,
  level = 12L,
  f32_compress = FALSE
)
}
\arguments{
  \item{bin}{Sample file from \code{\link{parse_mzml}()} or \code{\link{parse_ion}()}.}
  \item{level}{Compression level, integer 0 (none) to 22 (max). Default 12.}
  \item{f32_compress}{If \code{TRUE}, compress intensity values to 32-bit float. Default \code{FALSE}.}
}
\value{
A raw vector of ion binary bytes. Write to disk with \code{writeBin()}.
}
\examples{
\donttest{
buf <- readBin("/path/to/file.mzML", "raw", file.info("/path/to/file.mzML")$size)
sample <- parse_mzml(buf)
out <- mzml_to_ion(sample, level = 12L)
writeBin(out, "/path/to/file.ion")
}
}
\seealso{
\code{\link{parse_mzml}()}, \code{\link{parse_ion}()}, \code{\link{ion_to_mzml}()}
}
