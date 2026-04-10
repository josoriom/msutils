\name{mzml_to_bin}
\alias{mzml_to_bin}
\title{Convert a parsed handle to binary bytes}
\description{
Converts a loaded file into binary msutils bytes.
}
\usage{
mzml_to_bin(bin, level = 12L, f32_compress = FALSE)
}
\arguments{
  \item{bin}{A file handle returned by \code{\link{parse_mzml}()} or \code{\link{parse_bin}()}.}
  \item{level}{Compression level from 0 to 22.}
  \item{f32_compress}{Logical flag to allow float32 compression when supported.}
}
\value{
A raw vector containing encoded binary bytes.
}
\details{
You will usually use a handle from \code{\link{parse_mzml}()}. A binary handle also works.
}
\examples{
\dontrun{
mzml_path <- "/path/to/file.mzML"
buf <- readBin(mzml_path, "raw", file.info(mzml_path)$size)
handle <- parse_mzml(buf)
out <- mzml_to_bin(handle, 12L, f32_compress = FALSE)
writeBin(out, "/path/to/file.b64")
}
}
\seealso{
\code{\link{parse_mzml}()}, \code{\link{parse_bin}()}, \code{\link{bin_to_mzml}()}
}
