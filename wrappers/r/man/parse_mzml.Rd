\name{parse_mzml}
\alias{parse_mzml}
\title{Parse mzML bytes}
\description{
Reads raw mzML bytes and returns a file handle.
}
\usage{
parse_mzml(data)
}
\arguments{
  \item{data}{A raw vector with the bytes from one mzML file.}
}
\value{
A file handle for use in other \code{msutils} functions.
}
\details{
This does not return a data frame. It returns a handle to the loaded file.
}
\examples{
\dontrun{
mzml_path <- "/path/to/file.mzML"
buf <- readBin(mzml_path, "raw", file.info(mzml_path)$size)
handle <- parse_mzml(buf)
}
}
\seealso{
\code{\link{parse_bin}()}, \code{\link{mzml_to_bin}()}, \code{\link{dispose}()}
}
