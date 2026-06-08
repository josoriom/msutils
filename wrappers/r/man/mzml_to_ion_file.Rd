\name{mzml_to_ion_file}
\alias{mzml_to_ion_file}
\title{Stream-convert an mzML file into an ion file}
\description{
Reads an mzML file on disk one scan at a time and writes an ion file directly,
so the whole file is never held in memory. Use this for large files.
}
\usage{
mzml_to_ion_file(input_path, output_path, level = 12L, f32_compress = FALSE,
                 section_on_disk = FALSE)
}
\arguments{
  \item{input_path}{Path to the source mzML file.}
  \item{output_path}{Path where the ion file is written.}
  \item{level}{Compression level, integer 0 (none) to 22 (max). Default 12.}
  \item{f32_compress}{If \code{TRUE}, store intensity values as 32-bit float. Default \code{FALSE}.}
  \item{section_on_disk}{If \code{TRUE}, keep intermediate tables on disk to cap memory. Default \code{FALSE}.}
}
\value{
Invisibly returns \code{output_path}.
}
\examples{
\donttest{
mzml_to_ion_file("/path/to/file.mzML", "/path/to/file.ion", level = 12L)
}
}
\seealso{
\code{\link{mzml_to_ion}()}, \code{\link{parse_ion}()}, \code{\link{ion_to_mzml}()}
}
