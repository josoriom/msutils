\name{bin_to_df}
\alias{bin_to_df}
\title{Convert a parsed handle to R objects}
\description{
Converts a loaded file into nested R lists and data frames.
}
\usage{
bin_to_df(bin)
}
\arguments{
  \item{bin}{A file handle returned by \code{\link{parse_mzml}()} or \code{\link{parse_bin}()}.}
}
\value{
A nested R object with run, spectra, chromatograms, and metadata.
}
\details{
If you pass a binary handle, the full mzML content may be rebuilt in memory first.
}
\seealso{
\code{\link{bin_to_json}()}, \code{\link{bin_to_mzml}()}, \code{\link{parse_bin}()}
}
