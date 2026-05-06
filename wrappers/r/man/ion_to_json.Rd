\name{ion_to_json}
\alias{ion_to_json}
\title{Get JSON from a sample}
\description{
Returns the full file content as a JSON string.
}
\usage{
ion_to_json(file)
}
\arguments{
  \item{file}{Sample file from \code{\link{parse_mzml}()} or \code{\link{parse_ion}()}.}
}
\value{
A JSON character string with all run metadata, spectra, and chromatograms.
}
\seealso{
\code{\link{ion_to_df}()}, \code{\link{ion_to_mzml}()}, \code{\link{parse_ion}()}
}
