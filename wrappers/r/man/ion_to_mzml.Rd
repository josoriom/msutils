\name{ion_to_mzml}
\alias{ion_to_mzml}
\title{Get mzML text from a sample}
\description{
Returns the full file content as an mzML character string.
}
\usage{
ion_to_mzml(
  bin
)
}
\arguments{
  \item{bin}{Sample file from \code{\link{parse_mzml}()} or \code{\link{parse_ion}()}.}
}
\value{
A character string with valid mzML content.
}
\seealso{
\code{\link{parse_ion}()}, \code{\link{parse_mzml}()}
}
