\name{ion_to_df}
\alias{ion_to_df}
\title{Get R data from a sample}
\description{
Returns the full file content as nested R lists and data frames.
}
\usage{
ion_to_df(
  bin
)
}
\arguments{
  \item{bin}{Sample file from \code{\link{parse_mzml}()} or \code{\link{parse_ion}()}.}
}
\value{
A nested R object. The top-level list contains \code{run}, which holds
\code{spectra} and \code{chromatograms} as data frames with metadata columns and
list-columns for \code{mz} and \code{intensity} arrays.
}
\seealso{
\code{\link{ion_to_json}()}, \code{\link{ion_to_mzml}()}, \code{\link{parse_ion}()}
}
