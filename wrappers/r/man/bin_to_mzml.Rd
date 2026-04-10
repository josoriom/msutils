\name{bin_to_mzml}
\alias{bin_to_mzml}
\title{Convert a parsed handle to mzML text}
\description{
Converts a loaded file into mzML text.
}
\usage{
bin_to_mzml(bin)
}
\arguments{
  \item{bin}{A file handle returned by \code{\link{parse_mzml}()} or \code{\link{parse_bin}()}.}
}
\value{
A character string containing mzML text.
}
\details{
If you pass a binary handle, the full mzML content may be rebuilt in memory first.
}
\seealso{
\code{\link{parse_bin}()}, \code{\link{bin_to_json}()}, \code{\link{bin_to_df}()}
}
