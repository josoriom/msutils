\name{bin_to_json}
\alias{bin_to_json}
\title{Convert a parsed handle to JSON}
\description{
Converts a loaded file into JSON text.
}
\usage{
bin_to_json(bin)
}
\arguments{
  \item{bin}{A file handle returned by \code{\link{parse_mzml}()} or \code{\link{parse_bin}()}.}
}
\value{
A JSON string.
}
\details{
If you pass a binary handle, the full mzML content may be rebuilt in memory first.
}
\seealso{
\code{\link{bin_to_df}()}, \code{\link{bin_to_mzml}()}, \code{\link{parse_bin}()}
}
