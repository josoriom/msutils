\name{dispose}
\alias{dispose}
\title{Free a sample file}
\description{
Releases native memory held by a sample file. Call this when you want to free
memory immediately rather than waiting for R garbage collection.
}
\usage{
dispose(
  bin
)
}
\arguments{
  \item{bin}{Sample file from \code{\link{parse_mzml}()} or \code{\link{parse_ion}()}.}
}
\value{
\code{NULL}, called for its side effect.
}
\seealso{
\code{\link{parse_mzml}()}, \code{\link{parse_ion}()}
}
