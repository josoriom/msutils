\name{parse_ion_remote}
\alias{parse_ion_remote}
\title{Read an ion file over HTTP}
\description{
Opens an ion file held on a web server and returns a sample file. Only the
byte ranges a call needs are downloaded, so a large file can be read without
copying all of it. The server must support range requests; if it returns the
whole body instead, this stops with an error.
\code{\link{parse_ion}()} calls this for you when its input starts with
\code{http://} or \code{https://}.
}
\usage{
parse_ion_remote(url, max_cache_size = 0)
}
\arguments{
  \item{url}{A single non-empty URL starting with \code{http://} or \code{https://}.}
  \item{max_cache_size}{Cache size in bytes. \code{0} sets no limit. Larger values speed up repeated reads at the cost of RAM.}
}
\value{
A Sample file for use in other \code{quantion} functions.
}
\examples{
\donttest{
sample <- parse_ion_remote("https://example.org/file.ion")
}
}
\seealso{
\code{\link{parse_ion}()}, \code{\link{parse_ion_path}()}, \code{\link{parse_ion_raw}()}, \code{\link{dispose}()}
}
