suppressMessages(library(quantion))

MASS <- 90.05550
FROM_RT <- 0.0
TO_RT <- 5.0
EIC_PPM <- 20.0
EIC_MZ <- 0.005

same_bits <- function(left, right) {
  if (length(left) != length(right)) return(FALSE)
  all(mapply(function(a, b) {
    identical(writeBin(a, raw(), size = 8), writeBin(b, raw(), size = 8))
  }, left, right))
}

.origin <- function(url) sub("^(https?://[^/]+).*$", "\\1", url)

.stats <- function(url, path) {
  text <- paste(readLines(paste0(.origin(url), path), warn = FALSE), collapse = "")
  numbers <- as.numeric(regmatches(text, gregexpr("[0-9]+", text))[[1]])
  list(bytes_sent = numbers[1], requests = numbers[2])
}

main <- function() {
  url <- Sys.getenv("QUANTION_REMOTE_URL", "")
  if (!nzchar(url)) {
    message("QUANTION_REMOTE_URL is not set")
    return(1L)
  }
  fixture <- "core/tests/fixtures/api/api.ion"

  invisible(.stats(url, "/reset"))
  remote <- parse_ion_remote(url)
  opening <- .stats(url, "/stats")

  invisible(.stats(url, "/reset"))
  local <- parse_ion_path(fixture)

  eic_remote <- calculate_eic(remote, MASS, FROM_RT, TO_RT, EIC_PPM, EIC_MZ)
  eic_local <- calculate_eic(local, MASS, FROM_RT, TO_RT, EIC_PPM, EIC_MZ)
  one_query <- .stats(url, "/stats")

  total <- file.info(fixture)$size
  same <- same_bits(eic_remote$x, eic_local$x) && same_bits(eic_remote$y, eic_local$y)

  cat(sprintf("file_bytes = %.0f\n", total))
  cat(sprintf("opening_bytes = %.0f\n", opening$bytes_sent))
  cat(sprintf("opening_requests = %.0f\n", opening$requests))
  cat(sprintf("query_bytes = %.0f\n", one_query$bytes_sent))
  cat(sprintf("query_requests = %.0f\n", one_query$requests))
  cat(sprintf("matches_local = %s\n", if (same) "yes" else "no"))
  0L
}

quit(status = main())
