find_sample_root <- function() {
  here <- normalizePath(getwd(), mustWork = FALSE)
  repeat {
    if (file.exists(file.path(here, "core", "tests", "fixtures", "api", "api.ion")))
      return(here)
    parent <- dirname(here)
    if (identical(parent, here)) return("")
    here <- parent
  }
}

sample_ion_path <- function() {
  root <- find_sample_root()
  if (!nzchar(root)) return("")
  file.path(root, "core", "tests", "fixtures", "api", "api.ion")
}

write_sample_bytes <- function(size) {
  path <- tempfile(fileext = ".bin")
  writeBin(as.raw(rep(1:255, length.out = size)), path)
  path
}

file_url <- function(path) paste0("file://", path)
