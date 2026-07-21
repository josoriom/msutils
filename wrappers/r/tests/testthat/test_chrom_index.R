find_repo_root <- function() {
  here <- normalizePath(getwd(), mustWork = FALSE)
  repeat {
    if (file.exists(file.path(here, "core", "tests", "fixtures", "api", "api.ion")))
      return(here)
    parent <- dirname(here)
    if (identical(parent, here)) return("")
    here <- parent
  }
}

repo_root <- find_repo_root()
sample_file <- if (nzchar(repo_root))
  file.path(repo_root, "core", "tests", "fixtures", "api", "api.ion") else ""

run_apart <- function(work, ...) {
  callr::r(work, args = list(...),
           env = c(QUANTION_ARTIFACTS_ROOT = file.path(repo_root, "artifacts")))
}

ask_get_peaks_from_chrom <- function(path, idx) {
  library(quantion)
  sample <- parse_ion_path(path)
  get_peaks_from_chrom(sample, data.frame(id = "a", idx = idx, rt = 1, range = 0.1))
}

test_that("get_peaks_from_chrom refuses a negative chromatogram index", {
  skip_if_not_installed("quantion")
  skip_if_not_installed("callr")
  skip_if(!nzchar(sample_file), "no api.ion sample next to the package")

  expect_error(run_apart(ask_get_peaks_from_chrom, sample_file, -1L), "non-negative")
})

test_that("get_peaks_from_chrom refuses an index that is not a number", {
  skip_if_not_installed("quantion")
  skip_if_not_installed("callr")
  skip_if(!nzchar(sample_file), "no api.ion sample next to the package")

  expect_error(run_apart(ask_get_peaks_from_chrom, sample_file, "no"), "non-negative")
})
