sample_file <- sample_ion_path()

run_apart <- function(work, ...) {
  callr::r(work, args = list(...),
           env = c(QUANTION_ARTIFACTS_ROOT = file.path(find_sample_root(), "artifacts")))
}

ask_short_mz <- function(path) {
  library(quantion)
  sample <- parse_ion_path(path)
  find_feature(sample, rt = c(1, 2), mz = 89.04768, window = c(0.1, 0.1), id = c("a", "b"))
}

ask_short_window <- function(path) {
  library(quantion)
  sample <- parse_ion_path(path)
  find_feature(sample, rt = c(1, 2), mz = c(89.04768, 89.04768), window = 0.1, id = c("a", "b"))
}

ask_long_window <- function(path) {
  library(quantion)
  sample <- parse_ion_path(path)
  find_feature(sample, rt = 1, mz = 89.04768, window = c(0.1, 0.2), id = "a")
}

test_that("find_feature refuses an mz vector shorter than rt", {
  skip_if_not_installed("quantion")
  skip_if_not_installed("callr")
  skip_if(!nzchar(sample_file), "no api.ion sample next to the package")

  expect_error(run_apart(ask_short_mz, sample_file), "same length")
})

test_that("find_feature refuses a window vector shorter than rt", {
  skip_if_not_installed("quantion")
  skip_if_not_installed("callr")
  skip_if(!nzchar(sample_file), "no api.ion sample next to the package")

  expect_error(run_apart(ask_short_window, sample_file), "same length")
})

test_that("find_feature refuses a window vector longer than rt", {
  skip_if_not_installed("quantion")
  skip_if_not_installed("callr")
  skip_if(!nzchar(sample_file), "no api.ion sample next to the package")

  expect_error(run_apart(ask_long_window, sample_file), "same length")
})
