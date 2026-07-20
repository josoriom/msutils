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

ask_find_feature <- function(path, count) {
  library(quantion)
  sample <- parse_ion_path(path)
  find_feature(
    sample,
    rt = rep(1, count),
    mz = rep(89.04768, count),
    window = rep(0.1, count),
    id = "only-one"
  )
}

ask_bridge_find_feature <- function(path, rt_count, id_count) {
  library(quantion)
  sample <- parse_ion_path(path)
  .Call("C_find_feature", sample,
        as.numeric(rep(1, rt_count)),
        as.numeric(rep(89.04768, rt_count)),
        as.numeric(rep(0.1, rt_count)),
        rep("only-one", id_count),
        10, 0.003, 20, 0.005, NULL, 1L,
        PACKAGE = "quantion")
}

ask_bridge_get_peaks_from_eic <- function(path, rt_count, id_count) {
  library(quantion)
  sample <- parse_ion_path(path)
  .Call("C_get_peaks_from_eic", sample,
        as.numeric(rep(1, rt_count)),
        as.numeric(rep(89.04768, rt_count)),
        as.numeric(rep(0.1, rt_count)),
        rep("only-one", id_count),
        0.5, 0.5, NULL, 1L,
        PACKAGE = "quantion")
}

test_that("find_feature refuses an id vector shorter than rt", {
  skip_if_not_installed("quantion")
  skip_if_not_installed("callr")
  skip_if(!nzchar(sample_file), "no api.ion sample next to the package")

  expect_error(run_apart(ask_find_feature, sample_file, 3), "same length")
})

test_that("the C bridge stops a short id vector before it reads past the end", {
  skip_if_not_installed("quantion")
  skip_if_not_installed("callr")
  skip_if(!nzchar(sample_file), "no api.ion sample next to the package")

  expect_error(run_apart(ask_bridge_find_feature, sample_file, 4e6, 1e5), "same length")
  expect_error(run_apart(ask_bridge_get_peaks_from_eic, sample_file, 4e6, 1e5), "same length")
})
