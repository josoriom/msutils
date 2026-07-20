make_artifacts_root <- function(versions, with_manifest) {
  arch <- quantion:::.quantion_find_arch_dir()
  name <- quantion:::.quantion_rust_basename()
  root <- tempfile("artifacts")
  for (version in versions) {
    dir.create(file.path(root, version, arch), recursive = TRUE, showWarnings = FALSE)
    file.create(file.path(root, version, arch, name))
    if (with_manifest) file.create(file.path(root, version, "manifest.json"))
  }
  root
}

library_in <- function(root, version) {
  file.path(root, version,
            quantion:::.quantion_find_arch_dir(),
            quantion:::.quantion_rust_basename())
}

resolve_under <- function(root, own_version) {
  before <- Sys.getenv("QUANTION_ARTIFACTS_ROOT", "")
  on.exit(Sys.setenv(QUANTION_ARTIFACTS_ROOT = before), add = TRUE)
  Sys.setenv(QUANTION_ARTIFACTS_ROOT = root)
  quantion:::.quantion_find_library(own_version)
}

test_that("the version scan reads versions as numbers, not as strings", {
  skip_if_not_installed("quantion")

  root <- make_artifacts_root(c("0.1.0", "0.9.0", "0.10.0"), with_manifest = TRUE)
  expect_identical(resolve_under(root, ""), library_in(root, "0.10.0"))
})

test_that("a directory with no manifest is not a released version", {
  skip_if_not_installed("quantion")

  root <- make_artifacts_root(c("0.9.0", "0.10.0"), with_manifest = FALSE)
  expect_identical(resolve_under(root, "2.0.0"), "")
})

test_that("the package's own version is taken without a manifest", {
  skip_if_not_installed("quantion")

  root <- make_artifacts_root("0.1.0", with_manifest = FALSE)
  expect_identical(resolve_under(root, "0.1.0"), library_in(root, "0.1.0"))
})

test_that("QUANTION_LIB wins over every artifacts folder", {
  skip_if_not_installed("quantion")

  before <- Sys.getenv("QUANTION_LIB", "")
  on.exit(Sys.setenv(QUANTION_LIB = before), add = TRUE)
  Sys.setenv(QUANTION_LIB = "/chosen/by/hand.so")

  root <- make_artifacts_root("0.10.0", with_manifest = TRUE)
  expect_identical(resolve_under(root, "0.1.0"), "/chosen/by/hand.so")
})
