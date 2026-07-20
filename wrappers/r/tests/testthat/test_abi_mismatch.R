build_wrong_version_library <- function() {
  folder <- tempfile("wrong_abi")
  dir.create(folder)
  source_file <- file.path(folder, "wrong_abi.c")
  writeLines(
    c("#include <stdint.h>",
      "#include <stddef.h>",
      "uint32_t quantion_abi_version(void) { return 999; }",
      "size_t quantion_sizeof_peak_options(void) { return 0; }"),
    source_file
  )
  built <- suppressWarnings(system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "SHLIB", shQuote(source_file)),
    stdout = FALSE, stderr = FALSE
  ))
  made <- file.path(folder, paste0("wrong_abi", .Platform$dynlib.ext))
  if (built != 0 || !file.exists(made)) "" else made
}

load_the_package <- function() {
  library(quantion)
  TRUE
}

test_that("a native library with the wrong ABI version says so", {
  skip_if_not_installed("quantion")
  skip_if_not_installed("callr")

  wrong_library <- build_wrong_version_library()
  skip_if(!nzchar(wrong_library), "no C compiler to build a wrong-version library")

  expect_error(
    callr::r(load_the_package, env = c(QUANTION_LIB = wrong_library)),
    "ABI version mismatch"
  )
})
