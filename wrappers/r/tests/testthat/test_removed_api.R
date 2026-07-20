test_that("get_ion_image is gone from the package", {
  skip_if_not_installed("quantion")

  expect_false("get_ion_image" %in% getNamespaceExports("quantion"))
  expect_false(exists("get_ion_image", envir = asNamespace("quantion"), inherits = FALSE))
})
