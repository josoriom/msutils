rows <- function(...) matrix(c(...), ncol = 2, byrow = TRUE)

source_holding <- function(blocks) {
  source <- new.env(parent = emptyenv())
  source$held <- blocks
  source
}

source_over_file <- function(path) {
  source <- new.env(parent = emptyenv())
  source$url <- file_url(path)
  source$store <- .Call("C_store_new", PACKAGE = "quantion")
  source$held <- list()
  source$total <- file.info(path)$size
  source
}

test_that("the gap is the same constant whatever the total size is", {
  skip_if_not_installed("quantion")

  expect_identical(quantion:::.gap_for(0), 131072)
  expect_identical(quantion:::.gap_for(-1), 131072)
  expect_identical(quantion:::.gap_for(80000), 131072)
  expect_identical(quantion:::.gap_for(8e9), 131072)
})

test_that("merging drops the ranges that ask for no bytes", {
  skip_if_not_installed("quantion")

  expect_equal(quantion:::.merge_ranges(rows(0, 0, 10, 5), gap = 0), rows(10, 5))
})

test_that("merging joins two ranges that sit closer than the gap", {
  skip_if_not_installed("quantion")

  expect_equal(quantion:::.merge_ranges(rows(0, 10, 100, 10), gap = 100), rows(0, 110))
})

test_that("merging keeps two ranges that sit further apart than the gap", {
  skip_if_not_installed("quantion")

  expect_equal(
    quantion:::.merge_ranges(rows(0, 10, 5000, 10), gap = 100),
    rows(0, 10, 5000, 10)
  )
})

test_that("merging swallows a range that sits inside another", {
  skip_if_not_installed("quantion")

  expect_equal(quantion:::.merge_ranges(rows(0, 100, 10, 20), gap = 0), rows(0, 100))
})

test_that("a span counts as held only when one block covers all of it", {
  skip_if_not_installed("quantion")

  source <- source_holding(list(c(0, 1024)))

  expect_true(quantion:::.already_held(source, 0, 1024))
  expect_true(quantion:::.already_held(source, 100, 924))
  expect_false(quantion:::.already_held(source, 100, 925))
  expect_false(quantion:::.already_held(source, 2000, 10))
})

test_that("a span spread over two touching blocks is not held", {
  skip_if_not_installed("quantion")

  source <- source_holding(list(c(0, 100), c(100, 100)))

  expect_false(quantion:::.already_held(source, 50, 100))
})

test_that("prefetch records the offset and the delivered length of every block", {
  skip_if_not_installed("quantion")

  source <- source_over_file(write_sample_bytes(4096))
  quantion:::.prefetch(source, rows(0, 4096))

  expect_equal(source$held, list(c(0, 4096)))
})

test_that("prefetch does not ask again for bytes it already holds", {
  skip_if_not_installed("quantion")

  source <- source_over_file(write_sample_bytes(4096))
  quantion:::.prefetch(source, rows(0, 4096))
  quantion:::.prefetch(source, rows(10, 20))

  expect_length(source$held, 1)
})

test_that("prefetch asks for bytes it does not hold", {
  skip_if_not_installed("quantion")

  source <- source_over_file(write_sample_bytes(4096))

  expect_error(quantion:::.prefetch(source, rows(10, 20)), "does not support range requests")
})
