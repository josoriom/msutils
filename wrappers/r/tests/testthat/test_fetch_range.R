test_that("a range request gives back exactly the bytes it asked for", {
  skip_if_not_installed("quantion")

  path <- write_sample_bytes(5000)
  bytes <- quantion:::.fetch_range(file_url(path), 0, 5000)

  expect_true(is.raw(bytes))
  expect_length(bytes, 5000)
  expect_identical(bytes, readBin(path, "raw", 5000))
})

test_that("a source that sends the whole file for a range is refused", {
  skip_if_not_installed("quantion")

  path <- write_sample_bytes(5000)

  expect_error(
    quantion:::.fetch_range(file_url(path), 100, 256),
    "does not support range requests"
  )
})

test_that("the refusal names both the asked size and the delivered size", {
  skip_if_not_installed("quantion")

  path <- write_sample_bytes(5000)
  message <- tryCatch(
    quantion:::.fetch_range(file_url(path), 100, 256),
    error = function(e) conditionMessage(e)
  )

  expect_match(message, "5000 bytes")
  expect_match(message, "256 byte range")
})

test_that("a range request that cannot be made stops instead of returning bytes", {
  skip_if_not_installed("quantion")

  expect_error(suppressWarnings(quantion:::.fetch_range(file_url(tempfile()), 0, 8)))
})
