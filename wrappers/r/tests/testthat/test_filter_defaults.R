default_of <- function(work, name) eval(formals(work)[[name]], asNamespace("quantion"))

peak_finders <- function() list(
  find_peaks = quantion::find_peaks,
  get_peak = quantion::get_peak,
  find_feature = quantion::find_feature,
  get_features = quantion::get_features,
  get_peaks_from_eic = quantion::get_peaks_from_eic,
  get_peaks_from_chrom = quantion::get_peaks_from_chrom
)

shared_settings <- c("min_intensity", "min_peak_width_points", "auto_noise",
                     "auto_baseline", "min_snr", "min_r2", "shape")

test_that("get_peaks_from_chrom starts with the automatic noise and baseline on", {
  skip_if_not_installed("quantion")

  expect_true(default_of(quantion::get_peaks_from_chrom, "auto_noise"))
  expect_true(default_of(quantion::get_peaks_from_chrom, "auto_baseline"))
})

test_that("every peak finder starts from the same filter defaults", {
  skip_if_not_installed("quantion")

  finders <- peak_finders()
  for (name in shared_settings) {
    values <- lapply(finders, default_of, name = name)
    first <- values[[1]]
    for (value in values) expect_identical(value, first, info = name)
  }
})

test_that("the shared filter defaults are the ones the package publishes", {
  skip_if_not_installed("quantion")

  for (name in shared_settings) {
    expect_identical(
      default_of(quantion::get_peaks_from_chrom, name),
      quantion:::.DEFAULTS[[name]],
      info = name
    )
  }
})
