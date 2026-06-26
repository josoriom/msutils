.DEFAULTS <- list(
  eic_ppm_tol          = 10.0,
  eic_mz_tol           = 0.005,
  grid_start           = 40,
  grid_end             = 1000,
  grid_step            = 0.005,
  group_ppm_tol        = 5.0,
  group_mz_tol         = 0.0025,
  group_rt_tol         = 0.05,
  frequency            = 1L,
  min_intensity        = 500,
  min_peak_width_points = 3L,
  auto_noise           = TRUE,
  auto_baseline        = TRUE,
  min_snr              = 2,
  min_r2               = 0.0,
  shape                = "emg",
  kernel_size          = 0L
)

dispose <- function(bin) {
  if (typeof(bin) != "externalptr") stop("msutils: expected an external pointer")
  .Call("C_dispose_mzml", bin, PACKAGE = "msutils")
}

parse_mzml <- function(data) {
  stopifnot(is.raw(data))
  .Call("C_parse_mzml", data, PACKAGE="msutils")
}

mzml_to_ion <- function(bin, level = 12L, f32_compress = FALSE) {
  if (typeof(bin) != "externalptr") stop("msutils: expected an external pointer (MzML sample)")

  if (!is.numeric(level) || length(level) != 1 || is.na(level))
    stop("`level` must be a single number 0..22")
  lvl <- as.integer(level)
  if (lvl < 0L || lvl > 22L)
    stop("`level` must be between 0 and 22 (inclusive)")

  if (!is.logical(f32_compress) || length(f32_compress) != 1 || is.na(f32_compress))
    stop("`f32_compress` must be TRUE/FALSE")

  .Call("C_mzml_to_ion", bin, lvl, f32_compress, PACKAGE = "msutils")
}

mzml_to_ion_file <- function(input_path, output_path, level = 12L, f32_compress = FALSE, section_on_disk = FALSE) {
  if (!is.character(input_path) || length(input_path) != 1 || is.na(input_path))
    stop("`input_path` must be a single string")
  if (!is.character(output_path) || length(output_path) != 1 || is.na(output_path))
    stop("`output_path` must be a single string")

  if (!is.numeric(level) || length(level) != 1 || is.na(level))
    stop("`level` must be a single number 0..22")
  lvl <- as.integer(level)
  if (lvl < 0L || lvl > 22L)
    stop("`level` must be between 0 and 22 (inclusive)")

  if (!is.logical(f32_compress) || length(f32_compress) != 1 || is.na(f32_compress))
    stop("`f32_compress` must be TRUE/FALSE")
  if (!is.logical(section_on_disk) || length(section_on_disk) != 1 || is.na(section_on_disk))
    stop("`section_on_disk` must be TRUE/FALSE")

  .Call("C_mzml_to_ion_file", input_path, output_path, lvl, f32_compress, section_on_disk, PACKAGE = "msutils")
  invisible(output_path)
}

ion_to_json <- function(bin) {
  if (typeof(bin) != "externalptr") stop("msutils: expected an external pointer (MzML sample)")
  .Call("C_ion_to_json", bin, PACKAGE="msutils")
}

ion_to_mzml <- function(bin) {
  if (typeof(bin) != "externalptr") stop("msutils: expected an external pointer (MzML sample)")
  .Call("C_ion_to_mzml", bin, PACKAGE = "msutils")
}

get_peak <- function(
  x, y, rt, range,
  min_integral          = NaN,
  min_intensity         = .DEFAULTS$min_intensity,
  min_peak_width_points = .DEFAULTS$min_peak_width_points,
  noise                 = NaN,
  auto_noise            = .DEFAULTS$auto_noise,
  auto_baseline         = .DEFAULTS$auto_baseline,
  lambda                = 0L,
  max_iterations        = 0L,
  allow_overlap         = FALSE,
  min_snr               = .DEFAULTS$min_snr,
  min_r2                = .DEFAULTS$min_r2,
  shape                 = .DEFAULTS$shape,
  kernel_size           = .DEFAULTS$kernel_size
) {
  stopifnot(is.numeric(x), is.numeric(y))
  if (length(x) != length(y) || length(x) < 3) stop("x and y must have the same length (>= 3)")
  if (!is.logical(auto_noise) || length(auto_noise) != 1 || is.na(auto_noise)) stop("auto_noise must be logical TRUE/FALSE")
  if (!is.logical(allow_overlap) || length(allow_overlap) != 1 || is.na(allow_overlap)) stop("allow_overlap must be logical TRUE/FALSE")
  if (!is.logical(auto_baseline) || length(auto_baseline) != 1 || is.na(auto_baseline)) stop("auto_baseline must be logical TRUE/FALSE")
  opt <- list(
    min_integral=min_integral, min_intensity=min_intensity,
    min_peak_width_points=min_peak_width_points, noise=noise,
    auto_noise=auto_noise, auto_baseline=auto_baseline,
    lambda=lambda, max_iterations=max_iterations,
    allow_overlap=allow_overlap, min_snr=min_snr, min_r2=min_r2, shape=shape,
    kernel_size=kernel_size
  )
  out_json <- .Call("C_get_peak", as.numeric(x), as.numeric(y), as.numeric(rt), as.numeric(range), opt, PACKAGE="msutils")
  jsonlite::fromJSON(out_json, simplifyVector=TRUE)
}

get_peaks_from_eic <- function(
  bin, df, from=0.5, to=5, cores=1L,
  min_integral          = NaN,
  min_intensity         = .DEFAULTS$min_intensity,
  min_peak_width_points = .DEFAULTS$min_peak_width_points,
  noise                 = NaN,
  auto_noise            = .DEFAULTS$auto_noise,
  auto_baseline         = .DEFAULTS$auto_baseline,
  lambda                = 0L,
  max_iterations        = 0L,
  allow_overlap         = FALSE,
  min_snr               = .DEFAULTS$min_snr,
  min_r2                = .DEFAULTS$min_r2,
  shape                 = .DEFAULTS$shape,
  kernel_size           = .DEFAULTS$kernel_size
) {
  stopifnot(typeof(bin) == "externalptr")
  if (!is.data.frame(df)) stop("`df` must be a data.frame")
  req <- c("id","rt","mz","ranges"); miss <- setdiff(req, names(df))
  if (length(miss)) stop(paste0("missing columns: ", paste(miss, collapse=", ")))
  id <- as.character(df$id)
  rts <- suppressWarnings(as.numeric(df$rt))
  mzs <- suppressWarnings(as.numeric(df$mz))
  ranges <- suppressWarnings(as.numeric(df$ranges))
  n <- length(id)
  if (!(length(rts)==n && length(mzs)==n && length(ranges)==n)) stop("id, rt, mz, ranges must have the same length")
  id[is.na(id)] <- ""
  if (!is.logical(auto_noise) || length(auto_noise) != 1 || is.na(auto_noise)) stop("auto_noise must be logical TRUE/FALSE")
  if (!is.logical(allow_overlap) || length(allow_overlap) != 1 || is.na(allow_overlap)) stop("allow_overlap must be logical TRUE/FALSE")
  if (!is.logical(auto_baseline) || length(auto_baseline) != 1 || is.na(auto_baseline)) stop("auto_baseline must be logical TRUE/FALSE")
  opt <- list(
    min_integral=min_integral, min_intensity=min_intensity,
    min_peak_width_points=min_peak_width_points, noise=noise,
    auto_noise=auto_noise, auto_baseline=auto_baseline,
    lambda=lambda, max_iterations=max_iterations,
    allow_overlap=allow_overlap, min_snr=min_snr, min_r2=min_r2, shape=shape,
    kernel_size=kernel_size
  )
  cores <- .validate_cores(cores)
  out_json <- .Call("C_get_peaks_from_eic",
    bin, as.numeric(rts), as.numeric(mzs), as.numeric(ranges), as.character(id),
    as.numeric(from), as.numeric(to), opt, as.integer(cores),
    PACKAGE="msutils"
  )
  res <- jsonlite::fromJSON(out_json, simplifyVector=TRUE)
  if (!is.data.frame(res)) res <- as.data.frame(res)
  want <- c("id","mz","ort","rt","from","to","intensity","integral")
  present <- intersect(want, names(res)); extras <- setdiff(names(res), present)
  res <- res[, c(present, extras), drop=FALSE]
  rownames(res) <- NULL
  res
}

get_peaks_from_chrom <- function(
  bin, items, cores=1L,
  min_integral          = NaN,
  min_intensity         = .DEFAULTS$min_intensity,
  min_peak_width_points = .DEFAULTS$min_peak_width_points,
  noise                 = NaN,
  auto_noise            = FALSE,
  auto_baseline         = FALSE,
  lambda                = 0L,
  max_iterations        = 0L,
  allow_overlap         = FALSE,
  min_snr               = .DEFAULTS$min_snr,
  min_r2                = .DEFAULTS$min_r2,
  shape                 = .DEFAULTS$shape
) {
  stopifnot(typeof(bin) == "externalptr")
  if (is.null(items) || !(is.list(items) || is.data.frame(items))) stop("items must be a list/data.frame")
  idxs <- suppressWarnings(as.integer(if (!is.null(items$idx)) items$idx else items$index))
  rts  <- suppressWarnings(as.numeric(items$rt))
  wins <- suppressWarnings(as.numeric(if (!is.null(items$window)) items$window else items$range))
  if (length(idxs) != length(rts) || length(wins) != length(rts)) stop("idx, rt, range length mismatch")
  if (!is.logical(auto_noise) || length(auto_noise) != 1 || is.na(auto_noise)) stop("auto_noise must be logical TRUE/FALSE")
  if (!is.logical(allow_overlap) || length(allow_overlap) != 1 || is.na(allow_overlap)) stop("allow_overlap must be logical TRUE/FALSE")
  if (!is.logical(auto_baseline) || length(auto_baseline) != 1 || is.na(auto_baseline)) stop("auto_baseline must be logical TRUE/FALSE")
  opt <- list(
    min_integral=min_integral, min_intensity=min_intensity,
    min_peak_width_points=min_peak_width_points, noise=noise,
    auto_noise=auto_noise, auto_baseline=auto_baseline,
    lambda=lambda, max_iterations=max_iterations,
    allow_overlap=allow_overlap, min_snr=min_snr, min_r2=min_r2, shape=shape
  )
  cores <- .validate_cores(cores)
  out_json <- .Call("C_get_peaks_from_chrom",
    bin, idxs, rts, wins, opt, as.integer(cores), PACKAGE="msutils"
  )
  df <- jsonlite::fromJSON(out_json, simplifyVector=TRUE)
  if (!is.data.frame(df)) df <- as.data.frame(df)
  if (!"ort" %in% names(df)) stop("internal error: 'ort' missing from result")
  want <- c("index","id","ort","rt","from","to","intensity","integral","total_area")
  present <- want[want %in% names(df)]; extras <- setdiff(names(df), present)
  df <- df[, c(present, extras), drop=FALSE]
  if ("index" %in% names(df)) df <- df[order(df$index), , drop=FALSE]
  rownames(df) <- NULL
  df
}

calculate_eic <- function(bin, targets, from, to, ppm_tolerance=20, mz_tolerance=0.005) {
  if (typeof(bin) != "externalptr") stop("msutils: expected an external pointer (MzML sample)")
  stopifnot(is.numeric(targets), length(targets) == 1)
  .Call("C_calculate_eic",
    bin, as.numeric(targets), as.numeric(from), as.numeric(to),
    as.numeric(ppm_tolerance), as.numeric(mz_tolerance),
    PACKAGE="msutils"
  )
}

find_peaks <- function(
  x, y,
  min_integral          = NaN,
  min_intensity         = .DEFAULTS$min_intensity,
  min_peak_width_points = .DEFAULTS$min_peak_width_points,
  noise                 = NaN,
  auto_noise            = .DEFAULTS$auto_noise,
  auto_baseline         = .DEFAULTS$auto_baseline,
  lambda                = 0L,
  max_iterations        = 0L,
  allow_overlap         = FALSE,
  min_snr               = .DEFAULTS$min_snr,
  min_r2                = .DEFAULTS$min_r2,
  shape                 = .DEFAULTS$shape,
  kernel_size           = .DEFAULTS$kernel_size
) {
  stopifnot(is.numeric(x), is.numeric(y))
  if (length(x) != length(y) || length(x) < 3) stop("x and y must have the same length >= 3")
  if (!is.logical(auto_noise) || length(auto_noise) != 1 || is.na(auto_noise)) stop("auto_noise must be logical TRUE/FALSE")
  if (!is.logical(allow_overlap) || length(allow_overlap) != 1 || is.na(allow_overlap)) stop("allow_overlap must be logical TRUE/FALSE")
  if (!is.logical(auto_baseline) || length(auto_baseline) != 1 || is.na(auto_baseline)) stop("auto_baseline must be logical TRUE/FALSE")
  opt <- list(
    min_integral=min_integral, min_intensity=min_intensity,
    min_peak_width_points=min_peak_width_points, noise=noise,
    auto_noise=auto_noise, auto_baseline=auto_baseline,
    lambda=lambda, max_iterations=max_iterations,
    allow_overlap=allow_overlap, min_snr=min_snr, min_r2=min_r2, shape=shape,
    kernel_size=kernel_size
  )
  out_json <- .Call("C_find_peaks", as.numeric(x), as.numeric(y), opt, PACKAGE="msutils")
  jsonlite::fromJSON(out_json, simplifyVector=TRUE)
}

calculate_baseline <- function(y, lambda=0L, max_iterations=0L) {
  stopifnot(is.numeric(y))
  .Call("C_calculate_baseline",
        as.numeric(y),
        as.integer(lambda),
        as.integer(max_iterations),
        PACKAGE="msutils")
}

.shape_code <- function(shape) {
  s <- tolower(as.character(shape))
  if (s == "gaussian") return(0L)
  if (s == "emg")      return(1L)
  stop("shape must be 'gaussian' or 'emg'")
}

fit_peak <- function(x, y, rt, intensity, shape = "emg") {
  stopifnot(is.numeric(x), is.numeric(y), length(x) == length(y), length(x) >= 5)
  out_json <- .Call("C_fit_peak",
        as.numeric(x), as.numeric(y),
        as.numeric(rt), as.numeric(intensity), .shape_code(shape),
        PACKAGE="msutils")
  jsonlite::fromJSON(out_json, simplifyVector = TRUE)
}

draw_peak <- function(x, params) {
  stopifnot(is.numeric(x), is.list(params))
  tail <- if (is.null(params$tail)) 0 else params$tail
  .Call("C_draw_peak",
        as.numeric(x), .shape_code(params$shape),
        as.numeric(params$height), as.numeric(params$center),
        as.numeric(params$fwhm), as.numeric(tail),
        PACKAGE="msutils")
}

find_feature <- function(
  bin,
  rt, mz, window, id = NULL,
  scan_ppm, scan_mz, eic_ppm, eic_mz,
  cores                 = 1L,
  min_integral          = NaN,
  min_intensity         = .DEFAULTS$min_intensity,
  min_peak_width_points = .DEFAULTS$min_peak_width_points,
  noise                 = NaN,
  auto_noise            = .DEFAULTS$auto_noise,
  auto_baseline         = .DEFAULTS$auto_baseline,
  lambda                = 0L,
  max_iterations        = 0L,
  allow_overlap         = FALSE,
  min_snr               = .DEFAULTS$min_snr,
  min_r2                = .DEFAULTS$min_r2,
  shape                 = .DEFAULTS$shape
) {
  stopifnot(typeof(bin) == "externalptr")
  if (!is.numeric(rt)) stop("rt must be numeric")
  if (!is.numeric(mz)) stop("mz must be numeric")
  if (!is.numeric(window)) stop("window must be numeric")

  if (!is.numeric(scan_ppm) || length(scan_ppm) != 1) stop("scan_ppm must be a single numeric")
  if (!is.numeric(scan_mz)  || length(scan_mz)  != 1) stop("scan_mz must be a single numeric")
  if (!is.numeric(eic_ppm)  || length(eic_ppm)  != 1) stop("eic_ppm must be a single numeric")
  if (!is.numeric(eic_mz)   || length(eic_mz)   != 1) stop("eic_mz must be a single numeric")

  if (!is.logical(auto_noise) || length(auto_noise) != 1 || is.na(auto_noise)) stop("auto_noise must be logical TRUE/FALSE")
  if (!is.logical(allow_overlap) || length(allow_overlap) != 1 || is.na(allow_overlap)) stop("allow_overlap must be logical TRUE/FALSE")
  if (!is.logical(auto_baseline) || length(auto_baseline) != 1 || is.na(auto_baseline)) stop("auto_baseline must be logical TRUE/FALSE")

  if (is.null(id)) {
    id <- rep("", length(rt))
  } else {
    id <- as.character(id)
  }

  opt <- list(
    min_integral=min_integral,
    min_intensity=min_intensity,
    min_peak_width_points=min_peak_width_points,
    noise=noise,
    auto_noise=auto_noise,
    auto_baseline=auto_baseline,
    lambda=lambda,
    max_iterations=max_iterations,
    allow_overlap=allow_overlap,
    min_snr=min_snr, min_r2=min_r2, shape=shape
  )

  cores <- .validate_cores(cores)

  out_json <- .Call(
    "C_find_feature",
    bin,
    as.numeric(rt), as.numeric(mz), as.numeric(window),
    id,
    as.numeric(scan_ppm), as.numeric(scan_mz),
    as.numeric(eic_ppm), as.numeric(eic_mz),
    opt, as.integer(cores),
    PACKAGE = "msutils"
  )

  res <- jsonlite::fromJSON(out_json, simplifyVector = TRUE)
  if (!is.data.frame(res)) res <- as.data.frame(res)
  want <- c("id","mz","ort","rt","from","to","intensity","integral")
  present <- intersect(want, names(res)); extras <- setdiff(names(res), present)
  res <- res[, c(present, extras), drop = FALSE]
  rownames(res) <- NULL
  res
}

find_features <- function(
  data, from = 0, to = 10,
  ppm_tolerance         = .DEFAULTS$eic_ppm_tol,
  mz_tolerance          = .DEFAULTS$eic_mz_tol,
  grid_start            = .DEFAULTS$grid_start,
  grid_end              = .DEFAULTS$grid_end,
  grid_step_ppm         = .DEFAULTS$grid_step,
  cores                 = 1L,
  use_gpu               = FALSE,
  batch_size            = 0L,
  min_integral          = NaN,
  min_intensity         = .DEFAULTS$min_intensity,
  min_peak_width_points = .DEFAULTS$min_peak_width_points,
  noise                 = NaN,
  auto_noise            = .DEFAULTS$auto_noise,
  auto_baseline         = .DEFAULTS$auto_baseline,
  lambda                = 0L,
  max_iterations        = 0L,
  allow_overlap         = FALSE,
  min_snr               = .DEFAULTS$min_snr,
  min_r2                = .DEFAULTS$min_r2,
  shape                 = .DEFAULTS$shape
) {
  stopifnot(typeof(data) == "externalptr")
  if (!is.logical(use_gpu) || length(use_gpu) != 1 || is.na(use_gpu)) stop("use_gpu must be logical TRUE/FALSE")
  cores <- .validate_cores(cores)

  opt <- list(
    min_integral=min_integral, min_intensity=min_intensity,
    min_peak_width_points=min_peak_width_points, noise=noise,
    auto_noise=auto_noise, auto_baseline=auto_baseline,
    lambda=lambda, max_iterations=max_iterations,
    allow_overlap=allow_overlap, min_snr=min_snr, min_r2=min_r2, shape=shape
  )

  out_json <- .Call(
    "C_find_features", data,
    as.numeric(from), as.numeric(to),
    as.numeric(ppm_tolerance), as.numeric(mz_tolerance),
    as.numeric(grid_start), as.numeric(grid_end), as.numeric(grid_step_ppm),
    opt, as.integer(cores), as.logical(use_gpu), as.integer(batch_size),
    PACKAGE = "msutils"
  )

  df <- jsonlite::fromJSON(out_json, simplifyVector = TRUE)
  if (!is.data.frame(df)) df <- as.data.frame(df)
  want <- c("mz","rt","from","to","intensity","integral","n_points")
  present <- intersect(want, names(df))
  df <- df[, c(present, setdiff(names(df), present)), drop = FALSE]
  rownames(df) <- NULL
  df
}

parse_ion <- function(bin, max_cache_size = 0) {
  if (!is.raw(bin)) stop("`bin` must be a raw vector (ion bytes)")
  if (!is.numeric(max_cache_size) || length(max_cache_size) != 1 || is.na(max_cache_size) || max_cache_size < 0)
    stop("`max_cache_size` must be a single non-negative number")
  .Call("C_parse_ion", bin, as.numeric(max_cache_size), PACKAGE = "msutils")
}

parse_ion_url <- function(url, max_cache_size = 0) {
  if (!is.character(url) || length(url) != 1 || is.na(url) || nchar(url) == 0)
    stop("`url` must be a single non-empty string")
  if (!startsWith(tolower(url), "http://") && !startsWith(tolower(url), "https://"))
    stop("`url` must start with http:// or https://")
  if (!is.numeric(max_cache_size) || length(max_cache_size) != 1 || is.na(max_cache_size) || max_cache_size < 0)
    stop("`max_cache_size` must be a single non-negative number")

  .Call("C_parse_ion_url", url, as.numeric(max_cache_size), PACKAGE = "msutils")
}

.QUERY_RT_RANGE   <- 0L
.QUERY_CLOSEST_RT <- 1L
.QUERY_MZ_RANGE   <- 2L
.QUERY_CLOSEST_MZ <- 3L

get_scans <- function(bin, rt_from = NULL, rt_to = NULL, rt = NULL,
                      mz_from = NULL, mz_to = NULL, mz = NULL,
                      level = 1L) {
  if (typeof(bin) != "externalptr") stop("msutils: expected an external pointer (MzML sample)")
  if (!is.numeric(level) || length(level) != 1) stop("level must be a single numeric")

  provided <- c(!is.null(rt_from) && !is.null(rt_to),
                !is.null(rt),
                !is.null(mz_from) && !is.null(mz_to),
                !is.null(mz))
  if (sum(provided) != 1) stop("get_scans: exactly one of {rt_from+rt_to, rt, mz_from+mz_to, mz} must be provided")

  if (provided[1]) {
    a <- as.numeric(rt_from); b <- as.numeric(rt_to)
    if (!is.finite(a) || !is.finite(b)) stop("rt_from and rt_to must be finite numerics")
    query_type <- .QUERY_RT_RANGE
  } else if (provided[2]) {
    a <- as.numeric(rt)
    if (!is.finite(a)) stop("rt must be a single finite numeric")
    query_type <- .QUERY_CLOSEST_RT; b <- NaN
  } else if (provided[3]) {
    a <- as.numeric(mz_from); b <- as.numeric(mz_to)
    if (!is.finite(a) || !is.finite(b)) stop("mz_from and mz_to must be finite numerics")
    query_type <- .QUERY_MZ_RANGE
  } else {
    a <- as.numeric(mz)
    if (!is.finite(a) || a <= 0) stop("mz must be a single positive finite numeric")
    query_type <- .QUERY_CLOSEST_MZ; b <- NaN
  }

  out_json <- .Call("C_get_scans",
    bin, as.integer(query_type), a, b, as.integer(level),
    PACKAGE = "msutils"
  )
  if (query_type %in% c(.QUERY_CLOSEST_RT, .QUERY_CLOSEST_MZ)) {
    scans <- jsonlite::fromJSON(out_json, simplifyVector = FALSE)
    if (length(scans) == 0) NULL else scans[[1]]
  } else {
    jsonlite::fromJSON(out_json, simplifyVector = TRUE)
  }
}

get_features <- function(
  dir_path, from = 0, to = 10,
  eic_ppm_tol           = .DEFAULTS$eic_ppm_tol,
  eic_mz_tol            = .DEFAULTS$eic_mz_tol,
  grid_start            = .DEFAULTS$grid_start,
  grid_end              = .DEFAULTS$grid_end,
  grid_step             = .DEFAULTS$grid_step,
  group_ppm_tol         = .DEFAULTS$group_ppm_tol,
  group_mz_tol          = .DEFAULTS$group_mz_tol,
  group_rt_tol          = .DEFAULTS$group_rt_tol,
  frequency             = .DEFAULTS$frequency,
  cores                 = 1L,
  use_gpu               = FALSE,
  batch_size            = 0L,
  min_integral          = NaN,
  min_intensity         = .DEFAULTS$min_intensity,
  min_peak_width_points = .DEFAULTS$min_peak_width_points,
  noise                 = NaN,
  auto_noise            = .DEFAULTS$auto_noise,
  auto_baseline         = .DEFAULTS$auto_baseline,
  lambda                = 0L,
  max_iterations        = 0L,
  allow_overlap         = FALSE,
  min_snr               = .DEFAULTS$min_snr,
  min_r2                = .DEFAULTS$min_r2,
  shape                 = .DEFAULTS$shape
) {
  if (!is.character(dir_path) || length(dir_path) != 1) stop("dir_path must be a single string")
  if (!is.logical(use_gpu) || length(use_gpu) != 1 || is.na(use_gpu)) stop("use_gpu must be logical TRUE/FALSE")
  cores <- .validate_cores(cores)

  opt <- list(
    min_integral=min_integral, min_intensity=min_intensity,
    min_peak_width_points=min_peak_width_points, noise=noise,
    auto_noise=auto_noise, auto_baseline=auto_baseline,
    lambda=lambda, max_iterations=max_iterations,
    allow_overlap=allow_overlap, min_snr=min_snr, min_r2=min_r2, shape=shape
  )

  out_json <- .Call("C_get_features",
    dir_path,
    as.numeric(from), as.numeric(to),
    as.numeric(eic_ppm_tol), as.numeric(eic_mz_tol),
    as.numeric(grid_start), as.numeric(grid_end), as.numeric(grid_step),
    as.numeric(group_ppm_tol), as.numeric(group_mz_tol), as.numeric(group_rt_tol),
    as.integer(frequency),
    opt, as.integer(cores), as.logical(use_gpu), as.integer(batch_size),
    PACKAGE = "msutils"
  )

  df <- jsonlite::fromJSON(out_json, simplifyVector = TRUE)
  if (!is.data.frame(df)) df <- as.data.frame(df)
  rownames(df) <- NULL
  df
}

.validate_cores <- function(cores) {
  if (!is.numeric(cores) || length(cores) != 1L || is.na(cores)) stop("cores must be a single number")

  if (cores != floor(cores)) stop("cores must be a single number")

  cores <- as.integer(cores)
  if (is.na(cores) || cores < 1L) stop("cores must be a single number")

  maxc <- NA_integer_
  if (requireNamespace("parallel", quietly = TRUE)) {
    maxc <- suppressWarnings(tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA_integer_))
    if (is.na(maxc) || maxc < 1L) {
      maxc <- suppressWarnings(tryCatch(parallel::detectCores(logical = TRUE), error = function(e) NA_integer_))
    }
  }
  if (!is.na(maxc) && cores > maxc) {
    warning(sprintf("cores (%d) exceeds detected cores (%d)", cores, maxc))
  }

  cores
}

ion_to_df <- function(bin) {
  if (typeof(bin) != "externalptr") stop("msutils: expected an external pointer")
  x <- jsonlite::fromJSON(ion_to_json(bin), simplifyVector = TRUE)
  if (!is.null(x$Err)) stop(x$Err)

  root <- if (!is.null(x$Ok)) x$Ok else x
  run <- root$run
  if (is.null(run)) return(root)

  process_node <- function(meta_list, data_node) {
    if (is.null(meta_list$spectra) && is.null(meta_list$chromatograms)) return(data_node)
    meta <- if (!is.null(meta_list$spectra)) meta_list$spectra else meta_list$chromatograms
    df <- as.data.frame(meta, stringsAsFactors = FALSE)

    for (col in names(data_node)) {
      if (is.list(data_node[[col]])) df[[col]] <- I(data_node[[col]])
      else df[[col]] <- data_node[[col]]
    }
    df
  }

  run$spectra <- process_node(run$spectrum_list, run$spectra)
  run$chromatograms <- process_node(run$chromatogram_list, run$chromatograms)

  run$spectrum_list <- run$chromatogram_list <- root$spectra <- root$chromatograms <- NULL
  root$run <- run
  root
}
