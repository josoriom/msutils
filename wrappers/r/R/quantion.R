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
  if (typeof(bin) != "externalptr") stop("quantion: expected an external pointer")
  .Call("C_dispose_mzml", bin, PACKAGE = "quantion")
}

parse_mzml <- function(data) {
  stopifnot(is.raw(data))
  .Call("C_parse_mzml", data, PACKAGE="quantion")
}

mzml_to_ion <- function(bin, level = 12L, f32_compress = FALSE) {
  if (typeof(bin) != "externalptr") stop("quantion: expected an external pointer (MzML sample)")

  if (!is.numeric(level) || length(level) != 1 || is.na(level))
    stop("`level` must be a single number 0..22")
  lvl <- as.integer(level)
  if (lvl < 0L || lvl > 22L)
    stop("`level` must be between 0 and 22 (inclusive)")

  if (!is.logical(f32_compress) || length(f32_compress) != 1 || is.na(f32_compress))
    stop("`f32_compress` must be TRUE/FALSE")

  .Call("C_mzml_to_ion", bin, lvl, f32_compress, PACKAGE = "quantion")
}

mzml_to_ion_file <- function(input_path, output_path, level = 12L, f32_compress = FALSE) {
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

  .Call("C_mzml_to_ion_file", input_path, output_path, lvl, f32_compress, PACKAGE = "quantion")
  invisible(output_path)
}

ion_to_mzml <- function(bin) {
  if (typeof(bin) != "externalptr") stop("quantion: expected an external pointer (MzML sample)")
  .Call("C_ion_to_mzml", bin, PACKAGE = "quantion")
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
  filter_supplied <- !(missing(min_integral) && missing(min_intensity) &&
    missing(min_peak_width_points) && missing(noise) && missing(auto_noise) &&
    missing(auto_baseline) && missing(lambda) && missing(max_iterations) &&
    missing(allow_overlap) && missing(min_snr) && missing(min_r2) && missing(shape) &&
    missing(kernel_size))
  opt <- if (filter_supplied) list(
    min_integral=min_integral, min_intensity=min_intensity,
    min_peak_width_points=min_peak_width_points, noise=noise,
    auto_noise=auto_noise, auto_baseline=auto_baseline,
    lambda=lambda, max_iterations=max_iterations,
    allow_overlap=allow_overlap, min_snr=min_snr, min_r2=min_r2, shape=shape,
    kernel_size=kernel_size
  ) else NULL
  out_json <- .Call("C_get_peak", as.numeric(x), as.numeric(y), as.numeric(rt), as.numeric(range), opt, PACKAGE="quantion")
  out_json
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
  filter_supplied <- !(missing(min_integral) && missing(min_intensity) &&
    missing(min_peak_width_points) && missing(noise) && missing(auto_noise) &&
    missing(auto_baseline) && missing(lambda) && missing(max_iterations) &&
    missing(allow_overlap) && missing(min_snr) && missing(min_r2) && missing(shape) &&
    missing(kernel_size))
  opt <- if (filter_supplied) list(
    min_integral=min_integral, min_intensity=min_intensity,
    min_peak_width_points=min_peak_width_points, noise=noise,
    auto_noise=auto_noise, auto_baseline=auto_baseline,
    lambda=lambda, max_iterations=max_iterations,
    allow_overlap=allow_overlap, min_snr=min_snr, min_r2=min_r2, shape=shape,
    kernel_size=kernel_size
  ) else NULL
  cores <- .validate_cores(cores)
  out_json <- .Call("C_get_peaks_from_eic",
    bin, as.numeric(rts), as.numeric(mzs), as.numeric(ranges), as.character(id),
    as.numeric(from), as.numeric(to), opt, as.integer(cores),
    PACKAGE="quantion"
  )
  res <- out_json
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
  if (is.null(items) || !(is.list(items) || is.data.frame(items))) stop("items must be a list/data.frame")
  idxs <- suppressWarnings(as.integer(if (!is.null(items$idx)) items$idx else items$index))
  rts  <- suppressWarnings(as.numeric(items$rt))
  wins <- suppressWarnings(as.numeric(if (!is.null(items$window)) items$window else items$range))
  if (length(idxs) != length(rts) || length(wins) != length(rts)) stop("idx, rt, range length mismatch")
  if (any(is.na(idxs)) || any(idxs < 0)) stop("idx must be a non-negative integer for every item")
  if (!is.logical(auto_noise) || length(auto_noise) != 1 || is.na(auto_noise)) stop("auto_noise must be logical TRUE/FALSE")
  if (!is.logical(allow_overlap) || length(allow_overlap) != 1 || is.na(allow_overlap)) stop("allow_overlap must be logical TRUE/FALSE")
  if (!is.logical(auto_baseline) || length(auto_baseline) != 1 || is.na(auto_baseline)) stop("auto_baseline must be logical TRUE/FALSE")
  filter_supplied <- !(missing(min_integral) && missing(min_intensity) &&
    missing(min_peak_width_points) && missing(noise) && missing(auto_noise) &&
    missing(auto_baseline) && missing(lambda) && missing(max_iterations) &&
    missing(allow_overlap) && missing(min_snr) && missing(min_r2) && missing(shape))
  opt <- if (filter_supplied) list(
    min_integral=min_integral, min_intensity=min_intensity,
    min_peak_width_points=min_peak_width_points, noise=noise,
    auto_noise=auto_noise, auto_baseline=auto_baseline,
    lambda=lambda, max_iterations=max_iterations,
    allow_overlap=allow_overlap, min_snr=min_snr, min_r2=min_r2, shape=shape
  ) else NULL
  cores <- .validate_cores(cores)
  out_json <- .Call("C_get_peaks_from_chrom",
    bin, idxs, rts, wins, opt, as.integer(cores), PACKAGE="quantion"
  )
  df <- out_json
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
  if (typeof(bin) != "externalptr") stop("quantion: expected an external pointer (MzML sample)")
  stopifnot(is.numeric(targets), length(targets) == 1)
  source <- attr(bin, "quantion_source")
  if (!is.null(source)) {
    wanted <- .Call("C_plan_eic", bin, as.numeric(targets), as.numeric(from),
                    as.numeric(to), as.numeric(ppm_tolerance), as.numeric(mz_tolerance),
                    PACKAGE = "quantion")
    .prefetch(source, wanted)
  }
  .Call("C_calculate_eic",
    bin, as.numeric(targets), as.numeric(from), as.numeric(to),
    as.numeric(ppm_tolerance), as.numeric(mz_tolerance),
    PACKAGE="quantion"
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
  filter_supplied <- !(missing(min_integral) && missing(min_intensity) &&
    missing(min_peak_width_points) && missing(noise) && missing(auto_noise) &&
    missing(auto_baseline) && missing(lambda) && missing(max_iterations) &&
    missing(allow_overlap) && missing(min_snr) && missing(min_r2) && missing(shape) &&
    missing(kernel_size))
  opt <- if (filter_supplied) list(
    min_integral=min_integral, min_intensity=min_intensity,
    min_peak_width_points=min_peak_width_points, noise=noise,
    auto_noise=auto_noise, auto_baseline=auto_baseline,
    lambda=lambda, max_iterations=max_iterations,
    allow_overlap=allow_overlap, min_snr=min_snr, min_r2=min_r2, shape=shape,
    kernel_size=kernel_size
  ) else NULL
  out_json <- .Call("C_find_peaks", as.numeric(x), as.numeric(y), opt, PACKAGE="quantion")
  out_json
}

calculate_baseline <- function(y, lambda=0L, max_iterations=0L) {
  stopifnot(is.numeric(y))
  .Call("C_calculate_baseline",
        as.numeric(y),
        as.integer(lambda),
        as.integer(max_iterations),
        PACKAGE="quantion")
}

find_noise_level <- function(y) .Call("C_find_noise_level", as.numeric(y), PACKAGE="quantion")

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
        PACKAGE="quantion")
  out_json
}

draw_peak <- function(x, params) {
  stopifnot(is.numeric(x), is.list(params))
  tail <- if (is.null(params$tail)) 0 else params$tail
  .Call("C_draw_peak",
        as.numeric(x), .shape_code(params$shape),
        as.numeric(params$height), as.numeric(params$center),
        as.numeric(params$fwhm), as.numeric(tail),
        PACKAGE="quantion")
}

find_feature <- function(
  bin,
  rt, mz, window, id = NULL,
  scan_ppm = 10, scan_mz = 0.003, eic_ppm = 20, eic_mz = 0.005,
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
    if (length(id) != length(rt)) stop("id must have the same length as rt")
  }
  if (length(mz) != length(rt) || length(window) != length(rt))
    stop("rt, mz and window must have the same length")

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
    PACKAGE = "quantion"
  )

  res <- out_json
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
    opt, as.integer(cores),
    PACKAGE = "quantion"
  )

  df <- out_json
  if (!is.data.frame(df)) df <- as.data.frame(df)
  want <- c("mz","rt","from","to","intensity","integral","n_points")
  present <- intersect(want, names(df))
  df <- df[, c(present, setdiff(names(df), present)), drop = FALSE]
  rownames(df) <- NULL
  df
}

parse_ion <- function(input, max_cache_size = 0) {
  if (is.character(input)) {
    looks_like_url <- length(input) == 1 && !is.na(input) && grepl("^https?://", input)
    if (looks_like_url) return(parse_ion_remote(input, max_cache_size))
    return(parse_ion_path(input, max_cache_size))
  }
  parse_ion_raw(input, max_cache_size)
}

parse_ion_raw <- function(input, max_cache_size = 0) {
  if (!is.raw(input)) stop("`input` must be a raw vector of ion bytes")
  if (!is.numeric(max_cache_size) || length(max_cache_size) != 1 || is.na(max_cache_size) || max_cache_size < 0)
    stop("`max_cache_size` must be a single non-negative number")
  .Call("C_parse_ion", input, as.numeric(max_cache_size), PACKAGE = "quantion")
}

parse_ion_path <- function(path, max_cache_size = 0) {
  if (!is.character(path) || length(path) != 1 || is.na(path) || nchar(path) == 0)
    stop("`path` must be a single non-empty string")
  full_path <- path.expand(path)
  if (!file.exists(full_path)) stop("`path` does not point to an existing file: ", path)
  if (!is.numeric(max_cache_size) || length(max_cache_size) != 1 || is.na(max_cache_size) || max_cache_size < 0)
    stop("`max_cache_size` must be a single non-negative number")
  .Call("C_parse_ion_path", full_path, as.numeric(max_cache_size), PACKAGE = "quantion")
}

.fetch_range <- function(url, offset, span) {
  destination <- tempfile()
  on.exit(unlink(destination), add = TRUE)
  last <- offset + span - 1
  status <- utils::download.file(
    url, destination, method = "libcurl", quiet = TRUE, mode = "wb",
    headers = c(Range = sprintf("bytes=%.0f-%.0f", offset, last))
  )
  if (!identical(status, 0L)) stop("quantion: the range request failed for ", url)
  bytes <- readBin(destination, "raw", file.info(destination)$size)
  if (length(bytes) != span)
    stop(sprintf("quantion: the server returned %d bytes for a %d byte range; it does not support range requests",
                 length(bytes), span))
  bytes
}

.LARGEST_GAP <- 131072

.gap_for <- function(total) {
  min(.LARGEST_GAP, total / 8)
}

.merge_ranges <- function(ranges, gap = .LARGEST_GAP) {
  wanted <- ranges[ranges[, 2] > 0, , drop = FALSE]
  if (nrow(wanted) == 0) return(wanted)
  wanted <- wanted[order(wanted[, 1]), , drop = FALSE]

  merged <- list(c(wanted[1, 1], wanted[1, 2]))
  for (row in seq_len(nrow(wanted))[-1]) {
    last <- merged[[length(merged)]]
    if (wanted[row, 1] - (last[1] + last[2]) <= gap) {
      end <- max(last[1] + last[2], wanted[row, 1] + wanted[row, 2])
      merged[[length(merged)]] <- c(last[1], end - last[1])
    } else {
      merged[[length(merged) + 1]] <- c(wanted[row, 1], wanted[row, 2])
    }
  }
  do.call(rbind, merged)
}

.already_held <- function(source, offset, span) {
  for (held in source$held) {
    if (offset >= held[1] && offset + span <= held[1] + held[2]) return(TRUE)
  }
  FALSE
}

.prefetch <- function(source, ranges) {
  wanted <- ranges[!apply(ranges, 1, function(row) .already_held(source, row[1], row[2])), , drop = FALSE]
  if (nrow(wanted) == 0) return(invisible(NULL))

  spans <- .merge_ranges(wanted, .gap_for(source$total))
  for (row in seq_len(nrow(spans))) {
    offset <- spans[row, 1]
    span <- spans[row, 2]
    if (.already_held(source, offset, span)) next
    bytes <- .fetch_range(source$url, offset, span)
    .Call("C_store_add", source$store, offset, bytes, PACKAGE = "quantion")
    source$held[[length(source$held) + 1]] <- c(offset, length(bytes))
  }
  invisible(NULL)
}

parse_ion_remote <- function(url, max_cache_size = 0) {
  if (!is.character(url) || length(url) != 1 || is.na(url) || nchar(url) == 0)
    stop("`url` must be a single non-empty string")
  if (!startsWith(tolower(url), "http://") && !startsWith(tolower(url), "https://"))
    stop("`url` must start with http:// or https://")

  source <- new.env(parent = emptyenv())
  source$url <- url
  source$store <- .Call("C_store_new", PACKAGE = "quantion")
  source$held <- list()
  source$total <- 0

  header <- .fetch_range(url, 0, 1024)
  .Call("C_store_add", source$store, 0, header, PACKAGE = "quantion")
  source$held[[1]] <- c(0, length(header))

  planned <- .Call("C_plan_open", header, PACKAGE = "quantion")
  source$total <- if (nrow(planned)) max(planned[, 1] + planned[, 2]) else 0
  .prefetch(source, planned)

  handle <- .Call("C_parse_ion_source", source$store, as.numeric(max_cache_size),
                  PACKAGE = "quantion")
  attr(handle, "quantion_source") <- source
  handle
}

.QUERY_RT_RANGE   <- 0L
.QUERY_CLOSEST_RT <- 1L
.QUERY_MZ_RANGE   <- 2L
.QUERY_CLOSEST_MZ <- 3L

get_scans <- function(bin, rt_from = NULL, rt_to = NULL, rt = NULL,
                      mz_from = NULL, mz_to = NULL, mz = NULL,
                      level = 1L) {
  if (typeof(bin) != "externalptr") stop("quantion: expected an external pointer (MzML sample)")
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

  source <- attr(bin, "quantion_source")
  if (!is.null(source)) {
    wanted <- .Call("C_plan_scans", bin, as.integer(query_type), a, b,
                    as.integer(level), PACKAGE = "quantion")
    .prefetch(source, wanted)
  }

  scans <- .Call("C_get_scans",
    bin, as.integer(query_type), a, b, as.integer(level),
    PACKAGE = "quantion"
  )
  if (query_type %in% c(.QUERY_CLOSEST_RT, .QUERY_CLOSEST_MZ)) {
    if (nrow(scans) == 0) {
      NULL
    } else {
      list(
        rt = scans$rt[[1]],
        mz = as.list(scans$mz[[1]]),
        intensity = as.list(scans$intensity[[1]]),
        metadata = as.list(scans$metadata[1, ])
      )
    }
  } else {
    scans
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
    opt, as.integer(cores),
    PACKAGE = "quantion"
  )

  df <- out_json
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

