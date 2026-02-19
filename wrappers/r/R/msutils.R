.pack_opts <- function(options) {
  if (is.null(options)) return(NULL)
  if (is.raw(options)) {
    if (length(options) != 64L) stop("options raw blob must be length 64")
    return(options)
  }
  if (!is.list(options) || is.null(names(options))) stop("options must be a named list in snake_case")
  allow <- c(
    "integral_threshold","intensity_threshold","width_threshold","noise",
    "auto_noise","auto_baseline","baseline_window","baseline_window_factor",
    "allow_overlap","window_size","sn_ratio"
  )
  bad <- setdiff(names(options), allow)
  if (length(bad)) stop(paste0("unrecognized option(s): ", paste(bad, collapse = ", ")))
  num <- function(k, d) { v <- options[[k]]; if (is.null(v) || is.na(v[1])) d else as.numeric(v[1]) }
  int <- function(k, d) { v <- options[[k]]; if (is.null(v) || is.na(v[1])) d else as.integer(v[1]) }
  lgc <- function(k, d=FALSE) {
    v <- options[[k]]
    if (is.null(v)) return(as.logical(d))
    if (!is.logical(v) || length(v) != 1 || is.na(v)) stop(paste0("option '", k, "' must be logical TRUE/FALSE"))
    v
  }
  raw64 <- raw(64)
  con <- rawConnection(raw64, "r+"); on.exit(close(con), add=TRUE)
  w <- function(val, off, what=c("double","int")) {
    seek(con, off)
    if (what[1] == "double") writeBin(as.double(val), con, size=8, endian="little")
    else writeBin(as.integer(val), con, size=4, endian="little")
  }
  w(num("integral_threshold",        NaN),  0,  "double")
  w(num("intensity_threshold",       NaN),  8,  "double")
  w(int("width_threshold",           0L ), 16, "int")
  w(num("noise",                     NaN), 24, "double")
  w(if (lgc("auto_noise",          FALSE)) 1L else 0L, 32, "int")
  w(if (lgc("auto_baseline",       FALSE)) 1L else 0L, 36, "int")
  w(int("baseline_window",          0L ), 40, "int")
  w(int("baseline_window_factor",   0L ), 44, "int")
  w(if (lgc("allow_overlap",       FALSE)) 1L else 0L, 48, "int")
  w(int("window_size",              0L ), 52, "int")
  w(num("sn_ratio",                 NaN), 56, "double")
  rawConnectionValue(con)
}

dispose <- function(bin) {
  if (typeof(bin) != "externalptr") stop("msutils: expected an external pointer")
  .Call("C_dispose_mzml", bin, PACKAGE = "msutils")
}

parse_mzml <- function(data) {
  stopifnot(is.raw(data))
  .Call("C_parse_mzml", data, PACKAGE="msutils")
}

mzml_to_bin <- function(bin, level = 12L, f32_compress = FALSE) {
  if (typeof(bin) != "externalptr") stop("msutils: expected an external pointer (MzML handle)")

  if (!is.numeric(level) || length(level) != 1 || is.na(level))
    stop("`level` must be a single number 0..22")
  lvl <- as.integer(level)
  if (lvl < 0L || lvl > 22L)
    stop("`level` must be between 0 and 22 (inclusive)")

  if (!is.logical(f32_compress) || length(f32_compress) != 1 || is.na(f32_compress))
    stop("`f32_compress` must be TRUE/FALSE")

  .Call("C_mzml_to_bin", bin, lvl, f32_compress, PACKAGE = "msutils")
}

bin_to_json <- function(bin) {
  if (typeof(bin) != "externalptr") stop("msutils: expected an external pointer (MzML handle)")
  .Call("C_bin_to_json", bin, PACKAGE="msutils")
}

bin_to_mzml <- function(bin) {
  if (typeof(bin) != "externalptr") stop("msutils: expected an external pointer (MzML handle)")
  .Call("C_bin_to_mzml", bin, PACKAGE = "msutils")
}

get_peak <- function(
  x, y, rt, range,
  integral_threshold=NaN, intensity_threshold=NaN, width_threshold=0L,
  noise=NaN, auto_noise=FALSE, auto_baseline=FALSE,
  baseline_window=0L, baseline_window_factor=0L,
  allow_overlap=FALSE, window_size=0L, sn_ratio=NaN
) {
  stopifnot(is.numeric(x), is.numeric(y))
  if (length(x) != length(y) || length(x) < 3) stop("x and y must have the same length (>= 3)")
  if (!is.logical(auto_noise) || length(auto_noise) != 1 || is.na(auto_noise)) stop("auto_noise must be logical TRUE/FALSE")
  if (!is.logical(allow_overlap) || length(allow_overlap) != 1 || is.na(allow_overlap)) stop("allow_overlap must be logical TRUE/FALSE")
  if (!is.logical(auto_baseline) || length(auto_baseline) != 1 || is.na(auto_baseline)) stop("auto_baseline must be logical TRUE/FALSE")
  opt <- .pack_opts(list(
    integral_threshold=integral_threshold, intensity_threshold=intensity_threshold,
    width_threshold=width_threshold, noise=noise,
    auto_noise=auto_noise, auto_baseline=auto_baseline,
    baseline_window=baseline_window, baseline_window_factor=baseline_window_factor,
    allow_overlap=allow_overlap, window_size=window_size, sn_ratio=sn_ratio
  ))
  out_json <- .Call("C_get_peak", as.numeric(x), as.numeric(y), as.numeric(rt), as.numeric(range), opt, PACKAGE="msutils")
  jsonlite::fromJSON(out_json, simplifyVector=TRUE)
}

get_peaks_from_eic <- function(
  bin, df, from=0.5, to=5, cores=1L,
  integral_threshold=NaN, intensity_threshold=NaN, width_threshold=0L,
  noise=NaN, auto_noise=FALSE, auto_baseline=FALSE,
  baseline_window=0L, baseline_window_factor=0L,
  allow_overlap=FALSE, window_size=0L, sn_ratio=NaN
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
  opt <- .pack_opts(list(
    integral_threshold=integral_threshold, intensity_threshold=intensity_threshold,
    width_threshold=width_threshold, noise=noise,
    auto_noise=auto_noise, auto_baseline=auto_baseline,
    baseline_window=baseline_window, baseline_window_factor=baseline_window_factor,
    allow_overlap=allow_overlap, window_size=window_size, sn_ratio=sn_ratio
  ))
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
  integral_threshold=NaN, intensity_threshold=NaN, width_threshold=0L,
  noise=NaN, auto_noise=FALSE, auto_baseline=FALSE,
  baseline_window=0L, baseline_window_factor=0L,
  allow_overlap=FALSE, window_size=0L, sn_ratio=NaN
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
  opt <- .pack_opts(list(
    integral_threshold=integral_threshold, intensity_threshold=intensity_threshold,
    width_threshold=width_threshold, noise=noise,
    auto_noise=auto_noise, auto_baseline=auto_baseline,
    baseline_window=baseline_window, baseline_window_factor=baseline_window_factor,
    allow_overlap=allow_overlap, window_size=window_size, sn_ratio=sn_ratio
  ))
  cores <- .validate_cores(cores)
  out_json <- .Call("C_get_peaks_from_chrom",
    bin, idxs, rts, wins, opt, as.integer(cores), PACKAGE="msutils"
  )
  df <- jsonlite::fromJSON(out_json, simplifyVector=TRUE)
  if (!is.data.frame(df)) df <- as.data.frame(df)
  if (!"ort" %in% names(df)) stop("internal error: 'ort' missing from result")
  want <- c("index","id","ort","rt","from","to","intensity","integral", "total_area")
  present <- want[want %in% names(df)]; extras <- setdiff(names(df), present)
  df <- df[, c(present, extras), drop=FALSE]
  if ("index" %in% names(df)) df <- df[order(df$index), , drop=FALSE]
  rownames(df) <- NULL
  df
}

calculate_eic <- function(bin, targets, from, to, ppm_tolerance=20, mz_tolerance=0.005) {
  if (typeof(bin) != "externalptr") stop("msutils: expected an external pointer (MzML handle)")
  stopifnot(is.numeric(targets), length(targets) == 1)
  .Call("C_calculate_eic",
    bin, as.numeric(targets), as.numeric(from), as.numeric(to),
    as.numeric(ppm_tolerance), as.numeric(mz_tolerance),
    PACKAGE="msutils"
  )
}

find_peaks <- function(
  x, y,
  integral_threshold=NaN, intensity_threshold=NaN, width_threshold=0L,
  noise=NaN, auto_noise=FALSE, auto_baseline=FALSE,
  baseline_window=0L, baseline_window_factor=0L,
  allow_overlap=FALSE, window_size=0L, sn_ratio=NaN
) {
  stopifnot(is.numeric(x), is.numeric(y))
  if (length(x) != length(y) || length(x) < 3) stop("x and y must have the same length >= 3")
  if (!is.logical(auto_noise) || length(auto_noise) != 1 || is.na(auto_noise)) stop("auto_noise must be logical TRUE/FALSE")
  if (!is.logical(allow_overlap) || length(allow_overlap) != 1 || is.na(allow_overlap)) stop("allow_overlap must be logical TRUE/FALSE")
  if (!is.logical(auto_baseline) || length(auto_baseline) != 1 || is.na(auto_baseline)) stop("auto_baseline must be logical TRUE/FALSE")
  opt <- .pack_opts(list(
    integral_threshold=integral_threshold, intensity_threshold=intensity_threshold,
    width_threshold=width_threshold, noise=noise,
    auto_noise=auto_noise, auto_baseline=auto_baseline,
    baseline_window=baseline_window, baseline_window_factor=baseline_window_factor,
    allow_overlap=allow_overlap, window_size=window_size, sn_ratio=sn_ratio
  ))
  out_json <- .Call("C_find_peaks", as.numeric(x), as.numeric(y), opt, PACKAGE="msutils")
  jsonlite::fromJSON(out_json, simplifyVector=TRUE)
}

calculate_baseline <- function(y, baseline_window=15L, baseline_window_factor=1L) {
  stopifnot(is.numeric(y))
  .Call("C_calculate_baseline",
        as.numeric(y),
        as.integer(baseline_window),
        as.integer(baseline_window_factor),
        PACKAGE="msutils")
}

find_feature <- function(
  bin,
  rt, mz, window, id = NULL,
  scan_ppm, scan_mz, eic_ppm, eic_mz,
  cores = 1L,
  integral_threshold = NaN, intensity_threshold = NaN, width_threshold = 0L,
  noise = NaN, auto_noise = FALSE, auto_baseline = FALSE,
  baseline_window = 0L, baseline_window_factor = 0L,
  allow_overlap = FALSE, window_size = 0L, sn_ratio = NaN
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

  # FIX: Always pass a character vector so C GetHandle/R_alloc logic works
  if (is.null(id)) {
    id <- rep("", length(rt))
  } else {
    id <- as.character(id)
  }

  opt <- .pack_opts(list(
    integral_threshold = integral_threshold,
    intensity_threshold = intensity_threshold,
    width_threshold = width_threshold,
    noise = noise,
    auto_noise = auto_noise,
    auto_baseline = auto_baseline,
    baseline_window = baseline_window,
    baseline_window_factor = baseline_window_factor,
    allow_overlap = allow_overlap,
    window_size = window_size,
    sn_ratio = sn_ratio
  ))

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
  ppm_tolerance = NaN, mz_tolerance = NaN,
  grid_start = NaN, grid_end = NaN, grid_step_ppm = 0L,
  cores = 1L,
  integral_threshold = NaN, intensity_threshold = NaN, width_threshold = 0L,
  noise = NaN, auto_noise = FALSE, auto_baseline = FALSE,
  baseline_window = 0L, baseline_window_factor = 0L,
  allow_overlap = FALSE, window_size = 0L, sn_ratio = NaN
) {
  stopifnot(typeof(data) == "externalptr")
  cores <- .validate_cores(cores)

  if (!is.logical(auto_noise) || length(auto_noise) != 1 || is.na(auto_noise)) stop("auto_noise must be logical TRUE/FALSE")
  if (!is.logical(allow_overlap) || length(allow_overlap) != 1 || is.na(allow_overlap)) stop("allow_overlap must be logical TRUE/FALSE")
  if (!is.logical(auto_baseline) || length(auto_baseline) != 1 || is.na(auto_baseline)) stop("auto_baseline must be logical TRUE/FALSE")

  opt <- .pack_opts(list(
    integral_threshold = integral_threshold,
    intensity_threshold = intensity_threshold,
    width_threshold = width_threshold,
    noise = noise,
    auto_noise = auto_noise,
    auto_baseline = auto_baseline,
    baseline_window = baseline_window,
    baseline_window_factor = baseline_window_factor,
    allow_overlap = allow_overlap,
    window_size = window_size,
    sn_ratio = sn_ratio
  ))

  out_json <- .Call(
    "C_find_features",
    data,
    as.numeric(from), as.numeric(to),
    as.numeric(ppm_tolerance), as.numeric(mz_tolerance),
    as.numeric(grid_start), as.numeric(grid_end), as.numeric(grid_step_ppm),
    opt, as.integer(cores),
    PACKAGE = "msutils"
  )

  df <- jsonlite::fromJSON(out_json, simplifyVector = TRUE)
  if (!is.data.frame(df)) df <- as.data.frame(df)
  want <- c("mz","rt","from","to","intensity","integral","ratio","np")
  present <- intersect(want, names(df))
  df <- df[, c(present, setdiff(names(df), present)), drop = FALSE]
  rownames(df) <- NULL
  df
}

parse_bin <- function(bin) {
  if (!is.raw(bin)) stop("`bin` must be a raw vector (BINZ bytes)")
  .Call("C_parse_bin", bin, PACKAGE = "msutils")
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
  if (!is.na(maxc) && cores > maxc) stop("cores must be a single number")

  cores
}

bin_to_df <- function(bin) {
  if (typeof(bin) != "externalptr") stop("msutils: expected an external pointer")
  x <- jsonlite::fromJSON(bin_to_json(bin), simplifyVector = TRUE)
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