.quantion_find_arch_dir <- function() {
  sys  <- tolower(Sys.info()[["sysname"]])
  mach <- tolower(Sys.info()[["machine"]])
  if (is.na(mach) || !nzchar(mach)) mach <- ""
  if (sys == "darwin") {
    if (grepl("arm64|aarch64", mach)) "macos-arm64" else "macos-x86_64"
  } else if (sys == "linux") {
    if (grepl("arm64|aarch64", mach)) "linux-arm64" else "linux-x86_64"
  } else if (sys == "windows") {
    if (grepl("arm64|aarch64", mach)) stop("quantion: Windows ARM64 is not supported yet")
    "windows-x86_64"
  } else {
    stop(sprintf("Unsupported OS: %s (%s)", sys, mach))
  }
}

.quantion_rust_basename <- function() {
  sys <- tolower(Sys.info()[["sysname"]])
  if (sys == "windows") "libquantion.dll"
  else if (sys == "darwin") "libquantion.dylib"
  else if (sys == "linux") "libquantion.so"
  else stop("Unsupported OS")
}

.quantion_own_version <- function(libname, pkgname) {
  path <- file.path(libname, pkgname, "DESCRIPTION")
  if (!file.exists(path)) return("")
  as.character(read.dcf(path, fields = "Version")[1, 1])
}

.quantion_artifacts_root <- function(pkgname = "quantion") {
  from_env <- Sys.getenv("QUANTION_ARTIFACTS_ROOT", "")
  if (nzchar(from_env)) return(from_env)
  here <- system.file(package = pkgname)
  if (!nzchar(here)) here <- getwd()
  repeat {
    candidate <- file.path(here, "artifacts")
    if (dir.exists(candidate)) return(candidate)
    parent <- dirname(here)
    if (identical(parent, here)) return("")
    here <- parent
  }
}

.quantion_find_library <- function(own_version = "") {
  from_env <- Sys.getenv("QUANTION_LIB", "")
  if (nzchar(from_env)) return(from_env)

  arch <- .quantion_find_arch_dir()
  name <- .quantion_rust_basename()

  root <- .quantion_artifacts_root()
  if (nzchar(root)) {
    direct <- file.path(root, arch, name)
    if (file.exists(direct)) return(direct)

    if (nzchar(own_version)) {
      preferred <- file.path(root, own_version, arch, name)
      if (file.exists(preferred)) return(preferred)
    }

    found <- list.dirs(root, recursive = FALSE, full.names = FALSE)
    found <- found[vapply(found, function(v)
      file.exists(file.path(root, v, arch, name)) &&
      file.exists(file.path(root, v, "manifest.json")), logical(1))]
    if (length(found)) {
      best <- found[order(numeric_version(found, strict = FALSE), decreasing = TRUE, na.last = TRUE)][1]
      return(file.path(root, best, arch, name))
    }
  }

  bundled <- system.file("libs", arch, name, package = "quantion")
  if (nzchar(bundled) && file.exists(bundled)) return(bundled)
  ""
}

.onLoad <- function(libname, pkgname) {
  library.dynam("quantion", pkgname, libname)
  rust_path <- .quantion_find_library(.quantion_own_version(libname, pkgname))
  if (.Platform$OS.type == "windows" && nzchar(rust_path)) {
    rust_path <- normalizePath(rust_path, winslash = "\\", mustWork = FALSE)
  }
  if (!nzchar(rust_path) || !file.exists(rust_path)) {
    stop(sprintf(paste0("quantion: no library for '%s'. Build one with 'make %s', ",
                        "or set QUANTION_LIB to a library file, ",
                        "or QUANTION_ARTIFACTS_ROOT to the artifacts folder."),
                 .quantion_find_arch_dir(), .quantion_find_arch_dir()))
  }
  dyn.load(rust_path, local = FALSE, now = TRUE)
  .Call("C_bind_rust", rust_path, PACKAGE = "quantion")
  options(quantion.rust_dll_path = rust_path)
}

.onUnload <- function(libpath) {
  gc()
  try(.Call("C_unbind_rust", NULL, PACKAGE = "quantion"), silent = TRUE)
  p <- getOption("quantion.rust_dll_path", NULL)
  if (!is.null(p) && nzchar(p) && file.exists(p)) try(suppressWarnings(dyn.unload(p)), silent = TRUE)
  options(quantion.rust_dll_path = NULL)
}
