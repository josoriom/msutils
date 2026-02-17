.msutils_find_arch_dir <- function() {
  sys  <- tolower(Sys.info()[["sysname"]])
  mach <- tolower(Sys.info()[["machine"]])
  if (is.na(mach) || !nzchar(mach)) mach <- ""
  if (sys == "darwin") {
    if (grepl("arm64|aarch64", mach)) "macos-arm64" else "macos-x86_64"
  } else if (sys == "linux") {
    if (grepl("arm64|aarch64", mach)) "linux-arm64" else "linux-x86_64"
  } else if (sys == "windows") {
    if (grepl("arm64|aarch64", mach)) stop("msutils: Windows ARM64 is not supported yet")
    "windows-x86_64"
  } else {
    stop(sprintf("Unsupported OS: %s (%s)", sys, mach))
  }
}

.msutils_rust_basename <- function() {
  sys <- tolower(Sys.info()[["sysname"]])
  if (sys == "windows") "libmsutils.dll"
  else if (sys == "darwin") "libmsutils.dylib"
  else if (sys == "linux") "libmsutils.so"
  else stop("Unsupported OS")
}

.onLoad <- function(libname, pkgname) {
  library.dynam("msutils", pkgname, libname)
  arch <- .msutils_find_arch_dir()
  rust_path <- system.file("libs", arch, .msutils_rust_basename(), package = pkgname)
  if (.Platform$OS.type == "windows") rust_path <- normalizePath(rust_path, winslash = "\\", mustWork = FALSE)
  if (!nzchar(rust_path) || !file.exists(rust_path)) {
    have <- tryCatch(list.files(system.file("libs", package = pkgname), recursive = TRUE), error = function(e) character())
    stop(sprintf("msutils: Rust library not found for '%s'. Expected: %s\nAvailable under libs/: %s", arch, rust_path, paste(have, collapse = ", ")))
  }
  dyn.load(rust_path, local = FALSE, now = TRUE)
  print(rust_path)
  .Call("C_bind_rust", rust_path, PACKAGE = "msutils")
  options(msutils.rust_dll_path = rust_path)
}

.onUnload <- function(libpath) {
  p <- getOption("msutils.rust_dll_path", NULL)
  if (!is.null(p) && nzchar(p) && file.exists(p)) try(suppressWarnings(dyn.unload(p)), silent = TRUE)
  options(msutils.rust_dll_path = NULL)
}
