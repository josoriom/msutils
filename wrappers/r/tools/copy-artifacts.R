if (nzchar(Sys.getenv("QUANTION_ARTIFACTS_ROOT", "")) || nzchar(Sys.getenv("QUANTION_LIB", ""))) {
  message("quantion: the package will load the library at run time - nothing copied in")
  quit(status = 0L)
}

read_value <- function(name) {
  value <- tryCatch(read.dcf("DESCRIPTION")[1, name], error = function(e) "")
  if (length(value) == 0 || is.na(value) || !nzchar(value)) "" else trimws(value)
}

stop_install <- function(...) {
  message(...)
  quit(status = 1L)
}

find_platform <- function() {
  os <- tolower(Sys.info()[["sysname"]])
  cpu <- tolower(Sys.info()[["machine"]])
  if (is.na(cpu) || !nzchar(cpu)) cpu <- ""
  if (os == "darwin") {
    if (grepl("arm64|aarch64", cpu)) return("macos-arm64")
    if (grepl("x86_64|amd64|x86-64", cpu)) return("macos-x86_64")
  }
  if (os == "linux") {
    if (grepl("arm64|aarch64", cpu)) return("linux-arm64")
    if (grepl("x86_64|amd64|x86-64", cpu)) return("linux-x86_64")
  }
  if (os == "windows") {
    if (grepl("arm64|aarch64", cpu)) stop_install("quantion: Windows ARM64 is not supported yet")
    if (grepl("x86_64|amd64|x86-64", cpu)) return("windows-x86_64")
  }
  stop_install(sprintf("quantion: unsupported platform: %s (%s)", os, cpu))
}

find_platform_dir <- function(root, platform) {
  direct <- file.path(root, platform)
  if (dir.exists(direct)) return(direct)
  found <- list.dirs(root, recursive = FALSE, full.names = FALSE)
  found <- found[vapply(found, function(v) dir.exists(file.path(root, v, platform)), logical(1))]
  if (!length(found)) return("")
  best <- found[order(numeric_version(found, strict = FALSE), decreasing = TRUE, na.last = TRUE)][1]
  file.path(root, best, platform)
}

copy_files <- function(source, target, platform) {
  files <- list.files(source, full.names = TRUE)
  if (!length(files)) stop_install("quantion: no files found in ", source)
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(files, target, overwrite = TRUE)
  if (!all(ok)) stop_install("quantion: failed to copy files into ", target)
  message("quantion: staged ", platform, " into ", target)
}

find_local_files <- function(platform) {
  type <- read_value("RemoteType")
  path <- read_value("RemoteUrl")
  if (type != "local" || !nzchar(path)) return("")
  roots <- c(path, dirname(path), dirname(dirname(path)))
  for (root in roots) {
    source <- find_platform_dir(file.path(root, "artifacts"), platform)
    if (nzchar(source)) return(source)
  }
  ""
}

download_files <- function(platform, target) {
  sha <- read_value("RemoteSha")
  user <- read_value("RemoteUsername")
  repo <- read_value("RemoteRepo")
  if (!nzchar(sha) || !nzchar(user) || !nzchar(repo)) {
    stop_install("quantion: missing GitHub metadata and no local artifacts found")
  }
  temp <- tempfile("quantion-")
  dir.create(temp)
  on.exit(unlink(temp, recursive = TRUE), add = TRUE)
  archive <- file.path(temp, "repo.tar.gz")
  url <- sprintf("https://codeload.github.com/%s/%s/tar.gz/%s", user, repo, sha)
  status <- tryCatch(utils::download.file(url, archive, mode = "wb", quiet = TRUE), error = function(e) 1L)
  if (!identical(status, 0L) || !file.exists(archive)) stop_install("quantion: failed to download ", url)
  untar_ok <- tryCatch({
    utils::untar(archive, exdir = temp)
    TRUE
  }, error = function(e) FALSE)
  if (!untar_ok) stop_install("quantion: failed to unpack ", archive)
  roots <- list.dirs(temp, recursive = FALSE, full.names = TRUE)
  source <- ""
  for (root in roots) {
    path <- find_platform_dir(file.path(root, "artifacts"), platform)
    if (nzchar(path)) {
      source <- path
      break
    }
  }
  if (!nzchar(source)) stop_install("quantion: artifacts/", platform, " not found for ", sha)
  copy_files(source, target, platform)
}

repo_artifacts <- function() {
  here <- normalizePath(getwd(), mustWork = FALSE)
  repeat {
    candidate <- file.path(here, "artifacts")
    if (dir.exists(candidate)) return(candidate)
    parent <- dirname(here)
    if (identical(parent, here)) return("")
    here <- parent
  }
}

if (nzchar(repo_artifacts())) {
  message("quantion: building against ", repo_artifacts(), " - nothing copied into the package")
  quit(status = 0L)
}

platform <- find_platform()
target <- file.path("inst", "libs", platform)
env_root <- Sys.getenv("QUANTION_ARTIFACTS_ROOT", "")

if (nzchar(env_root)) {
  source <- find_platform_dir(env_root, platform)
  if (!nzchar(source)) stop_install("quantion: artifacts not found at ", env_root)
  copy_files(source, target, platform)
  quit(status = 0L)
}

source <- find_platform_dir(file.path("..", "..", "artifacts"), platform)
if (nzchar(source)) {
  copy_files(source, target, platform)
  quit(status = 0L)
}

source <- find_local_files(platform)
if (nzchar(source)) {
  copy_files(source, target, platform)
  quit(status = 0L)
}

download_files(platform, target)