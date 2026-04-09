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
    if (grepl("arm64|aarch64", cpu)) stop_install("msutils: Windows ARM64 is not supported yet")
    if (grepl("x86_64|amd64|x86-64", cpu)) return("windows-x86_64")
  }
  stop_install(sprintf("msutils: unsupported platform: %s (%s)", os, cpu))
}

copy_files <- function(source, target, platform) {
  files <- list.files(source, full.names = TRUE)
  if (!length(files)) stop_install("msutils: no files found in ", source)
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(files, target, overwrite = TRUE)
  if (!all(ok)) stop_install("msutils: failed to copy files into ", target)
  message("msutils: staged ", platform, " into ", target)
}

find_local_files <- function(platform) {
  type <- read_value("RemoteType")
  path <- read_value("RemoteUrl")
  if (type != "local" || !nzchar(path)) return("")
  roots <- c(path, dirname(path), dirname(dirname(path)))
  for (root in roots) {
    source <- file.path(root, "artifacts", platform)
    if (dir.exists(source)) return(source)
  }
  ""
}

download_files <- function(platform, target) {
  sha <- read_value("RemoteSha")
  user <- read_value("RemoteUsername")
  repo <- read_value("RemoteRepo")
  if (!nzchar(sha) || !nzchar(user) || !nzchar(repo)) {
    stop_install("msutils: missing GitHub metadata and no local artifacts found")
  }
  temp <- tempfile("msutils-")
  dir.create(temp)
  on.exit(unlink(temp, recursive = TRUE), add = TRUE)
  archive <- file.path(temp, "repo.tar.gz")
  url <- sprintf("https://codeload.github.com/%s/%s/tar.gz/%s", user, repo, sha)
  status <- tryCatch(utils::download.file(url, archive, mode = "wb", quiet = TRUE), error = function(e) 1L)
  if (!identical(status, 0L) || !file.exists(archive)) stop_install("msutils: failed to download ", url)
  untar_ok <- tryCatch({
    utils::untar(archive, exdir = temp)
    TRUE
  }, error = function(e) FALSE)
  if (!untar_ok) stop_install("msutils: failed to unpack ", archive)
  roots <- list.dirs(temp, recursive = FALSE, full.names = TRUE)
  source <- ""
  for (root in roots) {
    path <- file.path(root, "artifacts", platform)
    if (dir.exists(path)) {
      source <- path
      break
    }
  }
  if (!nzchar(source)) stop_install("msutils: artifacts/", platform, " not found for ", sha)
  copy_files(source, target, platform)
}

platform <- find_platform()
target <- file.path("inst", "libs", platform)
env_root <- Sys.getenv("MSUTILS_ARTIFACTS_ROOT", "")

if (nzchar(env_root)) {
  source <- file.path(env_root, platform)
  if (!dir.exists(source)) stop_install("msutils: artifacts not found at ", source)
  copy_files(source, target, platform)
  quit(status = 0L)
}

source <- file.path("..", "..", "artifacts", platform)
if (dir.exists(source)) {
  copy_files(source, target, platform)
  quit(status = 0L)
}

source <- find_local_files(platform)
if (nzchar(source)) {
  copy_files(source, target, platform)
  quit(status = 0L)
}

download_files(platform, target)