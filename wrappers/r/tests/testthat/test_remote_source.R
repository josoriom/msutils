serve_file <- function(port, path, honour_range) {
  body <- readBin(path, "raw", file.info(path)$size)
  total <- length(body)
  server <- serverSocket(port)
  on.exit(close(server), add = TRUE)
  repeat {
    con <- tryCatch(socketAccept(server, blocking = TRUE, open = "a+b", timeout = 30),
                    error = function(e) NULL)
    if (is.null(con)) break
    lines <- character(0)
    repeat {
      line <- readLines(con, n = 1, warn = FALSE)
      if (!length(line) || !nzchar(trimws(line))) break
      lines <- c(lines, trimws(line))
    }
    if (!length(lines)) {
      close(con)
      next
    }
    asked <- grep("^Range:", lines, ignore.case = TRUE, value = TRUE)
    first <- 0
    last <- total - 1
    partial <- FALSE
    if (honour_range && length(asked)) {
      numbers <- regmatches(asked[1], regexec("bytes=([0-9]+)-([0-9]*)", asked[1]))[[1]]
      if (length(numbers) == 3) {
        first <- as.numeric(numbers[2])
        last <- if (nzchar(numbers[3])) as.numeric(numbers[3]) else total - 1
        if (last > total - 1) last <- total - 1
        partial <- TRUE
      }
    }
    piece <- body[seq(first + 1, last + 1)]
    head_text <- paste0(
      if (partial) "HTTP/1.1 206 Partial Content\r\n" else "HTTP/1.1 200 OK\r\n",
      "Content-Length: ", length(piece), "\r\n",
      if (partial) paste0("Content-Range: bytes ", first, "-", last, "/", total, "\r\n") else "",
      "Accept-Ranges: bytes\r\n",
      "Connection: close\r\n\r\n"
    )
    writeBin(charToRaw(head_text), con)
    writeBin(piece, con)
    flush(con)
    close(con)
  }
}

find_free_port <- function() {
  for (attempt in seq_len(50)) {
    port <- sample(20000:60000, 1)
    free <- tryCatch({
      socket <- serverSocket(port)
      close(socket)
      TRUE
    }, error = function(e) FALSE)
    if (free) return(port)
  }
  0
}

with_served_file <- function(path, honour_range, work) {
  port <- find_free_port()
  if (port == 0) skip("no free port for a local server")

  child <- callr::r_bg(serve_file, args = list(port, path, honour_range), supervise = TRUE)
  on.exit(child$kill(), add = TRUE)

  url <- sprintf("http://127.0.0.1:%d/sample.ion", port)
  size <- file.info(path)$size
  ready <- FALSE
  for (attempt in seq_len(60)) {
    Sys.sleep(0.1)
    ready <- tryCatch({
      suppressWarnings(quantion:::.fetch_range(url, 0, size))
      TRUE
    }, error = function(e) FALSE)
    if (ready) break
  }
  if (!ready) skip("the local server did not answer")

  work(url)
}

sample_file <- sample_ion_path()

test_that("a remote sample is read over range requests", {
  skip_if_not_installed("quantion")
  skip_if_not_installed("callr")
  skip_if(!nzchar(sample_file), "no api.ion sample next to the package")

  with_served_file(sample_file, TRUE, function(url) {
    handle <- quantion::parse_ion_remote(url)
    expect_identical(typeof(handle), "externalptr")

    source <- attr(handle, "quantion_source")
    expect_false(is.null(source))
    expect_gt(length(source$held), 1)
  })
})

test_that("the total size of a remote sample comes from the plan, not from a guess", {
  skip_if_not_installed("quantion")
  skip_if_not_installed("callr")
  skip_if(!nzchar(sample_file), "no api.ion sample next to the package")

  with_served_file(sample_file, TRUE, function(url) {
    handle <- quantion::parse_ion_remote(url)
    source <- attr(handle, "quantion_source")

    expect_equal(source$total, file.info(sample_file)$size)
    expect_gt(source$total, 1024)
  })
})

test_that("every block a remote sample holds sits inside the file", {
  skip_if_not_installed("quantion")
  skip_if_not_installed("callr")
  skip_if(!nzchar(sample_file), "no api.ion sample next to the package")

  with_served_file(sample_file, TRUE, function(url) {
    handle <- quantion::parse_ion_remote(url)
    source <- attr(handle, "quantion_source")
    blocks <- do.call(rbind, source$held)

    expect_true(all(blocks[, 1] >= 0))
    expect_true(all(blocks[, 2] > 0))
    expect_true(all(blocks[, 1] + blocks[, 2] <= source$total))
  })
})

test_that("an extracted ion chromatogram brings in the bytes it still needs", {
  skip_if_not_installed("quantion")
  skip_if_not_installed("callr")
  skip_if(!nzchar(sample_file), "no api.ion sample next to the package")

  with_served_file(sample_file, TRUE, function(url) {
    handle <- quantion::parse_ion_remote(url)
    source <- attr(handle, "quantion_source")
    held_before <- length(source$held)

    eic <- quantion::calculate_eic(handle, 89.04768, 0.5, 1.5)

    expect_named(eic, c("x", "y"))
    expect_gt(length(eic$x), 0)
    expect_gt(length(source$held), held_before)
  })
})

test_that("a server that ignores range requests is reported", {
  skip_if_not_installed("quantion")
  skip_if_not_installed("callr")
  skip_if(!nzchar(sample_file), "no api.ion sample next to the package")

  with_served_file(sample_file, FALSE, function(url) {
    expect_error(quantion::parse_ion_remote(url), "does not support range requests")
  })
})

test_that("a remote read needs an http address", {
  skip_if_not_installed("quantion")

  expect_error(quantion::parse_ion_remote("ftp://host/sample.ion"), "http")
  expect_error(quantion::parse_ion_remote(""), "single non-empty string")
})
