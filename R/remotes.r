#' Clone or update the IMGTHLA repo on github
#'
#' @note Due to limitations in the \pkg{git2r} this will not work
#' behind a proxy.
#' @param local_path Local directory to clone to. Defaults to \file{~/local/db/IMGTHLA}.
#'
#' @return A S4 \code{\linkS4class{git_repository}} object.
#' @export
#' @examples
#' \dontrun{
#' repo <- fetch_IMGTHLA()
#' update_IMGTHLA()
#' }
fetch_IMGTHLA <- function(local_path = getOption("hlatools.local_repos")) {
  assertive.base::assert_all_are_true(
    requireNamespace("git2r", quietly = TRUE)
  )
  local_path <- normalizePath(file.path(local_path, "IMGTHLA"), mustWork = FALSE)
  if (!dir.exists(local_path)) {
    dir.create(local_path, recursive = TRUE)
  }
  # Currently not working behind proxy
  url <- "https://github.com/ANHIG/IMGTHLA"
  git2r::clone(url, local_path)
}

#' @rdname fetch_IMGTHLA
#' @return \code{invisible(NULL)}
#' @export
update_IMGTHLA <- function() {
  assertive.base::assert_all_are_true(
    requireNamespace("git2r", quietly = TRUE)
  )
  path <- file.path(getOption("hlatools.local_repos"), "IMGTHLA")
  tryCatch({
    repo <- git2r::repository(path)
    }, error =  function(e) {
      if (grepl("Unable to open repository", e$message))
        ## go and try cloning the IMGT/HLA repo
        return(fetch_IMGTHLA())
      else
        stop(e$message, call. = FALSE)
    })
  git2r::pull(repo)
  repo
}


#' Fetch and parse or update the IMGT/HLA hla.xml file
#'
#' @param remote [logical] Pull data from the IMGT/HLA ftp server or
#' \code{getOption("hlatools.local_repos")}
#'
#' @return An object of class (S3) \code{XMLInternalDocument}.
#' @export
#' @examples \dontrun{
#' doc <- read_hla_xml(remote = TRUE)
#' update_hla_xml()
#' }
read_hla_xml <- function(remote = FALSE) {
  tdir <- tempdir()
  if (remote) {
    ftpfile <- "ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/xml/hla.xml.zip"
    dlfile <- tempfile(tmpdir = tdir)
    download.file(url = ftpfile, destfile = dlfile, method = "libcurl")
    tfile <- unzip(zipfile = dlfile, exdir = tdir)
  } else {
    dbpath <- getOption("hlatools.local_repos")
    dbfile <- normalizePath(file.path(dbpath, "IMGTHLA", "xml", "hla.xml.zip"), mustWork = TRUE)
    tfile <- unzip(zipfile = dbfile, exdir = tdir)[1]
  }
  doc <- xml2::read_xml(tfile)
  unlink(tfile, force = TRUE)
  doc
}

#' @rdname read_hla_xml
#' @return \code{invisible(NULL)}
#' @export
update_hla_xml <- function() {
  assertive.properties::assert_is_not_null(getOption("hlatools.local_repos"))
  dbpath <- getOption("hlatools.local_repos")
  dlfile <- normalizePath(file.path(dbpath, "IMGTHLA", "xml", "hla.xml.zip"), mustWork = FALSE)
  ftpfile <- "ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/xml/hla.xml.zip"
  download.file(url = ftpfile, destfile = dlfile, method = "libcurl")
}


#' @keywords internal
checkout_db_version <- function(version = "Latest") {
  versions <- c("Latest", "3.29.0", "3.28.0", "3.27.0", "3.26.0", "3.25.0", "3.24.0", "3.23.0",
                "3.22.0", "3.21.0", "3.20.0", "3.19.0", "3.18.0", "3.17.0", "3.16.0", "3.15.0",
                "3.14.0", "3.13.0", "3.12.0", "3.11.0", "3.10.0", "3.9.0",  "3.8.0",   "3.7.0",
                "3.6.0",  "3.5.0",  "3.4.0",  "3.3.0",  "3.2.0",  "3.1.0",  "3.0.0")
  version <- match.arg(version, versions)
  assertive.base::assert_all_are_true(
    requireNamespace("git2r", quietly = TRUE)
  )
  if (version != "Latest") {
    version <- gsub("\\.", "", version)
  }
  path <- file.path(getOption("hlatools.local_repos"), "IMGTHLA")
  tryCatch({
    repo <- git2r::repository(path)
  }, error =  function(e) {
    if (grepl("Unable to open repository", e$message))
      ## go and try cloning the IMGT/HLA repo
      return(fetch_IMGTHLA())
    else
      stop(e$message, call. = FALSE)
  })

  refs <- git2r::references(repo)

  ## Check if "version" matches a single reference
  if (sum(na.omit(vapply(refs, function(r) strsplit(r@shorthand, "/")[[1]][2], FUN.VALUE = "") == version)) != 1) {
    stop("Version ", version, " not found")
  }

  git2r::checkout(repo, version)

  invisible(repo)
}



