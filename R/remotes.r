#' Clone or update the IMGTHLA repo on github
#'
#' @note Due to some limitation in the \pkg{git2r} this will not work
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
  url <- "https://github.com/jrob119/IMGTHLA"
  git2r::clone(url, local_path)
}

#' @rdname fetch_IMGTHLA
#' @return \code{invisible(NULL)}
#' @export
update_IMGTHLA <- function() {
  assertive.base::assert_all_are_true(
    requireNamespace("git2r", quietly = TRUE)
  )
  repo <- git2r::repository(file.path(getOption("hlatools.local_repos"), "IMGTHLA"))
  git2r::pull(repo)
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






