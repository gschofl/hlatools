#' Clone the IMGTHLA repo on github
#'
#' @param local_path Local directory to clone to. Defaults to \file{~/local/db/IMGTHLA}.
#'
#' @return A S4 \code{\linkS4class{git_repository}} object.
#' @export
#' @examples
#' \dontrun{
#' repo <- fetch_IMGTHLA()
#' }
fetch_IMGTHLA <- function(local_path = getOption("hlatools.local_repos")) {
  stopifnot(requireNamespace("git2r", quietly = TRUE))
  local_path <- normalizePath(file.path(local_path, "IMGTHLA"), mustWork = FALSE)
  if (!dir.exists(local_path)) {
    dir.create(local_path, recursive = TRUE)
  }
  # Currently not working behind proxy
  url <- "https://github.com/jrob119/IMGTHLA"
  git2r::clone(url, local_path)
}


#' Update a local IMGTHLA repo
#'
#' @return \code{invisible(NULL)}
#' @export
#' @examples
#' \dontrun{
#' update_IMGTHLA()
#' }
update_IMGTHLA <- function() {
  stopifnot(requireNamespace("git2r", quietly = TRUE))
  repo <- git2r::repository(file.path(getOption("hlatools.local_repos"), "IMGTHLA"))
  git2r::pull(repo)
}


#' Parse the IMGT/HLA hla.xml file
#'
#' @param remote [logical] Pull data from the IMGT/HLA ftp server or
#' \code{getOption("hlatools.local_repos")}
#'
#' @return An object of class (S3) \code{XMLInternalDocument}.
#' @export
#' @examples
#' doc <- read_hla_xml(remote = TRUE)
read_hla_xml <- function(remote = FALSE) {
  stopifnot(requireNamespace("XML", quietly = TRUE))

  tdir <- tempdir()

  if (remote) {
    ftpfile <- "ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/xml/hla.xml.zip"
    dlfile <- tempfile(tmpdir = tdir)
    download.file(url = ftpfile, destfile = dlfile, method = "libcurl")
    tfile <- unzip(zipfile = dlfile, exdir = tdir)
  } else {
    dbpath <- getOption("hlatools.local_repos")
    dbfile <- normalizePath(file.path(dbpath, "IMGTHLA", "xml", "hla.xml.zip"), mustWork = TRUE)
    tfile <- unzip(zipfile = dbfile, exdir = tdir)
  }

  doc <- XML::xmlInternalTreeParse(tfile)
  unlink(tfile, force = TRUE)
  doc
}

