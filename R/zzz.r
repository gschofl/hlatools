.onLoad <- function(libname, pkgname) {
  op <- options()
  op.hlatools <- list(
    hlatools.local_repos = "~/local/db"
  )
  toset <- !(names(op.hlatools) %in% names(op))
  if (any(toset)) {
    options(op.hlatools[toset])
  }
  ##
  check_IMGTHLA()
  ##
  invisible(NULL)
}

