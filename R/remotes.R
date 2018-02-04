#' Clone or pull from the [ANHIG/IMGTHLA](https://github.com/ANHIG/IMGTHLA) repo
#'
#' @note Due to limitations in the \pkg{git2r} this will not work
#' behind a proxy.
#' @param local_path Local directory to clone to. Defaults to \file{~/local/db/IMGTHLA}.
#'
#' @return A [git_repository-class] object.
#' @export
#' @examples
#' \dontrun{
#' check_IMGTHLA()
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
#' @return `invisible(NULL)`
#' @export
update_IMGTHLA <- function(local_path = getOption("hlatools.local_repos")) {
  assertive.base::assert_all_are_true(
    requireNamespace("git2r", quietly = TRUE)
  )
  repo_path <- file.path(local_path, "IMGTHLA")
  tryCatch({
    repo <- git2r::repository(repo_path)
    }, error =  function(e) {
      if (grepl("Unable to open repository", e$message))
        ## go and try cloning the IPD-IMGT/HLA repo
        return(fetch_IMGTHLA())
      else
        stop(e$message, call. = FALSE)
    })
  git2r::pull(repo)
  repo
}

#' @rdname fetch_IMGTHLA
#' @return A diagnostic message.
#' @export
check_IMGTHLA <- function(local_path = getOption("hlatools.local_repos")) {
  assertive.base::assert_all_are_true(
    requireNamespace("git2r", quietly = TRUE)
  )
  ##
  repo_path <- file.path(local_path, "IMGTHLA")
  repo <- tryCatch({
    git2r::repository(repo_path)
  }, error =  function(e) {
    if (grepl("No such file or directory", e$message))
      structure("Use fetch_IMGTHLA() to create a local IPD-IMGT/HLA repository", class = "msg")
    else
      structure(e$message, class = "msg")
  })
  ##
  remote_shas <- tryCatch({
    git2r::remote_ls("https://github.com/ANHIG/IMGTHLA")
  }, error = function(e) {
    structure("Remote ANHIG/IMGTHLA repository not reachable", class = "msg")
  })
  ##
  if (is(repo, "git_repository") && is(remote_shas, "character")) {
    rl          <- git2r::reflog(repo)
    ##
    msg <- if (rl[[1]]@sha == remote_shas[["HEAD"]]) {
      "Your local IPD-IMGT/HLA repository is up-to-date"
    } else {
      "Run update_IMGTHLA() to bring your local IPD-IMGT/HLA repository up-to-date"
    }
  } else if (is(repo, "msg")) {
    msg <- repo
  } else if (is(remote_shas, "msg")) {
    msg <- remote_shas
  }

  packageStartupMessage(msg)
}

#' @keywords internal
checkout_db_version <- function(version = "Latest") {
  versions <- c("Latest", "3.31.0", "3.30.0", "3.29.0", "3.28.0", "3.27.0",
                "3.26.0", "3.25.0", "3.24.0", "3.23.0", "3.22.0", "3.21.0",
                "3.20.0", "3.19.0", "3.18.0", "3.17.0", "3.16.0", "3.15.0",
                "3.14.0", "3.13.0", "3.12.0", "3.11.0", "3.10.0", "3.9.0",
                "3.8.0",  "3.7.0",  "3.6.0",  "3.5.0",  "3.4.0",  "3.3.0",
                "3.2.0",  "3.1.0",  "3.0.0")
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
      ## go and try cloning the IPD-IMGT/HLA repo
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



