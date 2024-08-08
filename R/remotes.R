#' Clone or pull the [ANHIG/IMGTHLA](https://github.com/ANHIG/IMGTHLA) repository
#'
#' @param db_path Local directory to clone into. Defaults to \file{~/local/db/IMGTHLA}.
#' @param branch The name of the branch to checkout. Defaults to the remote's default branch.
#' @param progress Show progress. Default is `TRUE`.
#'
#' @return A [git_repository-class] object.
#' @export
#' @examples
#' \dontrun{
#' check_IMGTHLA()
#' repo <- clone_IMGTHLA()
#' repo <- pull_IMGTHLA()
#' }
clone_IMGTHLA <- function(db_path = getOption("hlatools.local_repos"),
                          branch = NULL,
                          progress = TRUE) {
  assertive.base::assert_all_are_true(
    requireNamespace("git2r", quietly = TRUE)
  )
  repo_path <- normalizePath(file.path(db_path, "IMGTHLA"), mustWork = FALSE)
  if (!dir.exists(repo_path)) {
    dir.create(repo_path, recursive = TRUE)
  }
  repo <- git2r::clone("https://github.com/ANHIG/IMGTHLA", repo_path,
                       branch = branch, progress = progress)
  invisible(repo)
}

#' @export
fetch_IMGTHLA <- function(db_path = getOption("hlatools.local_repos"),
                          branch = NULL,
                          progress = TRUE) {
  .Deprecated("clone_IMGTHLA")
  clone_IMGTHLA(db_path, branch, progress)
}

#' @rdname clone_IMGTHLA
#' @export
pull_IMGTHLA <- function(db_path = getOption("hlatools.local_repos")) {
  assertive.base::assert_all_are_true(
    requireNamespace("git2r", quietly = TRUE)
  )
  repo_path <- normalizePath(file.path(db_path, "IMGTHLA"), mustWork = FALSE)
  tryCatch({
    repo <- git2r::repository(repo_path)
    }, error =  function(e) {
      if (grepl("Unable to open repository", e$message))
        ## go and try cloning the IPD-IMGT/HLA repo
        return(clone_IMGTHLA())
      else
        stop(e$message, call. = FALSE)
    })
  git2r::pull(repo)
  invisible(repo)
}

#' @export
update_IMGTHLA <- function(db_path = getOption("hlatools.local_repos")) {
  .Deprecated("pull_IMGTHLA")
  pull_IMGTHLA(db_path)
}

#' @rdname clone_IMGTHLA
#' @return [character]; a diagnostic message.
#' @export
check_IMGTHLA <- function(db_path = getOption("hlatools.local_repos")) {
  assertive.base::assert_all_are_true(
    requireNamespace("git2r", quietly = TRUE)
  )
  ##
  repo_path <- normalizePath(file.path(db_path, "IMGTHLA"), mustWork = FALSE)
  repo <- tryCatch({
    git2r::repository(repo_path)
  }, error =  function(e) {
    if (grepl("No such file or directory", e$message))
      structure("Use clone_IMGTHLA() to create a local IPD-IMGT/HLA repository", class = "msg")
    else
      structure(e$message, class = "msg")
  }, warning =  function(e) {
    if (grepl("No such file or directory", e$message))
      structure("Use clone_IMGTHLA() to create a local IPD-IMGT/HLA repository", class = "msg")
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
    rl <- git2r::reflog(repo)
    ##
    msg <- if (utils::packageVersion("git2r") < "0.22.1") {
      if (rl[[1]]@sha == remote_shas[["HEAD"]]) {
        "Your local IPD-IMGT/HLA repository is up-to-date"
      } else {
        "Run pull_IMGTHLA() to bring your local IPD-IMGT/HLA repository up-to-date"
      }
    } else {
      if (git2r::sha(rl[[1]]) == remote_shas[["HEAD"]]) {
        "Your local IPD-IMGT/HLA repository is up-to-date"
      } else {
        "Run pull_IMGTHLA() to bring your local IPD-IMGT/HLA repository up-to-date"
      }
    }
  } else if (is(repo, "msg")) {
    msg <- repo
  } else if (is(remote_shas, "msg")) {
    msg <- remote_shas
  }

  packageStartupMessage(msg)
}

#' @keywords internal
checkout_db_version <- function(db_version = "Latest",
                                db_path = getOption("hlatools.local_repos")) {
  db_versions <- c("Latest", paste("3", rev(0:56), "0", sep = "."))
  db_version <- match.arg(db_version, db_versions)
  assertive.base::assert_all_are_true(
    requireNamespace("git2r", quietly = TRUE)
  )
  if (db_version != "Latest") {
    db_version <- gsub("\\.", "", db_version)
  }
  repo_path <- normalizePath(file.path(db_path, "IMGTHLA"), mustWork = FALSE)
  tryCatch({
    repo <- git2r::repository(repo_path)
  }, error =  function(e) {
    if (grepl("Unable to open repository", e$message))
      ## go and try cloning the IPD-IMGT/HLA repo
      return(clone_IMGTHLA())
    else
      stop(e$message, call. = FALSE)
  })

  refs <- git2r::references(repo)

  ## Check if "db_version" matches a single reference
  shorthand_ <- if (utils::packageVersion("git2r") < "0.22.1") {
    function(r) strsplit1(r@shorthand, "/")[2]
  } else {
    function(r) strsplit1(r$shorthand, "/")[2]
  }
  if (sum(na.omit(vapply(refs, shorthand_, FUN.VALUE = "") == db_version)) != 1) {
    stop("Version ", db_version, " not found")
  }

  git2r::checkout(repo, db_version)

  invisible(repo)
}





