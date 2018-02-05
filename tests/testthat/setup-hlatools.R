# Directory to hold the IMGTHLA repo if hlatools.check_remotes == TRUE
tmp <- tempdir()

# Test remotes or use local repo
remote_tests <- identical(getOption("hlatools.check_remotes"), TRUE)

# Check if a local repo is available for testing
local_path <- normalizePath(file.path(getOption("hlatools.local_repos"), "IMGTHLA"), mustWork = FALSE)
local_repo <- dir.exists(local_path) && !is.null(tryCatch(git2r::repository(local_path), error = function(e) NULL))
