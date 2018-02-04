# Directory to hold the IMGTHLA repo if hlatools.check_remotes == TRUE
tmp <- tempdir()

# Test remotes or use local repo
no_remote_tests <- identical(getOption("hlatools.check_remotes"), FALSE)

# Check if a local repo is available for testing
local_repo    <- file.path(getOption("hlatools.local_repos"), "IMGTHLA")
no_local_repo <- !dir.exists(local_repo) || is.null(tryCatch(git2r::repository(local_repo), error = function(e) NULL))
