context("Remote data access")

test_that("clone_IMGTHLA works", {
  skip_if_not(remote_tests)

  repo <- clone_IMGTHLA(tmp)
  expect_is(repo, "git_repository")
  expect_equal(repo$path, file.path(tmp, "IMGTHLA", ".git"))
  expect_equal(git2r::repository_head(repo)$name, "Latest")
})

test_that("pull_IMGTHLA works", {
  skip_if_not(remote_tests)

  repo <- pull_IMGTHLA(tmp)
  expect_is(repo, "git_repository")
  expect_equal(repo$path, normalizePath(file.path(tmp, "IMGTHLA", ".git")))
})

test_that("pull_IMGTHLA works", {
  skip_if_not(local_repo)

  repo <- pull_IMGTHLA()
  expect_is(repo, "git_repository")
  expect_equal(repo$path, normalizePath(file.path(getOption("hlatools.local_repos"), "IMGTHLA", ".git")))
})

test_that("check_IMGTHLA works", {
  skip_if(remote_tests)

  wrong_path <- tempdir()
  expect_message(check_IMGTHLA(wrong_path), "Use clone_IMGTHLA() to create a local IPD-IMGT/HLA repository", fixed = TRUE)
})

test_that("check_IMGTHLA works", {
  skip_if_not(remote_tests)
  expect_message(check_IMGTHLA(tmp), "Your local IPD-IMGT/HLA repository is up-to-date", fixed = TRUE)
})

test_that("check_IMGTHLA works", {
  skip_if_not(local_repo)
  expect_message(check_IMGTHLA(), "Your local IPD-IMGT/HLA repository is up-to-date", fixed = TRUE)
})



