context("Remote data access")

test_that("fetch_IMGTHLA works", {
  skip_if(no_remote_tests)

  repo <- fetch_IMGTHLA(tmp)
  expect_is(repo, "git_repository")
  expect_equal(repo@path, file.path(tmp, "IMGTHLA"))
})

test_that("update_IMGTHLA works", {
  skip_if(no_remote_tests)

  repo <- update_IMGTHLA(tmp)
  expect_is(repo, "git_repository")
  expect_equal(repo@path, normalizePath(file.path(tmp, "IMGTHLA")))
})

test_that("update_IMGTHLA works", {
  skip_if(no_local_repo)

  repo <- update_IMGTHLA()
  expect_is(repo, "git_repository")
  expect_equal(repo@path, normalizePath(file.path(getOption("hlatools.local_repos"), "IMGTHLA")))
})

test_that("check_IMGTHLA works", {
  wrong_path <- tempdir()
  expect_message(check_IMGTHLA(wrong_path), "Use fetch_IMGTHLA() to create a local IPD-IMGT/HLA repository", fixed = TRUE)

  skip_if(no_remote_tests)
  expect_message(check_IMGTHLA(tmp), "Your local IPD-IMGT/HLA repository is up-to-date", fixed = TRUE)
})

test_that("check_IMGTHLA works", {
  skip_if(no_local_repo)
  expect_message(check_IMGTHLA(), "Your local IPD-IMGT/HLA repository is up-to-date", fixed = TRUE)
})



