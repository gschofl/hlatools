context("HLAGene")

dpa <- if (remote_tests) {
  HLAGene(locusname = "DPA1", db_version = "3.30.0", db_path = tmp)
} else if (local_repo) {
  HLAGene(locusname = "DPA1", db_version = "3.30.0", db_path = getOption("hlatools.local_repos"))
}

test_that("HLAGene global API", {
  expect_is(dpa, "HLAGene")
  ##
  expect_equal(db_version(dpa), numeric_version("3.30.0"))
  expect_equal(hlatools_version(dpa), utils::packageVersion("hlatools"))
  expect_equal(locusname(dpa), "HLA-DPA1")
  ##
  expect_is(features(dpa), "CompressedHLARangesList")
  expect_is(elementMetadata(dpa), "DataFrame")
  expect_is(sequences(dpa), "DNAStringSet")
})

test_that("HLAGene element-wise API", {
  ## subsetting
  x <- dpa["01:03:01:01"]
  expect_is(x, "HLAAllele")
  ##
  expect_length(allele_id(dpa), 56)
  expect_equal(allele_id(dpa)[1], "HLA00499", fixed = TRUE)
  ##
  expect_length(allele_name(dpa), 56)
  expect_equal(allele_name(dpa)[1], "HLA-DPA1*01:03:01:01", fixed = TRUE)
  expect_equal(names(dpa)[1], "HLA-DPA1*01:03:01:01", fixed = TRUE)
  ##
  expect_length(g_group(dpa), 56)
  expect_equal(g_group(dpa)[1], "DPA1*01:03:01G", fixed = TRUE)
  ##
  expect_length(p_group(dpa), 56)
  expect_equal(p_group(dpa)[1], "DPA1*01:03P", fixed = TRUE)
  ##
  expect_length(cwd_status(dpa), 56)
  expect_equal(cwd_status(dpa)[1], "Common", fixed = TRUE)
  ##
  expect_length(ethnicity(dpa), 56)
  expect_equal(ethnicity(dpa)[1], "Caucasoid:Oriental", fixed = TRUE)
  ##
  expect_length(sample_name(dpa), 56)
  expect_equal(sample_name(dpa)[1], "BOLETH:JM15:LB:LG2:PRIESS:QBL:T5-1:TUBO", fixed = TRUE)
  ##
  expect_true(all(is.logical(is_complete(dpa))))
  expect_equal(sum(is_complete(dpa)), 22)
  ##
  expect_true(all(is.logical(is_lsl(dpa))))
  expect_equal(sum(is_lsl(dpa)), 0)
  ## Partial subsetting
  y <- dpa["01:03"]
  expect_equal(length(y), 12)
})






