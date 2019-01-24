context("HLAAlleles")

doc <- if (remote_tests) {
  read_hla_xml(db_path = tmp, remote = FALSE)
} else if (local_repo) {
  read_hla_xml(db_path = getOption("hlatools.local_repos"), remote = FALSE)
}

test_that("read_hla_xml() works", {
  expect_is(doc, "xml_document")
})

dpa <- parse_hla_alleles(doc, "HLA-DPA1")

test_that("HLAAllele global API", {
  expect_is(dpa, "HLAAllele")
  ##
  expect_error(db_version(dpa))
  expect_error(hlatools_version(dpa))
  expect_equal(locusname(dpa), "HLA-DPA1")
  ##
  ## Features
  expect_is(features(dpa), "CompressedHLARangesList")
  expect_is(features(dpa)[[1]], "HLARanges")
  expect_equal(length(features(dpa)[[1]]), 9)
  ## Sequences
  expect_is(sequences(dpa), "DNAStringSet")
  expect_equal(length(sequences(dpa)), 85)
  expect_equal(names(sequences(dpa))[[1]], "HLA-DPA1*01:03:01:01")
  expect_equal(width(sequences(dpa))[[1]], 9775)
  ## Metadata
  expect_is(elementMetadata(dpa), "DataFrame")
  expect_equal(length(elementMetadata(dpa)), 14)
})

test_that("HLAAllele element-wise API", {
  ## Subsetting
  x <- dpa["HLA-DPA1*01:03:01:01"]
  expect_equal(length(x), 1)
  ##
  expect_equal(allele_id(x), "HLA00499")
  ##
  expect_equal(allele_name(x), "HLA-DPA1*01:03:01:01")
  expect_equal(names(x), "HLA-DPA1*01:03:01:01")
  ##
  expect_equal(g_group(x), "DPA1*01:03:01G")
  ##
  expect_equal(p_group(x), "DPA1*01:03P")
  ##
  expect_equal(cwd_status(x), "Common")
  ##
  expect_equal(ethnicity(x), "Caucasoid:Oriental")
  ##
  expect_equal(sample_name(x), "BOLETH:JM15:LB:LG2:PRIESS:QBL:T5-1:TUBO:WJR076")
  ##
  expect_true(is_complete(x))
  expect_false(is_lsl(x))
  ##
  expect_is(as(x, "data.table"), "tbl_dt")
})

test_that("HLAAllele subsetting", {
  ## Subsetting by allele name
  expect_equal(NROW(dpa["HLA-DPA1*01:03:01:01"]), 1)
  expect_equal(NROW(dpa["DPA1*01:03:01:01"]), 1)
  expect_equal(NROW(dpa["01:03:01:01"]), 1)
  ## Partial subsetting by allele name
  expect_equal(NROW(dpa["01:03:01"]), 19)
  expect_equal(NROW(dpa["DPA1*01:03"]), 25)
  expect_equal(NROW(dpa["HLA-DPA1*01"]), 40)
  ## Subsetting by number
  expect_equal(NROW(dpa[1]), 1)
  expect_equal(NROW(dpa[1:4]), 4)
  ## Combine
  x <- dpa[1:2]
  y <- dpa[3:4]
  expect_equal(NROW(c(x, y)), 4)
})

test_that("HLAAllele exon, intron, and utr", {
  x <- dpa["01:03:01:01"]
  ##
  e1 <- exon(x, 1)
  expect_is(e1, "DNAStringSet")
  expect_equal(width(e1), 100)
  expect_equal(Biostrings::toString(Biostrings::subseq(e1, 1, 10)), "ATGCGCCCTG")
  ##
  e2 <- exon(x, 2)
  expect_equal(Biostrings::toString(Biostrings::subseq(e2, 1, 10)), "CGGACCATGT")
  ##
  i12 <- intron(x, 1:2)
  expect_is(i12, "DNAStringSet")
  ##
  u <- utr(x, 2)
  expect_equal(Biostrings::toString(Biostrings::subseq(u, 1, 10)), "AATACTGTAA")
  expect_error(utr(x, 3))
})
