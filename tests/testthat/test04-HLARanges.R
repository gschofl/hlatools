context("HLARanges")

doc <- read_hla_xml(remote = FALSE)
dpa <- parse_hla_alleles(doc, "HLA-DPA1")
string0 <- "HLA-DPA1*01:03:01:01.1~1:UTR5:1:523:499.1:U:::|1:Exon1:524:623:499.2:E:C:1:|1:Intron1:624:4207:499.3:I:::|1:Exon2:4208:4453:499.4:E:C:3:|1:Intron2:4454:4793:499.5:I:::|1:Exon3:4794:5075:499.6:E:C:3:|1:Intron3:5076:5289:499.7:I:::|1:Exon4:5290:5444:499.8:E:C:3:|1:UTR3:5445:9775:499.9:U:::"
string1 <- "FAKE*01.1~1:Exon1:524:623:499.2:E:C:1:|1:Intron1:624:4207:499.3:I:::"

test_that("HLARanges methods work correctly", {
  ## Subsetting
  x <- dpa["01:03:01:01"]
  f <- features(x)[[1]]
  ##
  expect_length(getId(f), 9)
  expect_equal(names(f)[c(1, 9)], c("5' UTR", "3' UTR"))
  ## ranges
  r <- ranges(f)
  expect_is(r, "IRanges")
  expect_equal(start(r)[c(1, 9)], c(1, 5445))
  expect_equal(end(r)[c(1, 9)], c(523, 9775))
  expect_equal(width(r)[c(1, 9)], c(523, 4331))
})

test_that("HLARanges to string works", {
  ## Subsetting
  x <- dpa["01:03:01:01"]
  f <- features(x)[[1]]
  ##
  expect_equal(as(f, "character"), string0)
  ##
  expect_is(as(string1, "HLARanges"), "HLARanges")
})

