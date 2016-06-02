#' @importFrom XML xpathApply xmlValue xmlName xmlGetAttr
NULL

valid_hla_loci_ <- function() {
  c('HLA-A', 'HLA-B', 'HLA-C', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DRB')
}

expand_allele <- function(x, locus = NULL) {
  if (is.null(locus)) {
    ifelse(!grepl("^HLA-\\S+", x), paste0("HLA-", x), x)
  } else {
    locus <- sub("HLA-", "", toupper(locus))
    pattern1 <- paste0("^HLA-", locus, "[*]\\d\\d\\d?:.+$")
    pattern2 <- paste0("^", locus, "[*]\\d\\d\\d?:.+$")
    pattern3 <- "^\\d\\d\\d?:.+$"
    ifelse(grepl(pattern1, x),
           x,
           ifelse(grepl(pattern2, x),
                  paste0("HLA-", x),
                  ifelse(grepl(pattern3, x),
                         paste0("HLA-", locus, "*", x), x)))

  }
}

match_hla_locus <- function(locus) {
  locus <- expand_allele(locus)
  match.arg(locus, valid_hla_loci_())
}

`%||%` <- function(a, b) {
  if (length(a) == 0) b else a
}

`%|na|%` <- function(a, b) {
  ifelse(is.na(a), b, a)
}

`%|ch|%` <- function(a, b) {
  ifelse(!nzchar(a), b, a)
}

comma <- function(...) paste0(..., collapse = ", ")

colon <- function(...) paste0(..., collapse = ":")

`%.%` <- function(f, g) {
  f <- match.fun(f)
  g <- match.fun(g)
  function(...) f(g(...))
}

normalise_ranges <- function(x) {
  if (is(x, "GenomicRanges")) {
    x <- ranges(x)
  }
  assertive::assert_is_any_of(x, "IRanges")
  cs <- cumsum(width(x))
  IRanges(
    start = c(1L, cs[-length(cs)] + 1),
    end   = cs,
    names = names(x)
  )
}

xval <- function(doc, path, alt = NA_character_, as = 'character', ...) {
  v <- unlist(xpathApply(doc, path, xmlValue, ...)) %||% alt
  set_mode(v, as)
}

xname <- function(doc, path, alt = NA_character_, as = 'character', ...) {
  n <- unlist(xpathApply(doc, path, xmlName, ...)) %||% alt
  set_mode(n, as)
}

xattr <- function(doc, path, name, alt = NA_character_, as = 'character', ...) {
  a <- unlist(xpathApply(doc, path, xmlGetAttr, name = name, ...)) %||% alt
  set_mode(a, as)
}

xsize <- function(doc, path, ...) {
  length(xpathApply(doc, path, ...))
}

xset <- function(doc, path, ...) {
  xpathApply(doc, path, fun = NULL, ...)
}

set_mode <- function(x, as) {
  AS <- match.fun(paste0('as.', as))
  if (!is.null(x)) AS(x) else x
}

strsplitN <- function(x, split, n, from = "start", collapse = split, ...) {
  assertive::assert_is_vector(x)
  from <- match.arg(from, c("start", "end"))
  xs <- strsplit(x, split, ...)
  end <- vapply(xs, length, 0L)
  if (from == "end") {
    end <- end + 1L
    n <- lapply(end, `-`, n)
    n <- .mapply(`[<-`, list(x = n, i = lapply(n, `<`, 0),
                             value = 0L), NULL)
  } else {
    n <- lapply(rep.int(0, length(xs)), `+`, n)
    n <- .mapply(`[<-`, list(x = n, i = Map(`>`, n, end),
                             value = end), NULL)
  }
  n <- lapply(n, sort %.% unique)
  unlist(.mapply(function(x, n) paste0(x[n], collapse = collapse),
                 list(x = xs, n = n), NULL))
}

recycle_vector <- function(x, length.out) {
  if (length(x) == length.out) {
    x
  } else {
    ans <- vector(storage.mode(x), length.out)
    ans[] <- x
    ans
  }
}

collapse <- function(...) paste0(..., collapse = "|")



