#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom dplyr %>%
#' @usage lhs \%>\% rhs
NULL


valid_hla_loci_ <- function() {
  c(
    'HLA-A',    'HLA-B',    'HLA-C',
    'HLA-DPA1', 'HLA-DPB1',
    'HLA-DQA1', 'HLA-DQB1',
    'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5',
    'HLA-E',    'HLA-F',    'HLA-G',
    'HLA-DMA',  'HLA-DMB',
    'HLA-DOA',  'HLA-DOB',
    'HLA-DRA'
    # 'HLA-DRB2'
    # 'HLA-H',    'HLA-HFE',
    # 'HLA-J',    'HLA-K',    'HLA-L',
    # 'MICA',     'MICB',
    # 'HLA-P',    'HLA-T',
    # 'TAP1',     'TAP2',
    # 'HLA-V',    'HLA-W',    'HLA-Y'
  )
}

expand_hla_allele <- function(x, locus = NULL) {
  if (is.null(locus)) {
    ifelse(!grepl("^HLA-\\S+", x), paste0("HLA-", x), x)
  } else {
    locus <- sub("HLA-", "", toupper(locus))
    pattern1 <- paste0("^HLA-", locus, "[*]\\d\\d\\d?:?.*$")
    pattern2 <- paste0("^", locus, "[*]\\d\\d\\d?:?.*$")
    pattern3 <- "^\\d\\d\\d?:?.*$"
    ifelse(grepl(pattern1, x),
           x,
           ifelse(grepl(pattern2, x),
                  paste0("HLA-", x),
                  ifelse(grepl(pattern3, x),
                         paste0("HLA-", locus, "*", x), x)))

  }
}

## for backwards compatibility
expand_allele <- expand_hla_allele

match_hla_locus <- function(locusname) {
  locus <- expand_hla_allele(locusname)
  match.arg(locus, valid_hla_loci_())
}

#' (Partially) match allele in a HLAGene object.
#'
#' @param allele [character]; A vector of (partial) allele names.
#' @param x A [HLAAllele-class] or [HLAGene-class] object.
#' @param partially If `TRUE` match partial allele names.
#'
#' @return An [integer] vector
#' @export
#' @keywords internal
#' @examples
#' ##
match_alleles <- function(allele, x, partially = FALSE) {
  allele <- expand_hla_allele(x = allele, locus = locusname(x))

  ## first try matching complete allele names
  i <- match(allele, names(x), nomatch = NA)

  ## if there is no exact match and partially == TRUE try prefix matching
  if (partially && length(i <- i[!is.na(i)]) == 0) {
    for (a in allele) {
      i <- c(i, which(startsWith(names(x), a)))
    }
  }

  ## if there is still no match return an empty integer object
  if (length(i <- i[!is.na(i)]) == 0) {
    return(i)
  }

  i
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

collapse <- function(...) paste0(..., collapse = "|")

comma <- function(...) paste0(..., collapse = ", ")

colon <- function(...) paste0(..., collapse = ":")

slash <- function(...) paste0(..., collapse = "/")

`%.%` <- function(f, g) {
  f <- match.fun(f)
  g <- match.fun(g)
  function(...) f(g(...))
}

normalise_ranges <- function(x) {
  if (is(x, "GenomicRanges")) {
    x <- ranges(x)
  }
  assertive.types::assert_is_any_of(x, "IRanges")
  cs <- cumsum(width(x))
  IRanges::IRanges(
    start = c(1L, cs[-length(cs)] + 1),
    end   = cs,
    names = names(x)
  )
}

strsplit1 <- function(...) strsplit(...)[[1]]

strsplitN <- function(x, split, n, from = "start", collapse = split, ...) {
  assertive.properties::assert_is_vector(x)
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

#' Set remote data checking on of off
#' @param what `TRUE` or `FALSE`
#' @export
#' @keywords internal
check_remotes <- function(what) {
  if (missing(what)) {
    if (identical(getOption("hlatools.check_remotes"), TRUE)) {
      message("Checking remote data access is turned ON")
    } else if (identical(getOption("hlatools.check_remotes"), FALSE)) {
      message("Checking remote data access is turned OFF")
    } else {
      warning("Checking remote data access is not set")
    }
  } else {
    stopifnot(is.logical(what))
    options("hlatools.check_remotes" = what)
    if (identical(getOption("hlatools.check_remotes"), TRUE)) {
      message("Checking remote data access is turned ON")
    } else if (identical(getOption("hlatools.check_remotes"), FALSE)) {
      message("Checking remote data access is turned OFF")
    }
  }
  invisible(NULL)
}

