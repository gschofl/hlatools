#' @include HLARanges.R
NULL

#' R6 Class `"HLAAllele"``
#'
#' A container for data parsed from the IPD-IMGT/HLA hla.xml file.
#' Part of a [HLAGene][HLAGene_] object.
#'
#' @slot locusname A `character` string.
#' @slot sequence  A [DNAStringSet-class] object.
#' @slot metadata  A [DataFrame-class] object.
#' @slot features  A [HLARangesList-class] object.
#'
#' @importClassesFrom Biostrings DNAStringSet
#' @importClassesFrom S4Vectors DataFrame
#' @keywords classes internal
#' @seealso [read_hla_xml()], [parse_hla_alleles()], [HLARanges-class], [HLAGene][HLAGene_]
#' @export
#' @examples
#' showClass("HLAAllele")
setClass(
  Class = "HLAAllele",
  slots = list(
    locusname = "character",
    sequence  = "DNAStringSet",
    features  = "HLARangesList",
    metadata  = "DataFrame"
  )
)

#' Constructor for [HLAAllele-class] objects
#'
#' @note Don't run directly. This function is called by [parse_hla_alleles()].
#' @param nodes Allele nodes from a `hla.xml` object.
#' @param locusname A valid HLA locus name.
#' @param ncores The number of compute cores to use.
#'
#' @return A [HLAAllele-class] object
#' @seealso [parse_hla_alleles()]
#' @export
#' @examples
#' showClass("HLAAllele")
HLAAllele <- function(nodes, locusname, ncores = parallel::detectCores() - 2) {
  if (missing(nodes)) {
    return(new("HLAAllele"))
  }
  new("HLAAllele", nodes, locusname, ncores)
}

setMethod("initialize", signature(.Object = "HLAAllele"), function(.Object, nodes, locusname, ncores) {
  if (missing(nodes)) {
    .Object@locusname = NA_character_
    .Object@sequence  = Biostrings::DNAStringSet()
    .Object@features  = HLARangesList()
    .Object@metadata  = S4Vectors::DataFrame()
  } else {
    lcn               <- match_hla_locus(locusname)
    .Object@locusname = lcn
    .Object@sequence  = HLAAllele_parser$parse_sequence(nodes)
    .Object@features  = HLAAllele_parser$parse_features(nodes, ncores)
    .Object@metadata  = HLAAllele_parser$parse_metadata(nodes, lcn)
  }
  .Object
})

setMethod("locusname", "HLAAllele", function(x, ...) x@locusname)

setMethod("locusname<-", "HLAAllele", function(x, ..., value) {
  value <- match_hla_locus(value)
  x@locusname <- value
  x
})

setMethod("sequences", "HLAAllele", function(x, ...) x@sequence)

setMethod("sequences<-", "HLAAllele", function(x, ..., value) {
  x@sequence <- value
  x
})

setMethod("features", "HLAAllele", function(x, ...) x@features)

setMethod("features<-", "HLAAllele", function(x, ..., value) {
  x@features <- value
  x
})

setMethod("elementMetadata", "HLAAllele", function(x, use.names = FALSE, ...) x@metadata)

setMethod("elementMetadata<-", "HLAAllele", function(x, ..., value) {
  x@metadata <- value
  x
})

#' @describeIn HLAAllele Get sequence names from [HLAAllele-class]s.
setMethod("names", "HLAAllele", function(x) names(sequences(x)))

#' @describeIn HLAAllele Get number of [HLAAllele-class]s.
setMethod("length", "HLAAllele", function(x) length(sequences(x)))

setMethod("exon", "HLAAllele", function(x, exon = NULL, ...) {
  ftr <- features(x)
  if (!is.null(exon)) {
    ## convert exon number to feature order
    fod <- exon * 2
    rng <- ranges(ftr)[BiocGenerics::`%in%`(getOrder(ftr), fod)]
  } else {
    rng <- ranges(ftr)[getType(ftr) == "Exon"]
  }
  exon_seq <- sequences(x)[rng]
  ftr_name <- vapply(rng, colon %.% names, FUN.VALUE = "", USE.NAMES = FALSE)
  names(exon_seq) <- paste0(names(exon_seq), "|", ftr_name)
  exon_seq
})

setMethod("intron", "HLAAllele", function(x, intron = NULL, ...) {
  ftr <- features(x)
  if (!is.null(intron)) {
    ## convert intron number to feature order
    fod <- intron * 2 + 1
    rng <- ranges(ftr)[BiocGenerics::`%in%`(getOrder(ftr), fod)]
  } else {
    rng <- ranges(ftr)[getType(ftr) == "Intron"]
  }
  ## exclude UTRs
  rng <- rng[IRanges::LogicalList(lapply(rng, function(x) !grepl("UTR", names(x))))]
  intron_seq <- sequences(x)[rng]
  ftr_name <- vapply(rng, colon %.% names, FUN.VALUE = "", USE.NAMES = FALSE)
  names(intron_seq) <- paste0(names(intron_seq), "|", ftr_name)
  intron_seq
})

setMethod("utr", "HLAAllele", function(x, utr = NULL, ...) {
  ftr <- features(x)
  if (!is.null(utr)) {
    stopifnot(all(utr %in% 1:2))
    fod <- unique(unlist(IRanges::which(getType(ftr) == "UTR")))
    if (all(fod == 1) && all(utr == 2)) {
      stop("No 3' UTR available")
    }
    if (all(fod > 1) && all(utr == 1)) {
      stop("No 5' UTR available")
    }
    fod <- fod[utr]
    rng <- ranges(ftr)[BiocGenerics::`%in%`(getOrder(ftr), fod)]
  } else {
    rng <- ranges(ftr)[getType(ftr) == "UTR"]
  }
  utr_seq <- sequences(x)[rng]
  ftr_name <- vapply(rng, colon %.% names, FUN.VALUE = "", USE.NAMES = FALSE)
  names(utr_seq) <- paste0(names(utr_seq), "|", ftr_name)
  utr_seq
})

setMethod("noutr", "HLAAllele", function(x, ...) {
  rng0 <- IRanges::lapply(ranges(features(x, ...)), function(r0) {
    if (all(c("5' UTR", "3' UTR") %in% names(r0))) {
      IRanges::IRanges(start = end(r0["5' UTR"]) + 1L, end = start(r0["3' UTR"]) - 1L)
    } else if ("5' UTR" %in% names(r0)) {
      IRanges::IRanges(start = end(r0["5' UTR"]) + 1L, end = max(end(r0)))
    } else if ("3' UTR" %in% names(r0)) {
      IRanges::IRanges(start = min(start(r0)), end = start(r0["3' UTR"]) - 1L)
    } else {
      IRanges::IRanges(start = min(start(r0)), end =  max(end(r0)))
    }
  })
  noutr_seq <- sequences(x)[rng0]
  noutr_seq
})

#' @describeIn HLAAllele Combine [HLAAllele-class] objects.
setMethod("c", signature(x = "HLAAllele"), function(x, ..., recursive = FALSE) {
  args <- unname(list(x, ...))
  ans <- HLAAllele()
  if (length(args) == 2) {
    x <- args[[1]]
    y <- args[[2]]
    stopifnot(locusname(x) == locusname(y))
    locusname(ans) <- locusname(x)
    sequences(ans) <- c(sequences(x), sequences(y))
    features(ans)  <- suppressWarnings(c(features(x), features(y)))
    elementMetadata(ans) <- rbind(elementMetadata(x), elementMetadata(y))
  } else {
    stopifnot(unique(sapply(args, locusname)) == 1)
    locusname(ans) <- unique(sapply(args, locusname))
    sequences(ans) <- do.call("c", lapply(args, sequences))
    features(ans)  <- suppressWarnings(do.call("c", lapply(args, features)))
    elementMetadata(ans) <- do.call("rbind", lapply(args, elementMetadata))
  }
  ans
})

#' @describeIn HLAAllele Subset [HLAAllele-class] objects.
setMethod("[", signature(x = "HLAAllele", i = "numeric", j = "missing"), function(x, i, j, ..., drop = TRUE) {
  ans <- HLAAllele()
  locusname(ans) <- locusname(x)
  sequences(ans) <- sequences(x)[i]
  features(ans)  <- features(x)[i]
  elementMetadata(ans) <- elementMetadata(x)[i, ]
  ans
})

#' @describeIn HLAAllele Subset [HLAAllele-class] objects.
setMethod("[", signature(x = "HLAAllele", i = "logical", j = "missing"), function(x, i, j, ..., drop = TRUE) {
  ans <- HLAAllele()
  i <- which(i)
  if (length(i) == 0) {
    return(ans)
  }
  locusname(ans) <- locusname(x)
  sequences(ans) <- sequences(x)[i]
  features(ans)  <- features(x)[i]
  elementMetadata(ans) <- elementMetadata(x)[i, ]
  ans
})

#' @describeIn HLAAllele Subset [HLAAllele-class] objects.
setMethod("[", signature(x = "HLAAllele", i = "character", j = "missing"), function(x, i, j, ..., drop = TRUE) {
  ans <- HLAAllele()
  i   <- match_alleles(i, x, partially = TRUE)

  ## if there is no match return empty object
  if (length(i) == 0) {
    return(ans)
  }

  locusname(ans) <- locusname(x)
  sequences(ans) <- sequences(x)[i]
  features(ans)  <- features(x)[i]
  elementMetadata(ans) <- elementMetadata(x)[i, ]
  ans
})

setMethod("allele_id", signature(x = "HLAAllele"), function(x, ...) {
  elementMetadata(x)$allele_id
})

setMethod("allele_name", signature(x = "HLAAllele"), function(x, ...) {
  elementMetadata(x)$allele_name
})

setMethod("g_group", signature(x = "HLAAllele"), function(x, ...) {
  elementMetadata(x)$g_group
})

setMethod("p_group", signature(x = "HLAAllele"), function(x, ...) {
  elementMetadata(x)$p_group
})

setMethod("cwd_status", signature(x = "HLAAllele"), function(x, ...) {
  elementMetadata(x)$cwd_status
})

setMethod("ancestry", signature(x = "HLAAllele"), function(x, ...) {
  elementMetadata(x)$ancestry
})

setMethod("sample_name", signature(x = "HLAAllele"), function(x, ...) {
  elementMetadata(x)$sample
})

setMethod("is_complete", signature(x = "HLAAllele"), function(x, ...) {
  elementMetadata(x)$complete
})

setMethod("is_lsl", signature(x = "HLAAllele"), function(x, ...) {
  pttrn <- ".*DKMS-LSL.*$"
  grepl(pttrn, elementMetadata(x)$sample)
})

setAs(from = "HLAAllele", to = "data.table", function(from) {
  alleles <- strsplit(allele_name(from), "*", fixed = TRUE)
  fts <- features(from)
  dtplyr::lazy_dt(data.table(
    gene = vapply(alleles, `[[`, 1, FUN.VALUE = ""),
    allele_name = vapply(alleles, `[[`, 2, FUN.VALUE = ""),
    data_assigned = lubridate::ymd(hlatools::elementMetadata(from)$date_assigned),
    cwd_status = cwd_status(from),
    ancestry = gsub(":", "|", ancestry(from)),
    exons = vapply(fts, function(x) collapse(getOrder(x)[getType(x) == "Exon"] / 2), FUN.VALUE = "", USE.NAMES = FALSE),
    exon_status = vapply(fts, function(x) collapse(substr(getStatus(x), 1, 1)[getType(x) == "Exon"]), FUN.VALUE = "", USE.NAMES = FALSE)
  ))
})

#' @export
setMethod("as.data.table", signature(x = "HLAAllele"), function(x, ...) {
  as(x, "data.table")
})

show_HLAAllele <- function(x) {
  if (length(x) == 0) {
    cat("An empty object of class ", sQuote(class(x)), " for ", sQuote(locusname(x)), "\n", sep = "")
  } else {
    SeqLen <- width(sequences(x))
    FeatureType <- unname(vapply(getType(features(x)), paste0, collapse = ":", FUN.VALUE = "", USE.NAMES = FALSE))
    cat("An object of class ", sQuote(class(x)), " for ", sQuote(locusname(x)), "\n", sep = "")
    show(S4Vectors::DataFrame(elementMetadata(x)[, c("allele_name", "g_group", "p_group", "cwd_status")], SeqLen, FeatureType))
  }
}

setMethod("show", "HLAAllele", function(object) {
  show_HLAAllele(object)
})
