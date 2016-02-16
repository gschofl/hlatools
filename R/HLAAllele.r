#' @include HLARanges.r
NULL

#' Class \code{"HLAAllele"}
#'
#' A container for data parsed from the IMGT/HLA hla.xml file.
#' Part of a \code{\link{HLAGene}} object.
#'
#' @slot sequence A \code{\linkS4class{DNAStringSet}} object.
#' @slot metadata A \code{\linkS4class{DataFrame}} object.
#' @slot features A \code{\linkS4class{HLARangesList}} object.
#'
#' @importClassesFrom Biostrings DNAStringSet
#' @importClassesFrom S4Vectors DataFrame
#' @keywords classes internal
#' @seealso \code{\link{HLAGene}}.
#' @export
#' @examples
#' showClass("HLAAllele")
setClass(
  Class = "HLAAllele",
  slots = list(
    sequence = "DNAStringSet",
    features = "HLARangesList",
    metadata = "DataFrame"
  )
)

#' Constructor for \code{\linkS4class{HLAAllele}} objects
#'
#' @note Don't run directly. This function is called by \code{\link{parse_hla_alleles}}.
#' @param node An allele node from a hla.xml object.
#'
#' @return A \code{\linkS4class{HLAAllele}} object
#' @seealso \code{\link{parse_hla_alleles}}
#' @export
#' @examples
#' showClass("HLAAllele")
HLAAllele <- function(node) {
  if (missing(node)) {
    return(new("HLAAllele"))
  }
  new("HLAAllele", node)
}

setMethod("initialize", signature(.Object = "HLAAllele"), function(.Object, node) {
  if (missing(node)) {
    .Object@sequence = Biostrings::DNAStringSet()
    .Object@features = HLARangesList()
    .Object@metadata = DataFrame()
  } else {
    .Object@sequence = HLAAllele_parser$parse_sequence(node)
    .Object@features = HLAAllele_parser$parse_features(node)
    .Object@metadata = HLAAllele_parser$parse_metadata(node)
  }
  .Object
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

#' @export
setMethod("names", "HLAAllele", function(x) names(sequences(x)))

#' @export
setMethod("length", "HLAAllele", function(x) length(sequences(x)))

#' @describeIn HLAAllele Combine \code{HLAAllele} objects.
setMethod("c", signature(x = "HLAAllele"), function (x, ..., recursive = FALSE) {
  args <- unname(list(x, ...))
  ans <- HLAAllele()
  if (length(args) == 2) {
    x <- args[[1]]
    y <- args[[2]]
    sequences(ans) <- c(sequences(x), sequences(y))
    features(ans) <- suppressWarnings(c(features(x), features(y)))
    elementMetadata(ans) <- rbind(elementMetadata(x), elementMetadata(y))
  } else {
    sequences(ans) <- do.call("c", lapply(args, sequences))
    features(ans) <- suppressWarnings(do.call("c", lapply(args, features)))
    elementMetadata(ans) <- do.call("rbind", lapply(args, elementMetadata))
  }
  ans
})

#' @describeIn HLAAllele Subset \code{HLAAllele} objects.
setMethod("[", signature(x = "HLAAllele", i = "numeric", j = "missing"), function (x, i, j, ..., drop = TRUE) {
  ans <- HLAAllele()
  sequences(ans) <- sequences(x)[i]
  features(ans) <- features(x)[i]
  elementMetadata(ans) <- elementMetadata(x)[i, ]
  ans
})

#' @describeIn HLAAllele Subset \code{HLAAllele} objects.
setMethod("[", signature(x = "HLAAllele", i = "logical", j = "missing"), function (x, i, j, ..., drop = TRUE) {
  ans <- HLAAllele()
  i <- which(i)
  if (length(i) == 0) {
    return(ans)
  }
  sequences(ans) <- sequences(x)[i]
  features(ans) <- features(x)[i]
  elementMetadata(ans) <- elementMetadata(x)[i, ]
  ans
})

#' @describeIn HLAAllele Subset \code{HLAAllele} objects.
setMethod("[", signature(x = "HLAAllele", i = "character", j = "missing"), function (x, i, j, ..., drop = TRUE) {
  ans <- HLAAllele()
  i <- which( names(x) == expand_allele(i) )
  if (length(i) == 0) {
    return(ans)
  }
  sequences(ans) <- sequences(x)[i]
  features(ans) <- features(x)[i]
  elementMetadata(ans) <- elementMetadata(x)[i, ]
  ans
})

setMethod("allele_id", signature(x = "HLAAllele"), function(x, ...) {
  elementMetadata(x)$allele_id
})

setMethod("allele_name", signature(x = "HLAAllele"), function(x, ...) {
  elementMetadata(x)$allele_name
})

setMethod("cwd_status", signature(x = "HLAAllele"), function(x, ...) {
  elementMetadata(x)$cwd_status
})

setMethod("ethnicity", signature(x = "HLAAllele"), function(x, ...) {
  elementMetadata(x)$ethnicity
})

setMethod("is_complete", signature(x = "HLAAllele"), function(x, ...) {
  elementMetadata(x)$complete
})

setAs(from = "HLAAllele", to = "data.table", function(from) {
  alleles <- strsplit(allele_name(from), "*", fixed = TRUE)
  fts <- features(from)
  dplyr::tbl_dt(data.table(
    gene = vapply(alleles, `[[`, 1, FUN.VALUE = ""),
    allele_name = vapply(alleles, `[[`, 2, FUN.VALUE = ""),
    data_assigned = lubridate::ymd(hlatools::elementMetadata(from)$date_assigned),
    cwd_status = cwd_status(from),
    ethnicity = gsub(":", "|", ethnicity(from)),
    exons = vapply(fts, function(x) collapse(getOrder(x)[getType(x) == "Exon"]/2), FUN.VALUE = "", USE.NAMES = FALSE),
    exon_status = vapply(fts, function(x) collapse(substr(getStatus(x), 1, 1)[getType(x) == "Exon"]), FUN.VALUE = "", USE.NAMES = FALSE)
  ))
})

#' @export
setMethod("as.data.table", signature(x = "HLAAllele"), function(x, ...) {
  as(x, "data.table")
})

show_HLAAllele <- function(x) {
  SeqLen <- width(sequences(x))
  FeatureType <- unname(vapply(getType(features(x)), paste0, collapse = ":", FUN.VALUE = "", USE.NAMES = FALSE))
  cat("An object of class ", sQuote(class(x)), "\n", sep = "")
  show(DataFrame(elementMetadata(x)[, 1:4], SeqLen, FeatureType))
}

setMethod("show", "HLAAllele", function(object) {
  show_HLAAllele(object)
})

make_hla_allele_parser <- function() {
  ns <- c("x" = "http://hla.alleles.org/xml")
  list(
    # Parse nuclear sequence from an hla.xml allele node
    #
    # @param node A hla.xml allele node.
    #
    # @return A DNAStringSet object.
    # @keywords internal
    parse_sequence = function(node) {
      xpexpr <- './x:sequence/x:nucsequence'
      nucseq <- Biostrings::DNAStringSet(xval(node, xpexpr, namespaces = ns))
      names(nucseq) <- XML::xmlGetAttr(node, "name")
      nucseq
    },
    # Parse metadata from gr_unlistan hla.xml allele node
    #
    # @param node A hla.xml allele node.
    #
    # @return A DataFrame object.
    # @keywords internal
    parse_metadata = function(node) {
      DataFrame(
        allele_name   = XML::xmlGetAttr(node, "name"),
        allele_id     = XML::xmlGetAttr(node, "id"),
        date_assigned = XML::xmlGetAttr(node, "dateassigned"),
        cwd_status = xattr(node, "./x:cwd_catalogue", "cwd_status", namespaces = ns),
        ## Has 5'UTR or 3'UTR feature
        complete   = xval(node, "count(./x:sequence/x:feature[@featuretype=\"UTR\"])>0", as = "logical", namespaces = ns),
        pmid       = colon(xattr(node, "./x:citations/x:citation", "pubmed", namespaces = ns)),
        ethnicity  = colon(xval(node, "./x:sourcematerial/x:ethnicity/x:sample_ethnicity", namespaces = ns))
      )
    },
    # Parse features from an hla.xml allele node
    #
    # @param node A hla.xml allele node.
    #
    # @return A GRangesList object.
    # @keywords internal
    parse_features = function(node) {
      xp <- './x:sequence/x:feature[not(@featuretype="Protein")]'
      fset <- xset(node, xp, namespaces = ns)
      allele_name = xmlGetAttr(node, "name")
      hr <- HLARanges(
        seqnames = allele_name,
        ranges =  IRanges(
          start = xattr(node, paste0(xp, "/x:SequenceCoordinates"), "start", as = "integer", namespaces = ns),
          end   = xattr(node, paste0(xp, "/x:SequenceCoordinates"), "end", as = "integer", namespaces = ns),
          names = xattr(node, xp, "name", namespaces = ns)
        ),
        id     = xattr(node, xp, "id", namespaces = ns),
        order  = xattr(node, xp, "order", as = "integer", namespaces = ns),
        type   = xattr(node, xp, "featuretype", namespaces = ns),
        status = vapply(fset, function(x) xmlGetAttr(x, "status") %||% NA_character_, FUN.VALUE = ""),
        frame  = vapply(fset, function(x) {
          xattr(x, "./x:cDNACoordinates", "readingframe", as = "integer", namespaces = ns)
        }, FUN.VALUE = 0L)
      )
      hrl <- HLARangesList(hr)
      names(hrl) <- allele_name
      hrl
    }
  )
}

HLAAllele_parser <- make_hla_allele_parser()




