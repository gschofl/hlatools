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
#' @param nodes Allele nodes from a hla.xml object.
#' @param ncores The number of compute cores to use.
#'
#' @return A \code{\linkS4class{HLAAllele}} object
#' @seealso \code{\link{parse_hla_alleles}}
#' @export
#' @examples
#' showClass("HLAAllele")
HLAAllele <- function(nodes, ncores = parallel::detectCores()) {
  if (missing(nodes)) {
    return(new("HLAAllele"))
  }
  new("HLAAllele", nodes, ncores)
}

setMethod("initialize", signature(.Object = "HLAAllele"), function(.Object, nodes, ncores) {
  if (missing(nodes)) {
    .Object@sequence = Biostrings::DNAStringSet()
    .Object@features = HLARangesList()
    .Object@metadata = S4Vectors::DataFrame()
  } else {
    .Object@sequence = HLAAllele_parser$parse_sequence(nodes)
    .Object@features = HLAAllele_parser$parse_features(nodes, ncores)
    .Object@metadata = HLAAllele_parser$parse_metadata(nodes)
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

#' @describeIn HLAAllele Get sequence names from \code{HLAAllele}s.
setMethod("names", "HLAAllele", function(x) names(sequences(x)))

#' @describeIn HLAAllele Get number of \code{HLAAllele}s.
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
setMethod("[", signature(x = "HLAAllele", i = "numeric", j = "missing"), function(x, i, j, ..., drop = TRUE) {
  ans <- HLAAllele()
  sequences(ans) <- sequences(x)[i]
  features(ans) <- features(x)[i]
  elementMetadata(ans) <- elementMetadata(x)[i, ]
  ans
})

#' @describeIn HLAAllele Subset \code{HLAAllele} objects.
setMethod("[", signature(x = "HLAAllele", i = "logical", j = "missing"), function(x, i, j, ..., drop = TRUE) {
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
setMethod("[", signature(x = "HLAAllele", i = "character", j = "missing"), function(x, i, j, ..., drop = TRUE) {
  ans <- HLAAllele()
  i <- match(expand_hla_allele(i), names(x))
  if (length(i <- i[!is.na(i)]) == 0) {
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

setMethod("g_group", signature(x = "HLAAllele"), function(x, ...) {
  elementMetadata(x)$g_group
})

setMethod("p_group", signature(x = "HLAAllele"), function(x, ...) {
  elementMetadata(x)$p_group
})

setMethod("cwd_status", signature(x = "HLAAllele"), function(x, ...) {
  elementMetadata(x)$cwd_status
})

setMethod("ethnicity", signature(x = "HLAAllele"), function(x, ...) {
  elementMetadata(x)$ethnicity
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
  dtplyr::tbl_dt(data.table(
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
  show(S4Vectors::DataFrame(elementMetadata(x)[, c("allele_name", "g_group", "p_group", "cwd_status")], SeqLen, FeatureType))
}

setMethod("show", "HLAAllele", function(object) {
  show_HLAAllele(object)
})

make_hla_allele_parser <- function() {
  list(
    # Parse nuclear sequence from an hla.xml allele node
    #
    # @param nodes A hla.xml allele nodeset.
    #
    # @return A DNAStringSet object.
    # @keywords internal
    parse_sequence = function(nodes) {
      ns     <- xml2::xml_ns(nodes)
      xpath  <- './d1:sequence/d1:nucsequence'
      nucseq <- Biostrings::DNAStringSet(xml2::xml_text(xml2::xml_find_all(nodes, xpath, ns)))
      names(nucseq) <- xml2::xml_attr(nodes, "name")
      nucseq
    },
    # Parse metadata from an hla.xml allele node
    #
    # @param nodes A hla.xml allele nodeset.
    #
    # @return A DataFrame object.
    # @keywords internal
    parse_metadata = function(nodes) {
      ns      <- xml2::xml_ns(nodes)
      cit_idx <- xml2::xml_find_lgl(nodes, "boolean(d1:citations)", ns)
      smp_idx <- xml2::xml_find_lgl(nodes, "boolean(d1:sourcematerial/d1:samples)", ns)
      S4Vectors::DataFrame(
        ##
        ## Allele designation
        ##
        allele_name   = xml2::xml_attr(nodes, "name"),
        allele_id     = xml2::xml_attr(nodes, "id"),
        g_group       = xml2::xml_find_chr(nodes, "string(./d1:hla_g_group/@status)", ns),
        p_group       = xml2::xml_find_chr(nodes, "string(./d1:hla_p_group/@status)", ns),
        date_assigned = xml2::xml_attr(nodes, "dateassigned"),
        ##
        ## Release
        ##
        first_released = xml2::xml_find_chr(nodes, "string(./d1:releaseversions/@firstreleased)", ns),
        last_updated   = xml2::xml_find_chr(nodes, "string(./d1:releaseversions/@lastupdated)", ns),
        release_status = xml2::xml_find_chr(nodes, "string(./d1:releaseversions/@releasestatus)", ns),
        confirmed      = xml2::xml_find_lgl(nodes, "string(./d1:releaseversions/@confirmed)=\"Confirmed\"", ns),
        ##
        ## CWD status and Completeness (we consider as complete alleles for which both UTRs are present)
        ##
        cwd_status    = xml2::xml_find_chr(nodes, "string(./d1:cwd_catalogue/@cwd_status)", ns),
        complete      = xml2::xml_find_lgl(nodes, "count(./d1:sequence/d1:feature[@featuretype=\"UTR\"])=2", ns),
        ##
        ## Source (PubMed ID, Ethnicity, Sample/Cellline)
        ##
        pmid          = ifelse(cit_idx, vapply(xml2::xml_find_all(nodes[cit_idx], "./d1:citations", ns), function(node) {
          colon(xml2::xml_attr(xml2::xml_children(node), "pubmed"))
        }, FUN.VALUE = character(1)), NA_character_),
        ethnicity     = vapply(xml2::xml_find_all(nodes, "./d1:sourcematerial/d1:ethnicity", ns), function(node) {
          colon(xml2::xml_text(xml2::xml_children(node)))
        }, FUN.VALUE = character(1)),
        sample       = ifelse(smp_idx, vapply(xml2::xml_find_all(nodes, "./d1:sourcematerial/d1:samples", ns), function(node) {
          colon(xml2::xml_attr(xml2::xml_children(node), "name"))
        }, FUN.VALUE = character(1)), NA_character_)
      )
    },
    # Parse features from an hla.xml allele node
    #
    # @param node A hla.xml allele node.
    #
    # @return A GRangesList object.
    # @keywords internal
    parse_features = function(nodes, ncores) {
      ns       <- xml2::xml_ns(nodes)
      nodeset  <- xml2::xml_find_all(nodes, "./d1:sequence", ns)
      xpath    <- "./d1:feature[not(@featuretype=\"Protein\")]"
      seqnames <- xml2::xml_attr(nodes, "name")
      rs <- HLARangesList(parallel::mcMap(function(seqname, node) {
        HLARanges(
          seqnames = seqname,
          ranges   = IRanges::IRanges(
            start = as.integer(xml2::xml_text(xml2::xml_find_all(node, paste0(xpath, "/d1:SequenceCoordinates/@start"), ns))),
            end   = as.integer(xml2::xml_text(xml2::xml_find_all(node, paste0(xpath, "/d1:SequenceCoordinates/@end"), ns))),
            names = xml2::xml_text(xml2::xml_find_all(node, paste0(xpath, "/@name"), ns))
          ),
          id     = xml2::xml_text(xml2::xml_find_all(node, paste0(xpath, "/@id"), ns)),
          order  = as.integer(xml2::xml_text(xml2::xml_find_all(node, paste0(xpath, "/@order"), ns))),
          type   = xml2::xml_text(xml2::xml_find_all(node, paste0(xpath, "/@featuretype"), ns)),
          status = vapply(xml2::xml_find_all(node, xpath, ns), xml2::xml_attr, "status", FUN.VALUE = ""),
          frame  = vapply(xml2::xml_find_all(node, xpath, ns), function(node) {
            as.integer(xml2::xml_find_chr(node, "string(./d1:cDNACoordinates/@readingframe)", ns))
          }, FUN.VALUE = 0L)
        )
      }, seqname = seqnames, node = nodeset, mc.cores = ncores))
      rs
    }
  )
}

HLAAllele_parser <- make_hla_allele_parser()

