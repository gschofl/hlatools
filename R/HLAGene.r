#' @include HLAAllele.r
NULL


# Class: HLAGene ----------------------------------------------------------


#' Constructor for \code{\link[=HLAGene_]{HLAGene}} objects.
#'
#' @param locusname A valid HLA gene name.
#'
#' @return A \code{\link[=HLAGene_]{HLAGene}} object
#' @export
#' @examples
#' \dontrun{
#' HLAGene("DPB1")
#' }
HLAGene <- function(locusname) {
  HLAGene_$new(locusname)
}

#' Class \code{"HLAGene"}
#'
#' @docType class
#' @usage HLAgene(locusname)
#' @field locusname \code{[character]}; A valid HLA gene.
#' @field alleles \code{[HLAAllele]}; An \code{\linkS4class{HLAAllele}} object.
#' @field dm \code{[matrix]}; Genetic distances between exons 2.
#'
#' @keywords data internal
#' @importFrom XVector subseq subseq<-
#' @return Object of \code{\link{R6Class}} representing an HLA gene.
#' @section Methods:
#' \describe{
#'   \item{\code{x$new(locusname)}}{Create object of this class.}
#'   \item{\code{x$get_allele(allele)}}{Extract an allele.}
#'   \item{\code{x$get_closest_complete_neighbor(allele)}}{Get the complete allele that closest at exon 2 to the query allele}
#'   \item{\code{x$get_reference_sequence(allele)}}{}
#' }
HLAGene_ <- R6::R6Class(
  classname = "HLAGene",
  public = list(
    locusname = NULL, # [character]
    alleles   = NULL, # [HLAAllele]
    dm        = NULL, # [matrix]
    initialize = function(locusname, ncores = parallel::detectCores(), verbose = TRUE) {
      self$locusname <- match_hla_locus(locusname)
      self$alleles   <- parse_hla_alleles(read_hla_xml(), self$locusname, ncores)
      self$dm        <- calc_exon2_distance(self$alleles, verbose)
    }
  ),
  private = list(
    closest_complete_neighbor = function(allele) {
      allele <- expand_allele(allele, self$locusname)
      if (!allele %in% allele_name(self)) {
        stop("Allele ", dQuote(allele), " not found.", call. = FALSE)
      }
      nm_complete <- allele_name(self$alleles[which(is_complete(self$alleles))])
      if (allele %in% nm_complete) {
        return(allele)
      }
      dmi <- self$dm[allele, nm_complete]
      names(which.min(dmi))
    }
  )
)


# Methods: HLAGene --------------------------------------------------------


HLAGene_$set("public", "get_closest_complete_neighbor", function(allele) {
  allele <- expand_allele(allele, self$locusname)
  if (!allele %in% allele_name(self)) {
    stop("Allele ", dQuote(allele), " not found.", call. = FALSE)
  }
  nm_complete <- allele_name(self$alleles[which(is_complete(self$alleles))])
  if (allele %in% nm_complete) {
    return(allele)
  }
  dmi <- self$dm[allele, nm_complete]
  names(which.min(dmi))
})

HLAGene_$set("public", "get_allele", function(allele) {
  allele <- expand_allele(allele, self$locusname)
  self$alleles[allele]
})

HLAGene_$set("public", "get_reference_sequence", function(allele) {
  ref <- private$closest_complete_neighbor(allele)
  if (allele == ref) {
    ref <- self$get_allele(allele)
    sref <- Biostrings::BStringSet(sequences(ref)[[1]])
    names(sref) <- as(features(ref)[[1]], "character")
  } else {
    ref <- self$get_allele(ref)
    alt <- self$get_allele(allele)
    falt <- features(alt)[[1]]
    fref <- features(ref)[[1]]
    salt <- Biostrings::BString(sequences(alt)[[1]])
    sref <- Biostrings::BString(tolower(sequences(ref)[[1]]))
    i <- which(names(fref) %in% names(falt))
    Rref <- iter(ranges(fref)[i, ])
    Ralt <- iter(ranges(falt))
    while (hlatools::hasNext(Rref) && hlatools::hasNext(Ralt)) {
      rref <- nextElem(Rref)
      ralt <- nextElem(Ralt)
      subseq(sref, start(rref), end(rref)) <- subseq(salt, start(ralt), end(ralt))
    }
    sref <- Biostrings::BStringSet(sref)
    names(sref) <- as(merge(fref, falt), "character")
  }
  sref
})


# Formal Methods (S4) -----------------------------------------------------

setOldClass("HLAGene")

setMethod("allele_id", signature(x = "HLAGene"), function(x, ...) {
  elementMetadata(x$alleles)$allele_id
})

setMethod("allele_name", signature(x = "HLAGene"), function(x, ...) {
  elementMetadata(x$alleles)$allele_name
})

setMethod("cwd_status", signature(x = "HLAGene"), function(x, ...) {
  elementMetadata(x$alleles)$cwd_status
})

setMethod("is_complete", signature(x = "HLAGene"), function(x, ...) {
  elementMetadata(x$alleles)$complete
})

# Helpers -----------------------------------------------------------------


calc_exon2_distance <- function(x, verbose = TRUE) {
  stopifnot(requireNamespace("DECIPHER", quietly = TRUE))
  name <- "Exon 2"
  df <- as.data.frame(ranges(features(x)))
  df <- df[df$names == name, c("start", "end")]
  exon2 <- subseq(sequences(x), df$start, df$end)
  aln <- DECIPHER::AlignSeqs(exon2, verbose = verbose)
  DECIPHER::DistanceMatrix(aln, includeTerminalGaps = TRUE, verbose = verbose)
}
