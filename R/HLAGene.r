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
#' @field cons \code{[DNAStringSet]}; Strict consensus of completely sequenced alleles.
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
    cons      = NULL, # [DNAStringSet]
    initialize = function(locusname, ncores = parallel::detectCores(), verbose = TRUE) {
      self$locusname <- match_hla_locus(locusname)
      self$alleles   <- parse_hla_alleles(read_hla_xml(), self$locusname, ncores)
      self$dm        <- calc_exon2_distance(self$alleles, verbose)
      self$cons      <- calc_consensus_string(self$alleles, self$locusname, verbose)
    },
    print = function() {
      fmt0 <- "HLA locus <%s>\n"
      cat(sprintf(fmt0, self$locusname))
      print(self$alleles)
      invisible(self)
    },
    closest_complete_neighbor_ = function(allele) {
      allele <- expand_hla_allele(allele, self$locusname)
      if (!allele %in% allele_name(self)) {
        stop("Allele ", dQuote(allele), " not found.", call. = FALSE)
      }
      nm_complete <- allele_name(self$alleles[is_complete(self$alleles)])
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
  allele <- expand_hla_allele(allele, self$locusname)
  if (!allele %in% allele_name(self)) {
    stop("Allele ", dQuote(allele), " not found.", call. = FALSE)
  }
  nm_complete <- allele_name(self$alleles[is_complete(self$alleles)])
  if (allele %in% nm_complete) {
    return(allele)
  }
  dmi <- self$dm[allele, nm_complete]
  names(which.min(dmi))
})

HLAGene_$set("public", "get_allele", function(allele) {
  allele <- expand_hla_allele(allele, self$locusname)
  self$alleles[allele]
})

HLAGene_$set("public", "get_reference_sequence", function(allele) {
  ref <- self$closest_complete_neighbor_(allele)
  if (allele == ref) {
    ref <- self$get_allele(allele)
    sref <- Biostrings::BStringSet(sequences(ref)[[1]])
    names(sref) <- as(features(ref)[[1]], "character")
  } else {
    ref <- self$get_allele(ref)
    alt <- self$get_allele(allele)
    fref <- features(ref)[[1]]
    falt <- features(alt)[[1]]
    sref <- Biostrings::BString(tolower(sequences(ref)[[1]]))
    salt <- Biostrings::BString(toupper(sequences(alt)[[1]]))
    i <- match(names(falt), names(fref))
    Rref <- iter(ranges(fref)[i, ])
    Ralt <- iter(ranges(falt))
    while (hlatools::hasNext(Rref) && hlatools::hasNext(Ralt)) {
      rref <- nextElem(Rref)
      ralt <- nextElem(Ralt)
      subseq(sref, start(rref), end(rref)) <- subseq(salt, start(ralt), end(ralt))
    }
    sref <- Biostrings::BStringSet(sref)
    names(sref) <- as(merge(x = fref, y = falt), "character")
  }
  sref
})


HLAGene_$set("public", "get_all_reference_sequences", function(allele) {
  allele <- expand_hla_allele(allele, self$locusname)
  complete_alleles <- allele_name(self$alleles[is_complete(self$alleles)])
  if (allele %in% complete_alleles) {
    ref <- self$get_allele(allele)
    sref <- Biostrings::BStringSet(sequences(ref)[[1]])
    names(sref) <- as(features(ref)[[1]], "character")
    sref
  } else {
    alt  <- self$get_allele(allele)
    refs <- self$alleles[complete_alleles]
    falt <- features(alt)[[1]]
    frefs <- features(refs)
    salt <- Biostrings::BString(toupper(sequences(alt)[[1]]))
    srefs <- Biostrings::BStringSet(tolower(sequences(refs)))
    foreach(fref = frefs, sref = srefs, .combine = "c") %do% {
      i <- match(names(falt), names(fref))
      Rref <- iter(ranges(fref)[i, ])
      Ralt <- iter(ranges(falt))
      while (hlatools::hasNext(Rref) && hlatools::hasNext(Ralt)) {
        rref <- nextElem(Rref)
        ralt <- nextElem(Ralt)
        subseq(sref, start(rref), end(rref)) <- subseq(salt, start(ralt), end(ralt))
      }
      sref <- Biostrings::BStringSet(sref)
      names(sref) <- as(merge(x = fref, y = falt), "character")
      sref
    }
  }
})

HLAGene_$set("public", "get_reference_sequence_as_ft", function(allele) {
  allele <- expand_hla_allele(allele, self$locusname)
  x <- self$get_reference_sequence(allele)
  sprintf(
    "%s\n%s\n%s%s\n",
    hla_dat_id(x),
    hla_dat_de(x),
    hla_dat_ft(x),
    hla_dat_sq(x)
  )
})

HLAGene_$set("public", "get_all_reference_sequences_as_ft", function(allele) {
  allele <- expand_hla_allele(allele, self$locusname)
  x <- self$get_all_reference_sequences(allele)
  n <- length(x)
  rs <- character(n)
  for (i in seq_len(n)) {
    rs[i] <- sprintf(
      "%s\n%s\n%s%s\n",
      hlatools:::hla_dat_id(x[i]),
      hlatools:::hla_dat_de(x[i]),
      hlatools:::hla_dat_ft(x[i]),
      hlatools:::hla_dat_sq(x[i])
    )
  }
  paste0(rs, collapse = "")
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
  aln <- DECIPHER::AlignSeqs(exon2, iterations = 0, refinements = 0,
                             restrict = -500, verbose = verbose)
  DECIPHER::DistanceMatrix(aln, includeTerminalGaps = TRUE, verbose = verbose)
}

calc_consensus_string <- function(x, locusname = "", verbose = TRUE) {
  stopifnot(
    requireNamespace("DECIPHER", quietly = TRUE),
    requireNamespace("Biostrings", quietly = TRUE)
  )
  NUC <- c("A", "C", "G", "T")
  seqs <- sequences(x[is_complete(x)])
  aln <- DECIPHER::AlignSeqs(seqs, verbose = TRUE)
  cm <- t(Biostrings::consensusMatrix(aln))[, NUC]
  cm <- sweep(cm, 1, .rowSums(cm, NROW(cm), NCOL(cm)), `/`)
  cseq <- Biostrings::DNAStringSet(paste0(NUC[apply(cm, 1, which.max)], collapse = ""))
  names(cseq) <- trimws(paste0(locusname, " consensus"))
  cseq
}

hla_dat_id <- function(x) {
  paste0("ID   HLAxxxxx; SV 1; standard; DNA; HUM; ", width(x), " BP.")
}

hla_dat_de <- function(x) {
  s <- names(x)
  s <- strsplit(s, "~", fixed = TRUE)[[1]][1]
  s <- strsplitN(strsplit(s, "|", fixed = TRUE)[[1]], "[.]", 1)
  s <- if (length(s) > 1) {
    paste0("DE   ", s[1], "[", strsplitN(s[2], "*", 2, fixed = TRUE), "]")
  } else {
    paste0("DE   ", s)
  }
  paste0(s, ", Human MHC sequence")
}

hla_dat_ft <- function(x) {
  tbl <- as(names(x), "HLARanges")
  nm <- names(tbl)
  nm <- ifelse(grepl("UTR$", nm), "UTR", nm)
  status   <- getStatus(tbl)
  keys     <- strsplitN(nm, " ", 1)
  numbers  <- strsplitN(nm, " ", 2)
  keys     <- ifelse(keys == "UTR", keys, tolower(keys))
  numbers  <- ifelse(numbers == "UTR", "", numbers)
  #partials <- ifelse(is.na(status) | status == "Complete", FALSE, TRUE)
  foreach(k = keys, s = start(tbl), e = end(tbl), n = numbers, .combine = "paste0") %do% {
    l1 <- sprintf("FT   %-16s%s..%s\n", k, s, e)
    l2 <- if (nzchar(n)) {
      sprintf('FT                   /number="%s"\n', n)
    } else ""
    paste0(l1, l2)
  }
}

hla_dat_sq <- function(x) {
  l <- Biostrings::width(x)
  l1 <- paste0("SQ   Sequence ", l, " BP;\n")
  i <- seq(1, l, by = 60)
  if (i[length(i)] == l) {
    i <- i[-length(i)]
  }
  j <- unique(c(seq(61, l, by = 60), l))
  l2 <- paste0(vapply(Biostrings::substring(x, i, j), function(x) paste0("     ", paste0(
    substring(x, c(1, 11, 21, 31, 41, 51), c(10, 20, 30, 40, 50, 60)),
    collapse = " "
  )), FUN.VALUE = ""), collapse = "\n")
  l3 <- "\n//"
  paste0(l1, l2, l3)
}



