#' @include HLAAllele.r
NULL


# Class: HLAGene ----------------------------------------------------------


#' Constructor for \code{\link[=HLAGene_]{HLAGene}} objects.
#'
#' @param locusname A valid HLA gene name.
#' @param ... Passed on.
#'
#' @return A \code{\link[=HLAGene_]{HLAGene}} object
#' @export
#' @examples
#' \dontrun{
#' x <- HLAGene("DPB1")
#' db_version(x)
#' locusname(x)
#' length(x)
#'
#' ## extract alleles with complete sequences
#' x[is_complete(x)]
#'
#' ## extract all Common alleles
#' x[cwd_status(x) == "Common"]
#' }
HLAGene <- function(locusname, ...) {
  HLAGene_$new(locusname, ...)
}

#' Class \code{"HLAGene"}
#'
#' @docType class
#' @usage HLAgene(locusname, ncores = parallel::detectCores(), with_dist = FALSE)
#' @keywords data internal
#' @importFrom XVector subseq subseq<-
#' @import foreach
#' @return Object of \code{\link{R6Class}} representing an HLA gene.
#' @section Methods:
#' \describe{
#'   \item{\code{x$new(locusname, ncores = parallel::detectCores(), with_dist = FALSE)}}{Create an object of this class.}
#'   \item{\code{x$get_db_version()}}{Get the IMGT/HLA database version.}
#'   \item{\code{x$get_locusname()}}{Get the name of the locus.}
#'   \item{\code{x$get_alleles(allele)}}{Get alleles.}
#'   \item{\code{x$get_closest_complete_neighbor(allele)}}{Get the complete allele that is closest at exon 2 to the query allele.}
#'   \item{\code{x$get_reference_sequence(allele)}}{Get the (imputed) reference sequence for allele.}
#' }
HLAGene_ <- R6::R6Class(
  classname = "HLAGene",
  public = list(
    initialize = function(locusname, ncores = parallel::detectCores(), with_dist = FALSE) {
      doc <- read_hla_xml()
      private$dbv <- xml2::xml_attr(xml2::xml_find_all(doc, "//d1:alleles/d1:allele[1]/d1:releaseversions"), "currentrelease")
      private$lcn <- match_hla_locus(locusname)
      private$all <- parse_hla_alleles(doc, private$lcn, ncores)
      if (with_dist) {
        private$dmt <- calc_exon2_distance(private$all, verbose = TRUE)
        private$cns <- calc_consensus_string(private$all, private$lcn, verbose = TRUE)
      }
    },
    print = function() {
      fmt0 <- "IMGT/HLA database <%s>; Locus <%s>\n"
      cat(sprintf(fmt0, self$get_db_version(), self$get_locusname()))
      print(self$get_alleles())
      invisible(self)
    },
    ## getters and setters
    get_db_version = function() {
    private$dbv
    },
    get_locusname = function() {
      private$lcn
    },
    get_alleles = function(allele) {
      if (missing(allele)) {
        return(private$all)
      }
      allele <- if (is.character(allele)) {
        expand_hla_allele(allele, self$get_locusname())
      } else allele
      private$all[allele]
    },
    has_distances = function() {
      is.matrix(private$dmt)
    }
  ),
  private = list(
    dbv = NULL, # [character]; IMGT/HLA database version
    lcn = NULL, # [character]; locus name
    all = NULL, # [HLAAllele]; alleles
    dmt = NULL, # [matrix]; distance matrix based on exon 2 (it's the only one that is always present)
    cns = NULL  # [DNAStringSet]; consensus sequence based on all full-length alleles
  )
)

# Methods: HLAGene --------------------------------------------------------

HLAGene_$set("public", "get_closest_complete_neighbor", function(allele) {
  allele <- expand_hla_allele(allele, self$get_locusname())
  if (!allele %in% allele_name(self)) {
    stop("Allele ", dQuote(allele), " not found.", call. = FALSE)
  }
  nm_complete <- allele_name(self$get_alleles(is_complete(self)))
  if (allele %in% nm_complete) {
    return(allele)
  }
  if (!self$has_distances()) {
    private$dmt <- calc_exon2_distance(self$get_alleles(), verbose = TRUE)
  }
  dmi <- private$dmt[allele, nm_complete]
  names(which.min(dmi))
})

HLAGene_$set("public", "get_reference_sequence", function(allele) {
  ref <- self$get_closest_complete_neighbor(allele)
  if (allele == ref) {
    ref <- self$get_alleles(allele)
    sref <- Biostrings::BStringSet(sequences(ref)[[1]])
    names(sref) <- as(features(ref)[[1]], "character")
  } else {
    ref <- self$get_alleles(ref)
    alt <- self$get_alleles(allele)
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
  allele <- expand_hla_allele(allele, self$get_locusname())
  complete_alleles <- allele_name(self$get_alleles(is_complete(self)))
  if (allele %in% complete_alleles) {
    ref <- self$get_alleles(allele)
    sref <- Biostrings::BStringSet(sequences(ref)[[1]])
    names(sref) <- as(features(ref)[[1]], "character")
    sref
  } else {
    alt  <- self$get_alleles(allele)
    refs <- self$get_alleles(complete_alleles)
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
  allele <- expand_hla_allele(allele, self$get_locusname())
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
  allele <- expand_hla_allele(allele, self$get_locusname())
  x <- self$get_all_reference_sequences(allele)
  n <- length(x)
  rs <- character(n)
  for (i in seq_len(n)) {
    rs[i] <- sprintf(
      "%s\n%s\n%s%s\n",
      hla_dat_id(x[i]),
      hla_dat_de(x[i]),
      hla_dat_ft(x[i]),
      hla_dat_sq(x[i])
    )
  }
  paste0(rs, collapse = "")
})

# Formal Methods (S4) -----------------------------------------------------

setOldClass("HLAGene")

setMethod("sequences", signature(x = "HLAGene"), function(x, ...) {
  sequences(x$get_alleles())
})

setMethod("features", signature(x = "HLAGene"), function(x, ...) {
  features(x$get_alleles())
})

setMethod("elementMetadata", signature(x = "HLAGene"), function(x, use.names = FALSE, ...) {
  elementMetadata(x$get_alleles())
})

setMethod("allele_id", signature(x = "HLAGene"), function(x, ...) {
  elementMetadata(x)$allele_id
})

setMethod("allele_name", signature(x = "HLAGene"), function(x, ...) {
  elementMetadata(x)$allele_name
})

setMethod("g_group", signature(x = "HLAGene"), function(x, ...) {
  elementMetadata(x)$g_group
})

setMethod("p_group", signature(x = "HLAGene"), function(x, ...) {
  elementMetadata(x)$p_group
})

setMethod("cwd_status", signature(x = "HLAGene"), function(x, ...) {
  elementMetadata(x)$cwd_status
})

setMethod("ethnicity", signature(x = "HLAGene"), function(x, ...) {
  elementMetadata(x)$ethnicity
})

setMethod("sample_name", signature(x = "HLAGene"), function(x, ...) {
  elementMetadata(x)$sample
})

setMethod("is_complete", signature(x = "HLAGene"), function(x, ...) {
  elementMetadata(x)$complete
})

setMethod("is_lsl", signature(x = "HLAGene"), function(x, ...) {
  pttrn <- ".*DKMS-LSL.*$"
  grepl(pttrn, elementMetadata(x)$sample)
})


# S3 methods ----------------------------------------------------------------------------------

#' @export
locusname.HLAGene <- function(x) {
  x$get_locusname()
}

#' @export
db_version.HLAGene <- function(x) {
  x$get_db_version()
}

#' @export
`[.HLAGene` <- function(x, i, j, ..., drop = TRUE) {
  x$get_alleles(i)
}

#' @export
names.HLAGene <- function(x) {
  names(sequences(x))
}

#' @export
length.HLAGene <- function(x) {
  length(sequences(x))
}

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



