#' @importFrom R6 R6Class
NULL


# Class: NGSSample --------------------------------------------------------


#' Constructor for \code{\link[=NGSSample_]{NGSSample}} objects
#'
#' @param xmlpath Path to an NGSengine NGSSample XML file.
#' @param lims_donor_id [Optional] Sample ID.
#'
#' @return A \code{\link[=NGSSample_]{NGSSample}} object.
#' @export
#' @keywords internal
#' @examples
#' \dontrun{
#' NGSSample("path-to-NGSengine-xml")
#' }
NGSSample <- function(xmlpath, lims_donor_id = NULL) {
  NGSSample_$new(xmlpath, lims_donor_id)
}

#' Class \code{"NGSSample"}
#'
#' @docType class
#' @usage NGSSample(xmlpath)
#' @field path \code{[character]}.
#' @field metadata \code{[list]}.
#' @field region \code{[Biostrings::BStringSet]}.
#' @field genotypes \code{[list]}.
#'
#' @keywords data internal
#' @return Object of \code{\link{R6Class}} representing an NGSengine NGSSample XML.
#' @section Methods:
#' \describe{
#'   \item{\code{x$new(xmlpath)}}{Create object of this class.}
#'   \item{\code{x$allele1()}}{First allele.}
#'   \item{\code{x$allele2()}}{Second allele.}
#'   \item{\code{x$locusname()}}{}
#'   \item{\code{x$samplename()}}{}
#'   \item{\code{x$lims_donor_id()}}{}
#'   \item{\code{x$xml_path()}}{}
#' }
NGSSample_ <- R6::R6Class(
  classname = "NGSSample",
  public = list(
    path = NULL,      # [character]
    metadata = NULL,  # [list]
    region = NULL,    # [Biostrings::BStringSet]
    genotypes = NULL, # [list]
    initialize = function(xmlpath, lims_donor_id = NULL) {
      stopifnot(requireNamespace("XML", quietly = TRUE))
      doc <- XML::xmlInternalTreeParse(xmlpath)
      self$path <- dirname(xmlpath)

      source_  <- strsplitN(xval(doc, "/NGSSample/Source"), "[\\/]", 1, from = "end")
      barcode_ <- if (grepl("lbc\\d+", source_)) {
        m <- regexpr("lbc(\\d+)", source_, perl = TRUE)
        cs <- attr(m, "capture.start")
        cl <- attr(m, "capture.length")
        as.integer(substr(source_, cs, cs + cl - 1))
      } else NA_integer_
      self$metadata <- list(
        version = xattr(doc, "/NGSSample", "Version"),
        locus   = paste0("HLA-", xval(doc, "/NGSSample/Locus")),
        library = xval(doc, "/NGSSample/Library"),
        sample_name = xval(doc, "/NGSSample/SampleName"),
        source  = source_,
        barcode = barcode_,
        lims_donor_id = lims_donor_id
      )

      self$region <- Biostrings::BStringSet(c(
        seq0  = xval(doc, "/NGSSample/Phasing/region/sequence[1]"),
        seq1  = xval(doc, "/NGSSample/Phasing/region/sequence[2]"),
        p     = xval(doc, "/NGSSample/Phasing/region/p"),
        flags = xval(doc, "/NGSSample/Phasing/region/flags")
      ), use.names = TRUE)

      self$genotypes <- list(best_allele(doc, 1), best_allele(doc, 2))
    }
  )
)


# Methods: NGSSample ------------------------------------------------------


NGSSample_$set("public", "pacbio_reads", function() {
  if (dir.exists(file.path(self$path, "pacbio"))) {
    path <- file.path(self$path, "pacbio", self$metadata$source)
    tryCatch({
      normalizePath(path, mustWork = TRUE)
    }, error = function(e) {
      NULL
    })
  } else NULL
})

NGSSample_$set("public", "nanopore_reads", function() {
  if (dir.exists(file.path(self$path, "nanopore"))) {
    path <- file.path(self$path, "nanopore", self$metadata$source)
    tryCatch({
      normalizePath(path, mustWork = TRUE)
    }, error = function(e) {
      NULL
    })
  } else NULL
})

NGSSample_$set("public", "illumina_reads", function() {
  if (!is.null(ldid <- self$metadata$lims_donor_id)) {
    pattern <- paste0(self$metadata$lims_donor_id, "_.+_R[12]_\\d+.fastq")
    if (dir.exists(file.path(self$path, "illumina"))) {
      path <- dir(file.path(self$path, "illumina"), pattern, full.names = TRUE)
    }
    tryCatch({
      normalizePath(path, mustWork = TRUE)
    }, error = function(e) {
      NULL
    })
  } else NULL
})

NGSSample_$set("public", "allele1", function() {
  self$genotypes[[1]]$allele
})

NGSSample_$set("public", "allele2", function() {
  self$genotypes[[2]]$allele
})

NGSSample_$set("public", "locusname", function() {
  self$metadata$locus
})

NGSSample_$set("public", "samplename", function() {
  self$metadata$sample_name
})

NGSSample_$set("public", "lims_donor_id", function() {
  self$metadata$lims_donor_id
})

NGSSample_$set("public", "xml_path", function() {
  snm <- self$samplename()
  lnm <- sub("^HLA-", "", self$locusname())
  normalizePath(
    file.path(self$path, paste0(snm, ".", lnm, ".xml")),
    mustWork = TRUE
  )
})

# Helpers -----------------------------------------------------------------


best_allele <- function(doc, num = 1) {
  xp <- paste0("/NGSSample/Typing/Genotypes/Genotype/Alleles[", num, "]")
  ph <- as.integer(strsplit(xattr(doc, xp, "PhasingConfiguration"), " ")[[1]])
  xp2 <- paste0(xp, "/Allele[1]")
  allele <- xval(doc, xp2)
  anm <- paste0("HLA-", strsplit(allele, " ")[[1]][1])
  mm <- strsplit(strsplit(allele, " ")[[1]][-1], "=")
  mm <- setNames(as.integer(vapply(mm, `[`, 2, FUN.VALUE = "")), vapply(mm, `[`, 1, FUN.VALUE = ""))
  list(allele = anm, phasing = ph, mismatches = mm)
}
