#' @importFrom R6 R6Class
#' @importFrom Biostrings BStringSet
NULL

#' @export
NGSSample <- R6::R6Class(
  classname = "NGSSample",
  public = list(
    path = NULL,      # [character]
    metadata = NULL,  # [list]
    region = NULL,    # [Biostrings::BStringSet]
    genotypes = NULL, # [list]
    initialize = function(path) {
      stopifnot(requireNamespace("XML", quietly = TRUE))
      doc <- XML::xmlInternalTreeParse(path)
      self$path <- dirname(path)

      sample <- strsplitN(xval(doc, "/NGSSample/Source"), "\\", 1, from = "end", fixed = TRUE)
      bc <- if (grepl("lbc\\d+", sample)) {
        m <- regexpr("lbc(\\d+)", sample, perl = TRUE)
        cs <- attr(m, "capture.start")
        cl <- attr(m, "capture.length")
        as.integer(substr(sample, cs, cs + cl - 1))
      } else NA_integer_
      self$metadata <- list(
        locus = paste0("HLA-", xval(doc, "/NGSSample/Locus")),
        library = xval(doc, "/NGSSample/Library"),
        sample  = sample,
        bc = bc
      )

      self$region <- BStringSet(c(
        seq0  = xval(doc, "/NGSSample/Phasing/region/sequence[1]"),
        seq1  = xval(doc, "/NGSSample/Phasing/region/sequence[2]"),
        p     = xval(doc, "/NGSSample/Phasing/region/p"),
        flags = xval(doc, "/NGSSample/Phasing/region/flags")
      ), use.names = TRUE)

      self$genotypes <- list(best_allele(doc, 1), best_allele(doc, 2))
    }
  )
)

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

NGSSample$set("public", "pacbio_path", function() {
  path <- file.path(self$path, "pacbio", self$metadata$sample)
  normalizePath(path, mustWork = TRUE)
})




