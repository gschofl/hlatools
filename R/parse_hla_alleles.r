## declare "node" global so that "codetools" don't complain
## see "http://r.789695.n4.nabble.com/R-CMD-check-and-foreach-code-td4687660.html"
utils::globalVariables("node", package = "hlatools")

#' Parse all HLA alleles for a locus from hla.xml
#'
#' @param doc \file{hla.xml} as an \code{XMLInternalDocument} object.
#' @param locusname One of [\code{HLA-A, HLA-B, HLA-C, HLA-DPA1, HLA-DPB1,
#' HLA-DQA1, HLA-DQB1, HLA-DRB}]
#' @param ncores Number of cores
#'
#' @return A \code{\linkS4class{HLAAllele}} object.
#' @export
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @import foreach
#' @examples
#' \dontrun{
#' doc <- read_hla_xml()
#' dpb1 <- parse_hla_alleles(doc, "HLA-DPB1", 12)
#' }
parse_hla_alleles <- function (doc, locusname, ncores = detectCores()) {
  registerDoParallel(cl = ncores)
  drb_name  <- NULL
  locusname <- match_hla_locus(locusname)
  if (dkms::starts_with(locusname, "HLA-DRB")) {
    drb_name  <- locusname
    locusname <- "HLA-DRB"
  }
  ns <- c("x" = "http://hla.alleles.org/xml")
  xpexpr <- paste0("/x:alleles/x:allele/x:locus[@locusname='", locusname, "']/parent::node()")
  rs <- foreach(node = xset(doc, xpexpr, namespaces = ns), .combine = "c") %dopar% {
    tryCatch({
      HLAAllele(node)
    }, error = function(e) {
      warning(e$message, immediate. = TRUE)
      HLAAllele()
    })
  }
  ## nasty hack to work around the fact that all DRBs are put together
  ## in hla.xml
  if (!is.null(drb_name)) {
    rs <- rs[which(dkms::starts_with(allele_name(rs), drb_name))]
  }
  rs
}
