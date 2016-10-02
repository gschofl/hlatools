## declare "node" global so that "codetools" don't complain
## see "http://r.789695.n4.nabble.com/R-CMD-check-and-foreach-code-td4687660.html"
utils::globalVariables("node", package = "hlatools")

#' Parse all HLA alleles for a locus from hla.xml
#'
#' @param doc \file{hla.xml} as an \code{xml_document} object.
#' @param locusname One of [\code{HLA-A, HLA-B, HLA-C, HLA-DPA1, HLA-DPB1,
#' HLA-DQA1, HLA-DQB1, HLA-DRB}]
#' @param ncores he number of compute cores to use.
#'
#' @return A \code{\linkS4class{HLAAllele}} object.
#' @export
#' @examples
#' \dontrun{
#' doc <- read_hla_xml()
#' dpb1 <- parse_hla_alleles(doc, "HLA-DPB1")
#' }
parse_hla_alleles <- function(doc, locusname, ncores = parallel::detectCores()) {
  drb_name  <- NULL
  locusname <- match_hla_locus(locusname)
  if (startsWith(locusname, "HLA-DRB")) {
    drb_name  <- locusname
    locusname <- "HLA-DRB"
  }
  ns <- xml2::xml_ns(doc)
  xpath <- paste0("/d1:alleles/d1:allele/d1:locus[@locusname='", locusname, "']/parent::node()")
  rs <- HLAAllele(nodes = xml2::xml_find_all(doc, xpath, ns), ncores = ncores)
  ## nasty hack to work around the fact that all DRBs are put together
  ## in hla.xml
  if (!is.null(drb_name)) {
    rs <- rs[which(startsWith(allele_name(rs), drb_name))]
  }
  rs
}
