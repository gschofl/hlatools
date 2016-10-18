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
  locusname <- match_hla_locus(locusname)
  ns <- xml2::xml_ns(doc)
  xpath1 <- paste0("/d1:alleles/d1:allele/d1:locus[@locusname='", locusname, "']/parent::node()")
  nodes1 <- xml2::xml_find_all(doc, xpath1, ns)
  xpath2 <- paste0(".//d1:releaseversions[not(starts-with(@releasestatus, 'Allele Deleted'))]/parent::node()")
  nodes2 <- xml2::xml_find_all(nodes1, xpath2, ns)
  rs <- HLAAllele(nodes = nodes2, ncores = ncores)
  rs
}
