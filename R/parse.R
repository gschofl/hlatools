## declare "node" global so that "codetools" don't complain
## see "http://r.789695.n4.nabble.com/R-CMD-check-and-foreach-code-td4687660.html"
utils::globalVariables("node", package = "hlatools")

#' Fetch or update the IPD-IMGT/HLA hla.xml file
#'
#' @param db_path <[character]>; location of local IPD-IMGT/HLA repository.
#' @param remote <[logical]>; if `TRUE` pull data from the IPD-IMGT/HLA ftp server,
#' if `FALSE` retrieve data from `db_path`.
#'
#' @return An [xml_document][xml2::read_xml()].
#' @export
#' @examples \dontrun{
#' doc <- read_hla_xml(remote = TRUE)
#' update_hla_xml()
#' }
read_hla_xml <- function(db_path = getOption("hlatools.local_repos"), remote = FALSE) {
  tdir <- tempdir()
  if (remote) {
    ftpfile <- "ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/xml/hla.xml.zip"
    dlfile <- tempfile(tmpdir = tdir)
    download.file(url = ftpfile, destfile = dlfile, method = "libcurl")
    tfile <- unzip(zipfile = dlfile, exdir = tdir)
    on.exit(unlink(tfile, force = TRUE))
  } else {
    assertive.properties::assert_is_not_null(db_path)
    dbfile <- tryCatch(
      normalizePath(file.path(db_path, "IMGTHLA", "xml", "hla.xml.zip"), mustWork = TRUE),
      error = function(e) normalizePath(file.path(db_path, "IMGTHLA", "xml", "hla.xml.gz"), mustWork = TRUE))
    if (endsWith(dbfile, "zip")) {
      tfile  <- unzip(zipfile = dbfile, exdir = tdir)[1]
      on.exit(unlink(tfile, force = TRUE))
    } else if (endsWith(dbfile, "gz")) {
      tfile  <- gzfile(dbfile, open = "b")
      on.exit(close(tfile))
    }

  }

  xml2::read_xml(tfile)
}

#' @rdname read_hla_xml
#' @return `invisible(NULL)`
#' @export
update_hla_xml <- function(db_path = getOption("hlatools.local_repos")) {
  assertive.properties::assert_is_not_null(db_path)
  dlfile <- normalizePath(file.path(db_path, "IMGTHLA", "xml", "hla.xml.zip"), mustWork = FALSE)
  ftpfile <- "ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/xml/hla.xml.zip"
  download.file(url = ftpfile, destfile = dlfile, method = "libcurl")
}

#' Parse all HLA alleles for a locus from hla.xml
#'
#' @param doc \file{hla.xml} as an [xml2::xml_document] object.
#' @param locusname One of <`HLA-A, HLA-B, HLA-C, HLA-DQB1, HLA-DRB1, HLA-DPB1`>
#' @param ncores The number of compute cores to use.
#'
#' @return A [HLAAllele-class] object.
#' @seealso [read_hla_xml()], [parse_hla_alleles()], [HLARanges-class],
#' [HLAAllele-class], [HLAGene][HLAGene_]
#' @export
#' @examples
#' \dontrun{
#' doc  <- read_hla_xml()
#' dpb1 <- parse_hla_alleles(doc, "HLA-DPB1")
#' }
parse_hla_alleles <- function(doc, locusname, ncores = parallel::detectCores()) {
  locusname <- match_hla_locus(locusname)
  ns <- xml2::xml_ns(doc)
  xpath1 <- paste0("/d1:alleles/d1:allele/d1:locus[@locusname='", locusname, "']/parent::node()")
  nodes1 <- xml2::xml_find_all(doc, xpath1, ns)
  xpath2 <- paste0(".//d1:releaseversions[not(starts-with(@releasestatus, 'Allele Deleted'))]/parent::node()")
  nodes2 <- xml2::xml_find_all(nodes1, xpath2, ns)
  #slen_nodeset <- length(nodes2)
  rs <- HLAAllele(nodes = nodes2, locusname = locusname, ncores = ncores)
  rs
}

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
      ##
      cit_idx <- xml2::xml_find_lgl(nodes, "boolean(d1:citations)", ns)
      pmids   <- rep(NA_character_, length(cit_idx))
      pmids[cit_idx] <- vapply(xml2::xml_find_all(nodes[cit_idx], "./d1:citations", ns), function(node) {
        colon(xml2::xml_attr(xml2::xml_children(node), "pubmed"))
      }, FUN.VALUE = character(1))
      ##
      smp_idx <- xml2::xml_find_lgl(nodes, "boolean(d1:sourcematerial/d1:samples)", ns)
      samples <- rep(NA_character_, length(smp_idx))
      samples[smp_idx] <- vapply(xml2::xml_find_all(nodes[smp_idx], "./d1:sourcematerial/d1:samples", ns), function(node) {
        colon(xml2::xml_attr(xml2::xml_children(node), "name"))
      }, FUN.VALUE = character(1))
      ##
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
        pmid          = pmids,
        ethnicity     = vapply(xml2::xml_find_all(nodes, "./d1:sourcematerial/d1:ethnicity", ns), function(node) {
          colon(xml2::xml_text(xml2::xml_children(node)))
        }, FUN.VALUE = character(1)),
        sample       = samples
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

