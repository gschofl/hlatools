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
#' @param doc \file{hla.xml} as an XML document.
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
parse_hla_alleles <- function(doc, locusname, ncores = parallel::detectCores() - 2) {
  locusname <- match_hla_locus(locusname)
  dbv <- numeric_version(xml2::xml_attr(xml2::xml_find_all(doc, "//d1:alleles/d1:allele[1]/d1:releaseversions"), "currentrelease"))
  ## prior to release 3.26.0 all DRBs were lumped together
  if (dbv < "3.26.0" && startsWith(locusname, "HLA-DRB")) {
    locusname0 <- locusname
    locusname  <- "HLA-DRB"
  }
  ## MIC locus naming varies across IPD-IMGT/HLA versions:
  ## version <= 3.38.0           : HLA-MICA, HLA-MICB
  ## 3.39.0 <= version <= 3.40.0 : MICA, MICB
  ## 3.41.x <= version <= 3.42.0 : HLA-MICA, HLA-MICB
  ## version == 3.43.0           : (no MIC data in XML - IPD-IMGT/HLA bug!)
  ## version == 3.44.0           : HLA-MICA, HLA-MICB
  ## version >= 3.45.0           : MICA, MICB

  if (((package_version(dbv) >= "3.39.0" & package_version(dbv) <= "3.40.0") |
       package_version(dbv) >= "3.45.0") && startsWith(locusname, "HLA-MIC")) {
    locusname  <- sub("HLA-", "", locusname)
  }
  ns <- xml2::xml_ns(doc)
  if (dbv <= "3.55.0") {
    xpath1 <- paste0("/d1:alleles/d1:allele/d1:locus[@locusname='", locusname, "']/parent::node()")
  } else {
    xpath1 <- paste0("/d1:release/d1:alleles/d1:allele/d1:locus[@locusname='", locusname, "']/parent::node()")
  }
  nodes1 <- xml2::xml_find_all(doc, xpath1, ns)
  xpath2 <- paste0(".//d1:releaseversions[not(starts-with(@releasestatus,'Allele Deleted'))]/parent::node()")
  nodes2 <- xml2::xml_find_all(nodes1, xpath2, ns)
  xpath3 <- paste0(".//d1:releaseversions[not(starts-with(@releasestatus,'Deleted'))]/parent::node()")
  nodes3 <- xml2::xml_find_all(nodes2, xpath3, ns)
  if (locusname == "HLA-DRB") {
    xpath4 <- paste0("./self::node()[starts-with(@name,'", locusname0, "')]")
    nodes3 <- xml2::xml_find_all(nodes3, xpath4, ns)
    locusname <- locusname0
  }
  #slen_nodeset <- length(nodes2)
  rs <- HLAAllele(nodes = nodes3, locusname = locusname, ncores = ncores)
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
      xpath  <- "./d1:sequence/d1:nucsequence"
      nucseq <- Biostrings::DNAStringSet(xml2::xml_text(xml2::xml_find_all(nodes, xpath, ns)))
      nucseq_names <- xml2::xml_attr(nodes, "name")
      #message("L1: ", length(nucseq))
      #message("L2: ", length(nucseq_names))
      stopifnot(length(nucseq) == length(nucseq_names))
      names(nucseq) <- nucseq_names
      nucseq
    },
    # Parse metadata from an hla.xml allele node
    #
    # @param nodes A hla.xml allele nodeset.
    # @param locusname The locus name.
    #
    # @return A DataFrame object.
    # @keywords internal
    parse_metadata = function(nodes, locusname) {
      ns      <- xml2::xml_ns(nodes)
      ##
      cit_idx <- xml2::xml_find_lgl(nodes, "boolean(d1:citations)", ns)
      #message("L5: ", length(cit_idx))
      pmids   <- rep(NA_character_, length(cit_idx))
      pmids[cit_idx] <- vapply(xml2::xml_find_all(nodes[cit_idx], "./d1:citations", ns), function(node) {
        colon(xml2::xml_attr(xml2::xml_children(node), "pubmed"))
      }, FUN.VALUE = character(1))
      ##
      smp_idx <- xml2::xml_find_lgl(nodes, "boolean(d1:sourcematerial/d1:samples)", ns)
      #message("L6: ", length(smp_idx))
      samples <- rep(NA_character_, length(smp_idx))
      samples[smp_idx] <- vapply(xml2::xml_find_all(nodes[smp_idx], "./d1:sourcematerial/d1:samples", ns), function(node) {
        colon(xml2::xml_attr(xml2::xml_children(node), "name"))
      }, FUN.VALUE = character(1))
      ## With v3.52.0 "ethnicity" was changed to "ancestry"
      ancestry <- vapply(xml2::xml_find_all(nodes, "./d1:sourcematerial/d1:ancestry", ns), function(node) {
        colon(xml2::xml_text(xml2::xml_children(node)))
      }, FUN.VALUE = character(1))
      if (length(ancestry) == 0) {
        ancestry <-  vapply(xml2::xml_find_all(nodes, "./d1:sourcematerial/d1:ethnicity", ns), function(node) {
          colon(xml2::xml_text(xml2::xml_children(node)))
        }, FUN.VALUE = character(1))
      }
      ## Expected number of Exons and introns for completeness
      #nfeatures <- length(feature_orders_(locus = locusname))
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
        ## CWD status and Completeness (we consider as complete alleles for which
        ## a full set of nonUTR features is present)
        #
        #cwd_status    = xml2::xml_find_chr(nodes, "string(./d1:cwd_catalogue/@cwd_status)", ns),
        #complete      = xml2::xml_find_lgl(nodes,
        #  paste0("count(./d1:sequence/d1:feature[@featuretype=\"Exon\" or @featuretype=\"Intron\"])=", nfeatures), ns),
        ##
        ## Source (PubMed ID, Ancestry, Sample/Cellline)
        ##
        pmid          = pmids,
        ancestry      = ancestry,
        sample        = samples
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
      # message("L3: ", length(nodeset))
      # message("L4: ", length(nodeset))
      # stopifnot(length(nodeset) == length(nodeset))
      rs <- HLARangesList(parallel::mcMap(function(seqname, node) {
        #message(seqname, " => ", appendLF = FALSE)
        #message(length(node))
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
