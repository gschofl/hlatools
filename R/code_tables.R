#' Class: nmdp_tbl
#'
#' Constructor for an <`nmdp_tbl`> object.
#'
#' @source [NMDP](https://bioinformatics.bethematchclinical.org/HLA-Resources/Allele-Codes/Allele-Code-Lists/Allele-Code-List-in-Numerical-Order/)
#'
#' @return
#' A table of mappings between NMDP codes and allelic subtypes with the fields:
#' \itemize{
#'  \item "code": An NMDP code.
#'  \item "subtype": The '/'-separated subtypes into which an NMDP code expands.
#' }
#' @seealso [g_table], [allele_table]
#' @export
#' @examples
#' \dontrun{
#' nmdp_codes <- nmdp_table()
#' }
nmdp_table <- function() {
  u <- "https://bioinformatics.bethematchclinical.org/HLA/numeric.v3.zip"
  tmp <- tempfile(fileext = ".zip")
  on.exit(unlink(tmp))
  curl::curl_download(u, tmp)
  con <- unz(tmp, filename = "numer.v3.txt")
  rs <- dtplyr::tbl_dt(
    data.table::setDT(
      scan(con, what = list("character", "character"),
           nlines = -1, sep = "\t", skip = 3, quiet = TRUE)
    )
  )
  close(con)
  data.table::setnames(rs, names(rs), c("code", "subtype"))
  data.table::setkeyv(rs, "code")
  structure(rs, class = c("nmdp_tbl", class(rs)))
}

#' @export
print.nmdp_tbl <- function(x, ..., n = 5) {
  cat("NMDP codes: ", sep = "")
  NextMethod(n = n)
  cat("...\n", sep = "")
}

#' @export
generate_nmdp_lookup <- function(alleles, nmdp_tbl) {
  x <- data.table::data.table(f1 = hla_field1(alleles), f2 = hla_field2(alleles))
  x <- x[grepl("^[A-Z]+$", f2)]
  if (is.null(data.table::key(nmdp_tbl)) || data.table::key(nmdp_tbl) != "code") {
    data.table::setkeyv(nmdp_tbl, "code")
  }
  data.table::setkeyv(x, "f2")
  x <- nmdp_tbl[x]
  x[, `:=`(subtype, unlist(purrr::map2(f1, subtype, function(a, b) {
    if (!grepl(":", b, fixed = TRUE)) {
      slash(paste0(a, ":", strsplit(b, split = "/", fixed = TRUE)[[1]]))
    }
    else b
  })))]
  x[, `:=`(code, paste0(f1, ":", code))]
  x[, `:=`(f1, NULL)]
  data.table::setkeyv(x, "code")
  x
}

#' Class: g_tbl
#'
#' Constructor for a <`g_tbl`> object.
#'
#' @source [hla.alleles.org](http://hla.alleles.org/nomenclature/g_groups.html).
#'
#' @return
#' A table with the fields:
#' \itemize{
#'  \item "gene": A HLA gene.
#'  \item "code": The G Code.
#'  \item "subtype": The '/'-separated subtypes into which a G Code expands.
#' }
#' @seealso [nmdp_table], [allele_table]
#' @export
#' @examples
#' \dontrun{
#' g_codes <- g_table()
#' }
g_table <- function() {
  dbpath <- getOption("hlatools.local_repos")
  gfile <- normalizePath(file.path(dbpath, "IMGTHLA", "wmda", "hla_nom_g.txt"), mustWork = FALSE)
  if (!file.exists(gfile)) {
    con <- curl::curl("http://hla.alleles.org/wmda/hla_nom_g.txt")
    on.exit(close(con))
    tryCatch(open(con), error = function(e) {
      stop("Trying to access http://hla.alleles.org/wmda/: ", e$message, call. = FALSE)
    })
    if (readLines(con, n = 1) != "# file: hla_nom_g.txt") {
      warning("Possibly malformed file \"hla_nom_g.txt\" ",
              "downloaded from http://hla.alleles.org/wmda/.", immediate. = TRUE)
    }
  } else {
    con <- file(gfile, "r")
    on.exit(close(con))
  }
  rs <- read.csv(con, header = FALSE, colClasses = "character", sep = ";", comment.char = "#")
  rs <- dtplyr::tbl_dt(data.table::setDT(rs))
  data.table::setnames(rs, names(rs), c("gene", "subtype", "code"))
  rs <- rs[, gene := paste0("HLA-", sub("*", "", gene, fixed = TRUE))]
  rs <- rs[nzchar(rs[, code])][, .(gene, code, subtype)]
  data.table::setkeyv(rs, "gene")
  structure(rs, class = c("g_tbl", class(rs)))
}

#' @export
print.g_tbl <- function(x, ..., n = 5) {
  cat("G codes: ", sep = "")
  NextMethod(n = n)
  cat("...\n", sep = "")
}


#' Class: allele_tbl
#'
#' Constructor for a <`allele_tbl`> object.
#'
#' @param db_path <[character]>; location of local IPD-IMGT/HLA repository.
#' @param remote <[logical]>; if `TRUE` pull data from [hla.alleles.org](https://www.hla.alleles.org/wmda),
#' if `FALSE` retrieve data from `db_path`.
#' @source [hla.alleles.org](https://www.hla.alleles.org/wmda).
#' @return
#' A [tibble] with the fields:
#' \itemize{
#'  \item "gene": The HLA gene.
#'  \item "allele_name": HLA Allele name.
#'  \item "date_assigned": Date assigned.
#'  \item "date_deleted": Date deleted, if the name has now been abandoned or `NA`.
#'  \item "identical_to":  Allele that the deleted allele was shown to be identical to.
#'  \item "reason_for_deletion": Reason for the Allele to be deleted.
#' }
#' @seealso [nmdp_table], [g_table]
#' @export
#' @examples
#' \dontrun{
#' a_tbl <- allele_table()
#' }
allele_table <- function(db_path = getOption("hlatools.local_repos"), remote = FALSE) {
  datafile <- normalizePath(file.path(db_path, "IMGTHLA", "wmda", "hla_nom.txt"), mustWork = FALSE)
  if (!file.exists(datafile) || remote) {
    con <- base::url("http://hla.alleles.org/wmda/hla_nom.txt", open = "rb")
    on.exit(close(con))
    tryCatch(open(con), error = function(e) {
      stop("Trying to access http://hla.alleles.org/wmda/: ", e$message, call. = FALSE)
    })
    if (readr::read_lines(con, n_max = 1) != "# file: hla_nom.txt") {
      warning("Possibly malformed file \"hla_nom.txt\" ",
              "downloaded from http://hla.alleles.org/wmda/.", immediate. = TRUE)
    }
  } else {
    con <- base::file(datafile, "rb")
    on.exit(close(con))
  }
  col_names <- c("gene", "allele_name", "date_assigned", "date_deleted", "identical_to", "reason_for_deletion")
  col_types <- readr::cols(
    readr::col_character(),
    readr::col_character(),
    readr::col_date(format = "%Y%m%d"),
    readr::col_date(format = "%Y%m%d"),
    readr::col_character(),
    readr::col_character()
  )
  rs <- readr::read_delim(
    con, ";",  col_names = col_names, col_types = col_types, comment = "#"
  ) %>%
    dplyr::filter(endsWith(gene, "*")) %>%
    dplyr::mutate(gene = sub("*", "", gene, fixed = TRUE))
  structure(rs, class = c("allele_tbl", class(rs)))
}

#' @export
print.allele_tbl <- function(x, ..., n = 5) {
  cat("HLA allele codes: ", sep = "")
  NextMethod(n = n)
  cat("...\n", sep = "")
}




