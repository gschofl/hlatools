#' Class: nmdp_tbl
#'
#' Constructor for an <`nmdp_tbl`> object.
#'
#' @source [NMDP](https://bioinformatics.bethematchclinical.org/HLA-Resources/Allele-Codes/Allele-Code-Lists/Allele-Code-List-in-Numerical-Order/)
#' @return
#' A [tibble] of mappings between NMDP codes and allelic subtypes with the fields:
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
  tmp <- tempfile(fileext = ".zip")
  on.exit(unlink(tmp))
  url <- "https://hml.nmdp.org/mac/files/numer.v3.zip"
  utils::download.file(url, tmp)
  rs <- readr::read_delim(tmp, "\t",
                          col_names = c("code", "subtype"),
                          col_types = readr::cols(
                            readr::col_character(),
                            readr::col_character()
                          ),
                          skip = 3)
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
  x0 <- dplyr::tibble(f1 = hla_field1(alleles), f2 = hla_field2(alleles)) %>%
    dplyr::filter(grepl("^[A-Z]+$", f2))
  x1 <- dplyr::left_join(x0, nmdp_tbl, by = c("f2" = "code"))
  x2 <- dplyr::mutate(
    .data = x1,
    subtype = unlist(purrr::map2(f1, subtype, function(a, b) {
      if (!grepl(":", b, fixed = TRUE)) {
        slash(paste0(a, ":", strsplit1(b, split = "/", fixed = TRUE)))
      }
      else b
    })),
    code = paste0(f1, ":", f2)
  ) %>%
    dplyr::select(code, subtype)
  x2
}

#' Class: g_tbl
#'
#' Constructor for a <`g_tbl`> object.
#'
#' @param db_path <[character]>; location of local IPD-IMGT/HLA repository.
#' @param remote <[logical]>; if `TRUE` pull data from [hla.alleles.org](https://www.hla.alleles.org/wmda),
#' if `FALSE` retrieve data from `db_path`.
#' @source [hla.alleles.org](https://www.hla.alleles.org/nomenclature/g_groups.html).
#' @return
#' A [tibble] with the fields:
#' \itemize{
#'  \item "gene": The HLA gene.
#'  \item "code": The G group name.
#'  \item "subtype": The '/'-separated subtypes into which a G Code expands.
#' }
#' @seealso [nmdp_table], [allele_table]
#' @export
#' @examples
#' \dontrun{
#' g_codes <- g_table()
#' }
g_table <- function(db_path = getOption("hlatools.local_repos"), remote = FALSE) {
  datafile <- normalizePath(file.path(db_path, "IMGTHLA", "wmda", "hla_nom_g.txt"), mustWork = FALSE)
  if (!file.exists(datafile) || remote) {
    con <- base::url("http://hla.alleles.org/wmda/hla_nom_g.txt", open = "rb")
    on.exit(close(con))
    tryCatch(open(con), error = function(e) {
      stop("Trying to access http://hla.alleles.org/wmda/: ", e$message, call. = FALSE)
    })
    if (readr::read_lines(con, n_max = 1) != "# file: hla_nom_g.txt") {
      warning("Possibly malformed file \"hla_nom_g.txt\" ",
              "downloaded from http://hla.alleles.org/wmda/.", immediate. = TRUE)
    }
  } else {
    con <- base::file(datafile, "rb")
    on.exit(close(con))
  }
  col_names <- c("gene", "subtype", "code")
  col_types <- readr::cols(
    readr::col_character(),
    readr::col_character(),
    readr::col_character()
  )
  rs <- readr::read_delim(
    con, ";",  col_names = col_names, col_types = col_types, comment = "#"
  ) %>%
    dplyr::filter(endsWith(gene, "*")) %>%
    dplyr::mutate(gene = sub("*", "", gene, fixed = TRUE)) %>%
    dplyr::filter(!is.na(code)) %>%
    dplyr::select(gene, code, subtype)
  structure(rs, class = c("g_tbl", class(rs)))
}

#' @export
print.g_tbl <- function(x, ..., n = 5) {
  cat("G group codes: ", sep = "")
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
