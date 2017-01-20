#' Class: nmdp_tbl
#'
#' Constructor for an <\code{nmdp_tbl}> object.
#'
#' @source \url{https://bioinformatics.bethematchclinical.org/HLA-Resources/Allele-Codes/Allele-Code-Lists/Allele-Code-List-in-Numerical-Order/}.
#'
#' @return
#' A table of mappings between NMDP codes and allelic subtypes with the fields:
#' \itemize{
#'  \item "code": An NMDP code.
#'  \item "subtype": The '/'-separated subtypes into which an NMDP code expands.
#' }
#' @seealso \code{\link{g_table}}
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
#' Constructor for a <\code{g_tbl}> object.
#'
#' @source \url{http://hla.alleles.org/nomenclature/g_groups.html}.
#'
#' @return
#' A table with the fields:
#' \itemize{
#'  \item "gene": A HLA gene.
#'  \item "code": The G Code.
#'  \item "subtype": The '/'-separated subtypes into which a G Code expands.
#' }
#' @seealso \code{\link{nmdp_table}}
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
