#' @import methods
#' @include all-generics.R
NULL

#' Class `"HLARanges"`
#'
#' Extends [GRanges-class].
#'
#' @slot id character.
#' @slot order integer.
#' @slot type character.
#' @slot status character.
#' @slot frame integer.
#'
#' @importClassesFrom GenomicRanges GRanges GRangesList
#' @keywords classes internal
#' @seealso [read_hla_xml()], [parse_hla_alleles()], [HLAAllele-class], [HLAGene][HLAGene_]
#' @export
#' @examples
#' showClass("HLARanges")
setClass(
  Class = "HLARanges",
  contains = "GRanges",
  slots = representation(
    id       = "character",
    order    = "integer",
    type     = "character", # UTR, Exon, Intron
    status   = "character", # Complete, Partial, NA
    frame    = "integer"
  )
)

#' Constructor for [HLARanges-class] objects
#'
#' @param seqnames `character` vector of the sequence name.
#' @param ranges [IRanges-class] object containing the ranges.
#' @param id `character` vector of feature IDs.
#' @param order `integer` vector of feature order.
#' @param type `character` vector of feature type.
#' @param status `character` vector of completeness status
#' @param frame `integer` vector of reading order.
#' @param ... Arguments passed to [GenomicRanges::GRanges()].
#'
#' @return A [HLARanges-class] object.
#' @seealso [parse_hla_alleles()], [HLAAllele-class]
#' @export
#' @examples
#' showClass("HLARanges")
HLARanges <- function(seqnames = S4Vectors::Rle(), ranges = IRanges::IRanges(),
                      id = NA_character_, order = NA_integer_,
                      type = NA_character_, status = NA_character_,
                      frame = NA_integer_, ...) {
  len <- length(ranges)
  if (is(seqnames, "character")) {
    seqnames <- S4Vectors::Rle(seqnames, len)
  }
  gr <- GenomicRanges::GRanges(seqnames, ranges, strand = S4Vectors::Rle("+", len), ...)
  if (!all(i <- type %in% c("Exon", "Intron", "UTR")) && !is.na(type)) {
    stop("Unknown feature type(s) ", comma(sQuote(type[!i])))
  }
  new("HLARanges", gr,
      id     = recycle_vector(id, len),
      order  = recycle_vector(order, len),
      type   = recycle_vector(type, len),
      status = recycle_vector(status, len),
      frame  = recycle_vector(frame, len))
}

setMethod("getId", "HLARanges", function(x, ...) x@id)

setMethod("getOrder", "HLARanges", function(x, ...) x@order)

setMethod("getType", "HLARanges", function(x, ...) x@type)

setMethod("getStatus", "HLARanges", function(x, ...) x@status)

setMethod("getFrame", "HLARanges", function(x, ...) x@frame)

setMethod("merge", signature(x = "HLARanges", y = "HLARanges"), function(x, y, ...) {
  nmx <- names(x)
  nmy <- names(y)
  i <- which(nmx %in% nmy)
  j <- which(nmy %in% nmx)
  x[i] <- y[j]
  ranges(x) <- normalise_ranges(x)
  x
})

setMethod("show", "HLARanges", function(object) {
  callNextMethod()
})

### Re-exports ####

#' @importFrom GenomicRanges start
#' @export
GenomicRanges::start

#' @importFrom GenomicRanges end
#' @export
GenomicRanges::end

#' @importFrom GenomicRanges width
#' @export
GenomicRanges::width

#' @importFrom GenomicRanges seqnames
#' @export
GenomicRanges::seqnames

#' @importFrom GenomicRanges ranges
#' @export
GenomicRanges::ranges

#' @importFrom GenomicRanges "ranges<-"
#' @export
GenomicRanges::`ranges<-`

#' @importFrom GenomicRanges elementMetadata
#' @export
GenomicRanges::elementMetadata

#' @importFrom GenomicRanges "elementMetadata<-"
#' @export
GenomicRanges::`elementMetadata<-`

setMethod(GenomicRanges:::extraColumnSlotNames, "HLARanges", function(x) {
  c("id", "order", "type", "status", "frame")
})

#' Class `"HLARangesList"`
#'
#' List container for [HLARanges-class] objects.
#'
#' @importClassesFrom S4Vectors DataFrame List
#' @keywords classes
#' @export
#' @seealso [parse_hla_alleles()], [HLAAllele-class], [HLARanges-class]
#' @examples
#' showClass("HLARangesList")
setClass(
  Class = "HLARangesList",
  representation = representation("VIRTUAL"),
  prototype = prototype(elementType = "HLARanges"),
  contains = "List"
)

setClass(
  Class = "CompressedHLARangesList",
  representation(elementMetadata = "DataFrame"),
  prototype = prototype(unlistData = new("HLARanges")),
  contains = c("HLARangesList", "CompressedGRangesList")
)

#' Constructor for [HLARangesList-class] objects
#'
#' @param ... [HLARanges-class] objects.
#'
#' @return A [HLARangesList-class] object
#' @seealso [parse_hla_alleles()], [HLAAllele-class], [HLARanges-class]
#' @importFrom GenomicRanges merge
#' @export
#' @examples
#' showClass("HLARangesList")
HLARangesList <- function(...) {
  new("CompressedHLARangesList", GenomicRanges::GRangesList(..., compress = TRUE))
}

setMethod("getId", "HLARangesList", function(x, ...) {
  unlisted_x <- unlist(x)
  IRanges::relist(getId(unlisted_x), x)
})

setMethod("getOrder", "HLARangesList", function(x, ...) {
  unlisted_x <- unlist(x)
  IRanges::relist(getOrder(unlisted_x), x)
})

setMethod("getType", "HLARangesList", function(x, ...) {
  unlisted_x <- unlist(x)
  IRanges::relist(getType(unlisted_x), x)
})

setMethod("getStatus", "HLARangesList", function(x, ...) {
  unlisted_x <- unlist(x)
  IRanges::relist(getStatus(unlisted_x), x)
})

setMethod("getFrame", "HLARangesList", function(x, ...) {
  unlisted_x <- unlist(x)
  IRanges::relist(getFrame(unlisted_x), x)
})

HLARanges_to_string <- function(hr) {
  seqnm  <- as.character(unique(seqnames(hr)))
  seqnum <- seq_along(seqnm)
  seqnm2  <- paste0(paste0(seqnm, ".", seqnum), collapse = "|")
  seqnum2 <- ifelse(as.character(seqnames(hr)) %in% seqnm[1], 1, 2)
  fname  <- gsub("[ ']", "", names(hr))
  fname  <- ifelse(grepl("UTR$", fname), paste0(substr(fname, 2, 4), substr(fname, 1, 1)), fname)
  fstart <- start(hr)
  fend   <- end(hr)
  fid    <- getId(hr)
  ftype  <- toupper(substr(getType(hr), 1, 1))
  fstat  <- toupper(substr(getStatus(hr) %|na|% "", 1, 1))
  fframe <- getFrame(hr) %|na|% ""
  str <- paste0(
    paste0(seqnum2, ":", fname, ":", fstart, ":", fend, ":", fid, ":", ftype, ":", fstat, ":", fframe, ":"),
    collapse = "|")
  sprintf("%s~%s", seqnm2, str)
}

setAs("HLARanges", "character", function(from) HLARanges_to_string(from))

setAs("HLARangesList", "character", function(from) sapply(from, HLARanges_to_string))

#' @import data.table
string_to_HLARanges <- function(hstr) {
  s <- strsplit(hstr, "~", fixed = TRUE)[[1]]
  fstr1  <- strsplit(s[1], "|", fixed = TRUE)[[1]]
  seqmat <- data.table(do.call("rbind", strsplit(fstr1, "[.]")), key = "V2")
  fstr2  <- strsplit(s[2], "|", fixed = TRUE)[[1]]
  fmat <- data.table(do.call("rbind", strsplit(fstr2, ":")), key = "V1")
  m <- seqmat[fmat][, V2 := NULL]
  setnames(m, names(m), c("seqnames", "names", "start", "end", "id", "type", "status", "frame"))
  m[, names := paste0(substr(names, 1, nchar(names) - 1), " ", substr(names, nchar(names), nchar(names)))]
  m[, names := ifelse(grepl("^UTR", names), paste0(substr(names, nchar(names), nchar(names)), "' UTR"), names)]
  m[, order := as.integer(strsplitN(id, ".", 2, fixed = TRUE))]
  m[, type := ifelse(type == "E", "Exon", ifelse(type == "I", "Intron", ifelse(type == "U", "UTR", NA_character_)))]
  m[, status := ifelse(status == "C", "Complete", ifelse(type == "P", "Partial", NA_character_))]
  m[, frame := ifelse(nzchar(frame), frame, NA_character_)]
  setorder(m, "order")
  HLARanges(
    seqnames = S4Vectors::Rle(m$seqnames),
    ranges = IRanges::IRanges(
      start = as.integer(m$start),
      end = as.integer(m$end),
      names = m$names
    ),
    id = m$id,
    order = as.integer(m$order),
    type  = m$type,
    status = m$status,
    frame  = as.integer(m$frame)
  )
}

setAs("character", "HLARanges", function(from) string_to_HLARanges(from))

setAs("character", "HLARangesList", function(from) {
  HLARangesList(lapply(from, string_to_HLARanges))
})
