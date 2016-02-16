#' @include utils.r
NULL

#' Get or set sequence data
#'
#' @param x A \code{\linkS4class{HLAAllele}} object.
#' @param ... Further arguments passed to methods.
#'
#' @return A DNAStringSet object.
#' @export
#' @examples
#' ##
setGeneric("sequences", signature = "x", function(x, ...) standardGeneric("sequences"))
#' @param value A \code{DNAStringSet} object.
#' @rdname sequences
setGeneric("sequences<-", signature = "x", function(x, ..., value) standardGeneric("sequences<-"))

#' Get or set feature data
#'
#' @param x A \code{\linkS4class{HLAAllele}} object.
#' @param ... Further arguments passed to methods.
#'
#' @return A HLAFeatureList object.
#' @export
#' @examples
#' ##
setGeneric("features", signature = "x", function(x, ...) standardGeneric("features"))
#' @param value A \code{HLAFeatureList} object.
#' @rdname features
setGeneric("features<-", signature = "x", function(x, ..., value) standardGeneric("features<-"))

#' Access allele IDs
#'
#' @param x A \code{\linkS4class{HLAAllele}} object.
#' @param ... Further arguments passed to methods.
#'
#' @return A character vector of allele IDs.
#' @export
#' @examples
#' \dontrun{
#' dpb1 <- parse_hla_alleles(read_hla_xml(), "DPB1")
#' allele_id(dpb1[is_complete(dpb1), ])
#' }
setGeneric("allele_id", signature = "x", function(x, ...) standardGeneric("allele_id"))

#' Access allele names
#'
#' @param x A \code{\linkS4class{HLAAllele}} object.
#' @param ... Further arguments passed to methods.
#'
#' @return A character vector of allele names.
#' @export
#' @examples
#' \dontrun{
#' dpb1 <- parse_hla_alleles(read_hla_xml(), "DPB1")
#' allele_name(dpb1[is_complete(dpb1), ])
#' }
setGeneric("allele_name", signature = "x", function(x, ...) standardGeneric("allele_name"))

#' Access CWD status
#'
#' @param x A \code{\linkS4class{HLAAllele}} object.
#' @param ... Further arguments passed to methods.
#'
#' @return A character vector of CWD status codes.
#' @export
#' @examples
#' \dontrun{
#' dpb1 <- parse_hla_alleles(read_hla_xml(), "DPB1")
#' cwd_status(dpb1[is_complete(dpb1), ])
#' }
setGeneric("cwd_status", signature = "x", function(x, ...) standardGeneric("cwd_status"))

#' Access ethnicity status
#'
#' @param x A \code{\linkS4class{HLAAllele}} object.
#' @param ... Further arguments passed to methods.
#'
#' @return A character vector of colon-separated sample ethnicities.
#' @export
#' @examples
#' \dontrun{
#' dpb1 <- parse_hla_alleles(read_hla_xml(), "DPB1")
#' etnicity(dpb1[is_complete(dpb1), ])
#' }
setGeneric("ethnicity", signature = "x", function(x, ...) standardGeneric("ethnicity"))

#' Access completeness status
#'
#' @param x A \code{\linkS4class{HLAAllele}} object.
#' @param ... Further arguments passed to methods.
#'
#' @return A logical vector.
#' @export
#' @examples
#' \dontrun{
#' dpb1 <- parse_hla_alleles(read_hla_xml(), "DPB1")
#' dpb1[is_complete(dpb1), ]
#' }
setGeneric("is_complete", signature = "x", function(x, ...) standardGeneric("is_complete"))

### -------------------------------------------------------------------------
### HLARanges class
###

#' Get feature id
#'
#' @param x A \code{\linkS4class{HLARanges}} object.
#' @param ... Further arguments passed to methods.
#'
#' @return A character vector.
#' @export
#' @examples
#' ##
setGeneric("getId", signature = "x", function(x, ...) standardGeneric("getId"))

#' Get feature order
#'
#' @param x A \code{\linkS4class{HLARanges}} object.
#' @param ... Further arguments passed to methods.
#'
#' @return An integer vector.
#' @export
#' @examples
#' ##
setGeneric("getOrder", signature = "x", function(x, ...) standardGeneric("getOrder"))

#' Get feature type
#'
#' @param x A \code{\linkS4class{HLARanges}} object.
#' @param ... Further arguments passed to methods.
#'
#' @return A character vector.
#' @export
#' @examples
#' ##
setGeneric("getType", signature = "x", function(x, ...) standardGeneric("getType"))

#' Get feature status
#'
#' @param x A \code{\linkS4class{HLARanges}} object.
#' @param ... Further arguments passed to methods.
#'
#' @return A character vector.
#' @export
#' @examples
#' ##
setGeneric("getStatus", signature = "x", function(x, ...) standardGeneric("getStatus"))

#' Get feature reading frame
#'
#' @param x A \code{\linkS4class{HLARanges}} object.
#' @param ... Further arguments passed to methods.
#'
#' @return An integer vector.
#' @export
#' @examples
#' ##
setGeneric("getFrame", signature = "x", function(x, ...) standardGeneric("getFrame"))

