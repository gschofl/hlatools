% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HLARanges.R
\name{HLARanges}
\alias{HLARanges}
\title{Constructor for \linkS4class{HLARanges} objects}
\usage{
HLARanges(
  seqnames = S4Vectors::Rle(),
  ranges = IRanges::IRanges(),
  id = NA_character_,
  order = NA_integer_,
  type = NA_character_,
  status = NA_character_,
  frame = NA_integer_,
  ...
)
}
\arguments{
\item{seqnames}{\code{character} vector of the sequence name.}

\item{ranges}{\linkS4class{IRanges} object containing the ranges.}

\item{id}{\code{character} vector of feature IDs.}

\item{order}{\code{integer} vector of feature order.}

\item{type}{\code{character} vector of feature type.}

\item{status}{\code{character} vector of completeness status}

\item{frame}{\code{integer} vector of reading order.}

\item{...}{Arguments passed to \code{\link[GenomicRanges:GRanges-class]{GenomicRanges::GRanges()}}.}
}
\value{
A \linkS4class{HLARanges} object.
}
\description{
Constructor for \linkS4class{HLARanges} objects
}
\examples{
showClass("HLARanges")
}
\seealso{
\code{\link[=parse_hla_alleles]{parse_hla_alleles()}}, \linkS4class{HLAAllele}
}
