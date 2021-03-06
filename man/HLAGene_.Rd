% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HLAGene.R
\name{HLAGene_}
\alias{HLAGene_}
\title{R6 Class \code{"HLAGene"}}
\description{
R6 Class \code{"HLAGene"}

R6 Class \code{"HLAGene"}
}
\details{
A container for the allelic information of an HLA gene.
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{HLAGene_$new()}}
\item \href{#method-print}{\code{HLAGene_$print()}}
\item \href{#method-get_hlatools_version}{\code{HLAGene_$get_hlatools_version()}}
\item \href{#method-get_db_version}{\code{HLAGene_$get_db_version()}}
\item \href{#method-get_locusname}{\code{HLAGene_$get_locusname()}}
\item \href{#method-get_alleles}{\code{HLAGene_$get_alleles()}}
\item \href{#method-has_distances}{\code{HLAGene_$has_distances()}}
\item \href{#method-calculate_exon_distance_matrix}{\code{HLAGene_$calculate_exon_distance_matrix()}}
\item \href{#method-get_closest_complete_neighbor}{\code{HLAGene_$get_closest_complete_neighbor()}}
\item \href{#method-get_reference_sequence}{\code{HLAGene_$get_reference_sequence()}}
\item \href{#method-get_all_reference_sequences}{\code{HLAGene_$get_all_reference_sequences()}}
\item \href{#method-get_reference_sequence_as_ft}{\code{HLAGene_$get_reference_sequence_as_ft()}}
\item \href{#method-get_all_reference_sequences_as_ft}{\code{HLAGene_$get_all_reference_sequences_as_ft()}}
\item \href{#method-get_extended_reference_set}{\code{HLAGene_$get_extended_reference_set()}}
\item \href{#method-clone}{\code{HLAGene_$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new HLAgene object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HLAGene_$new(
  locusname,
  db_version,
  db_path = getOption("hlatools.local_repos"),
  ncores = parallel::detectCores() - 2,
  with_dist = FALSE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{locusname}}{Name of locus.}

\item{\code{db_version}}{Version number of IPD-IMGT/HLA release or "Latest".}

\item{\code{db_path}}{Local path to IPD-IMGT/HLA.}

\item{\code{ncores}}{How many cores are used for parsing the data.}

\item{\code{with_dist}}{Calulate Levenshtein distances between all alleles.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{HLAgene} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
Display a HLAgene object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HLAGene_$print()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_hlatools_version"></a>}}
\if{latex}{\out{\hypertarget{method-get_hlatools_version}{}}}
\subsection{Method \code{get_hlatools_version()}}{
The package version under which a HLAGene object was created.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HLAGene_$get_hlatools_version()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_db_version"></a>}}
\if{latex}{\out{\hypertarget{method-get_db_version}{}}}
\subsection{Method \code{get_db_version()}}{
The IPD-IMGT/HLA database version.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HLAGene_$get_db_version()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_locusname"></a>}}
\if{latex}{\out{\hypertarget{method-get_locusname}{}}}
\subsection{Method \code{get_locusname()}}{
The name of the HLA locus.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HLAGene_$get_locusname()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_alleles"></a>}}
\if{latex}{\out{\hypertarget{method-get_alleles}{}}}
\subsection{Method \code{get_alleles()}}{
A \link{HLAAllele} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HLAGene_$get_alleles(allele)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{allele}}{The allele code.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-has_distances"></a>}}
\if{latex}{\out{\hypertarget{method-has_distances}{}}}
\subsection{Method \code{has_distances()}}{
Has a distance matrix between all alleles been calculated?
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HLAGene_$has_distances()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-calculate_exon_distance_matrix"></a>}}
\if{latex}{\out{\hypertarget{method-calculate_exon_distance_matrix}{}}}
\subsection{Method \code{calculate_exon_distance_matrix()}}{
Calculate a distance matrix between all alleles based on the available exon sequences.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HLAGene_$calculate_exon_distance_matrix()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_closest_complete_neighbor"></a>}}
\if{latex}{\out{\hypertarget{method-get_closest_complete_neighbor}{}}}
\subsection{Method \code{get_closest_complete_neighbor()}}{
Get the full-length allele that is closest to the query allele, based on the
available exon sequences.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HLAGene_$get_closest_complete_neighbor(allele, partially = TRUE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{allele}}{The allele code.}

\item{\code{partially}}{If \code{TRUE} match partial allele names.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_reference_sequence"></a>}}
\if{latex}{\out{\hypertarget{method-get_reference_sequence}{}}}
\subsection{Method \code{get_reference_sequence()}}{
Get the reference sequence for query allele extended to full length by
filling up missing sequence information with the closest complete neighbour.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HLAGene_$get_reference_sequence(allele)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{allele}}{The allele code.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_all_reference_sequences"></a>}}
\if{latex}{\out{\hypertarget{method-get_all_reference_sequences}{}}}
\subsection{Method \code{get_all_reference_sequences()}}{
Get the reference sequences for all allele extended to full length by
filling up missing sequence information with the closest complete neighbours.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HLAGene_$get_all_reference_sequences(allele)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{allele}}{The allele code.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_reference_sequence_as_ft"></a>}}
\if{latex}{\out{\hypertarget{method-get_reference_sequence_as_ft}{}}}
\subsection{Method \code{get_reference_sequence_as_ft()}}{
Get the reference sequence for query allele extended to full length by
filling up missing sequence information with the closest complete neighbour
formatted as a GenBank Feature Table.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HLAGene_$get_reference_sequence_as_ft(allele)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{allele}}{The allele code.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_all_reference_sequences_as_ft"></a>}}
\if{latex}{\out{\hypertarget{method-get_all_reference_sequences_as_ft}{}}}
\subsection{Method \code{get_all_reference_sequences_as_ft()}}{
Get the reference sequences for all allele extended to full length by
filling up missing sequence information with the closest complete neighbours
formatted as a GenBank Feature Table.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HLAGene_$get_all_reference_sequences_as_ft(allele)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{allele}}{The allele code.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_extended_reference_set"></a>}}
\if{latex}{\out{\hypertarget{method-get_extended_reference_set}{}}}
\subsection{Method \code{get_extended_reference_set()}}{
TODO: Documentation
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HLAGene_$get_extended_reference_set(verbose = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{verbose}}{Emit progress messages.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HLAGene_$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
