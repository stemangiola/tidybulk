% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R, R/methods_SE.R
\docType{methods}
\name{test_gene_rank}
\alias{test_gene_rank}
\alias{test_gene_rank,spec_tbl_df-method}
\alias{test_gene_rank,tbl_df-method}
\alias{test_gene_rank,tidybulk-method}
\alias{test_gene_rank,SummarizedExperiment-method}
\alias{test_gene_rank,RangedSummarizedExperiment-method}
\title{analyse gene rank with GSEA}
\usage{
test_gene_rank(
  .data,
  .entrez,
  .arrange_desc,
  species,
  .sample = NULL,
  gene_sets = NULL,
  gene_set = NULL
)

\S4method{test_gene_rank}{spec_tbl_df}(
  .data,
  .entrez,
  .arrange_desc,
  species,
  .sample = NULL,
  gene_sets = c("h", "c1", "c2", "c3", "c4", "c5", "c6", "c7"),
  gene_set = NULL
)

\S4method{test_gene_rank}{tbl_df}(
  .data,
  .entrez,
  .arrange_desc,
  species,
  .sample = NULL,
  gene_sets = c("h", "c1", "c2", "c3", "c4", "c5", "c6", "c7"),
  gene_set = NULL
)

\S4method{test_gene_rank}{tidybulk}(
  .data,
  .entrez,
  .arrange_desc,
  species,
  .sample = NULL,
  gene_sets = c("h", "c1", "c2", "c3", "c4", "c5", "c6", "c7"),
  gene_set = NULL
)

\S4method{test_gene_rank}{SummarizedExperiment}(
  .data,
  .entrez,
  .arrange_desc,
  species,
  .sample = NULL,
  gene_sets = NULL,
  gene_set = NULL
)

\S4method{test_gene_rank}{RangedSummarizedExperiment}(
  .data,
  .entrez,
  .arrange_desc,
  species,
  .sample = NULL,
  gene_sets = NULL,
  gene_set = NULL
)
}
\arguments{
\item{.data}{A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))}

\item{.entrez}{The ENTREZ ID of the transcripts/genes}

\item{.arrange_desc}{A column name of the column to arrange in decreasing order}

\item{species}{A character. For example, human or mouse. MSigDB uses the latin species names (e.g., \"Mus musculus\", \"Homo sapiens\")}

\item{.sample}{The name of the sample column}

\item{gene_sets}{A character vector or a list. It can take one or more of the following built-in collections as a character vector: c("h", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "kegg_disease", "kegg_metabolism", "kegg_signaling"), to be used with EGSEA buildIdx. c1 is human specific. Alternatively, a list of user-supplied gene sets can be provided, to be used with EGSEA buildCustomIdx. In that case, each gene set is a character vector of Entrez IDs and the names of the list are the gene set names.}

\item{gene_set}{DEPRECATED. Use gene_sets instead.}
}
\value{
A consistent object (to the input)

A `spec_tbl_df` object

A `tbl_df` object

A `tidybulk` object

A `SummarizedExperiment` object

A `RangedSummarizedExperiment` object
}
\description{
test_gene_rank() takes as input a `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a `tbl` with the GSEA statistics
}
\details{
\lifecycle{maturing}

This wrapper execute gene enrichment analyses of the dataset using a list of transcripts and GSEA.
This wrapper uses clusterProfiler (DOI: doi.org/10.1089/omi.2011.0118) on the back-end.

Undelying method:
msigdbr::msigdbr(species = species) %>%

# Filter specific gene_sets  if specified. This was introduced to speed up examples executionS
when(
	!is.null(gene_sets ) ~ filter(., gs_collection %in% gene_sets ),
	~ (.)
) |>

# Execute calculation
nest(data = -gs_collection) |>
mutate(fit =
			 	map(
			 		data,
			 		~ 	clusterProfiler::GSEA(
			 			my_entrez_rank,
			 			TERM2GENE=.x |> select(gs_name, ncbi_gene),
			 			pvalueCutoff = 1
			 		)

			 	))
}
\examples{

print("Not run for build time.")

\dontrun{

df_entrez = tidybulk::se_mini
df_entrez = mutate(df_entrez, do_test = .feature \%in\% c("TNFRSF4", "PLCH2", "PADI4", "PAX7"))
df_entrez  = df_entrez |> test_differential_abundance(~ condition)


test_gene_rank(
	df_entrez,
		.sample = .sample,
	.entrez = entrez,
		species="Homo sapiens",
   gene_sets =c("C2"),
 .arrange_desc = logFC
	)
}

}
