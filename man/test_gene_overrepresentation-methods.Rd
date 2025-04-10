% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R, R/methods_SE.R
\docType{methods}
\name{test_gene_overrepresentation}
\alias{test_gene_overrepresentation}
\alias{test_gene_overrepresentation,spec_tbl_df-method}
\alias{test_gene_overrepresentation,tbl_df-method}
\alias{test_gene_overrepresentation,tidybulk-method}
\alias{test_gene_overrepresentation,SummarizedExperiment-method}
\alias{test_gene_overrepresentation,RangedSummarizedExperiment-method}
\title{analyse gene over-representation with GSEA}
\usage{
test_gene_overrepresentation(
  .data,
  .entrez,
  .do_test,
  species,
  .sample = NULL,
  gene_sets = NULL,
  gene_set = NULL
)

\S4method{test_gene_overrepresentation}{spec_tbl_df}(
  .data,
  .entrez,
  .do_test,
  species,
  .sample = NULL,
  gene_sets = NULL,
  gene_set = NULL
)

\S4method{test_gene_overrepresentation}{tbl_df}(
  .data,
  .entrez,
  .do_test,
  species,
  .sample = NULL,
  gene_sets = NULL,
  gene_set = NULL
)

\S4method{test_gene_overrepresentation}{tidybulk}(
  .data,
  .entrez,
  .do_test,
  species,
  .sample = NULL,
  gene_sets = NULL,
  gene_set = NULL
)

\S4method{test_gene_overrepresentation}{SummarizedExperiment}(
  .data,
  .entrez,
  .do_test,
  species,
  .sample = NULL,
  gene_sets = NULL,
  gene_set = NULL
)

\S4method{test_gene_overrepresentation}{RangedSummarizedExperiment}(
  .data,
  .entrez,
  .do_test,
  species,
  .sample = NULL,
  gene_sets = NULL,
  gene_set = NULL
)
}
\arguments{
\item{.data}{A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))}

\item{.entrez}{The ENTREZ ID of the transcripts/genes}

\item{.do_test}{A boolean column name symbol. It indicates the transcript to check}

\item{species}{A character. For example, human or mouse. MSigDB uses the latin species names (e.g., \"Mus musculus\", \"Homo sapiens\")}

\item{.sample}{The name of the sample column}

\item{gene_sets}{A character vector. The subset of MSigDB datasets you want to test against (e.g. \"C2\"). If NULL all gene sets are used (suggested). This argument was added to avoid time overflow of the examples.}

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
test_gene_overrepresentation() takes as input a `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a `tbl` with the GSEA statistics
}
\details{
`r lifecycle::badge("maturing")`

This wrapper execute gene enrichment analyses of the dataset using a list of transcripts and GSEA.
This wrapper uses clusterProfiler (DOI: doi.org/10.1089/omi.2011.0118) on the back-end.

# Get MSigDB data
msigdb_data = msigdbr::msigdbr(species = species)

# Filter for specific gene collections if provided
if (!is.null(gene_collections)) {
  msigdb_data = filter(msigdb_data, gs_collection %in% gene_collections)
}

# Process the data
msigdb_data |>
  nest(data = -gs_collection) |>
  mutate(test =
           map(
             data,
             ~ clusterProfiler::enricher(
               my_entrez_rank,
               TERM2GENE=.x |> select(gs_name, ncbi_gene),
               pvalueCutoff = 1
             ) |>
               as_tibble()
           ))
}
\examples{

print("Not run for build time.")

#se_mini = aggregate_duplicates(tidybulk::se_mini, .transcript = entrez)
#df_entrez = mutate(df_entrez, do_test = feature \%in\% c("TNFRSF4", "PLCH2", "PADI4", "PAX7"))

\dontrun{
	test_gene_overrepresentation(
		df_entrez,
		.sample = sample,
		.entrez = entrez,
		.do_test = do_test,
		species="Homo sapiens",
   gene_sets =c("C2")
	)
}

}
