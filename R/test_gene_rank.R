#' analyse gene rank with GSEA
#'
#' \lifecycle{maturing}
#'
#' @description test_gene_rank() takes as input a `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a `tbl` with the GSEA statistics
#'
#' @importFrom rlang enquo quo_is_missing
#' @importFrom magrittr not
#' @importFrom dplyr filter arrange select mutate
#' @importFrom SummarizedExperiment rowData
#'
#'
#' @name test_gene_rank
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .sample The name of the sample column
#' @param .entrez The ENTREZ ID of the transcripts/genes
#' @param .arrange_desc A column name of the column to arrange in decreasing order
#' @param species A character. For example, human or mouse. MSigDB uses the latin species names (e.g., \"Mus musculus\", \"Homo sapiens\")
#' @param gene_sets A character vector or a list. It can take one or more of the following built-in collections as a character vector: c("h", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "kegg_disease", "kegg_metabolism", "kegg_signaling"), to be used with EGSEA buildIdx. c1 is human specific. Alternatively, a list of user-supplied gene sets can be provided, to be used with EGSEA buildCustomIdx. In that case, each gene set is a character vector of Entrez IDs and the names of the list are the gene set names.
#'
#' @param gene_set DEPRECATED. Use gene_sets instead.
#'
#' @details This wrapper execute gene enrichment analyses of the dataset using a list of transcripts and GSEA.
#' This wrapper uses clusterProfiler (DOI: doi.org/10.1089/omi.2011.0118) on the back-end.
#'
#' Undelying method:
## Get gene sets signatures
#' msigdbr::msigdbr(species = species) %>%
#'
#'	# Filter specific gene_sets  if specified. This was introduced to speed up examples executionS
#'	when(
#'		!is.null(gene_sets ) ~ filter(., gs_collection %in% gene_sets ),
#'		~ (.)
#'	) |>
#'
#'	# Execute calculation
#'	nest(data = -gs_collection) |>
#'	mutate(fit =
#'				 	map(
#'				 		data,
#'				 		~ 	clusterProfiler::GSEA(
#'				 			my_entrez_rank,
#'				 			TERM2GENE=.x |> select(gs_name, ncbi_gene),
#'				 			pvalueCutoff = 1
#'				 		)
#'
#'				 	))
#'
#' @return A consistent object (to the input)
#'
#'
#'
#'
#' @examples
#'
#' print("Not run for build time.")
#'
#' \dontrun{
#'
#' df_entrez = tidybulk::se_mini
#' df_entrez = mutate(df_entrez, do_test = .feature %in% c("TNFRSF4", "PLCH2", "PADI4", "PAX7"))
#' df_entrez  = df_entrez |> test_differential_abundance(~ condition)
#'
#'
#'	test_gene_rank(
#'		df_entrez,
#' 		.sample = .sample,
#'		.entrez = entrez,
#' 		species="Homo sapiens",
#'    gene_sets =c("C2"),
#'  .arrange_desc = logFC
#' 	)
#' }
#'
#' @references
#' Mangiola, S., Molania, R., Dong, R., Doyle, M. A., & Papenfuss, A. T. (2021). tidybulk: an R tidy framework for modular transcriptomic data analysis. Genome Biology, 22(1), 42. doi:10.1186/s13059-020-02233-7
#'
#' Yu, G., Wang, L. G., Han, Y., & He, Q. Y. (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology, 16(5), 284-287. doi:10.1089/omi.2011.0118
#'
#' @docType methods
#' @rdname test_gene_rank-methods
#' @export
#'
#'
setGeneric("test_gene_rank", function(.data,
                                      .entrez,
                                      .arrange_desc,
                                      species,
                                      
                                      gene_sets  = NULL,
                                      gene_set = NULL  # DEPRECATED
)
  standardGeneric("test_gene_rank"))



# Set internal
.test_gene_rank_SE = 		function(.data,
                                .entrez,
                                .arrange_desc,
                                species,
                                
                                gene_sets = NULL,
                                gene_set = NULL  # DEPRECATED
)	{
  
  # Comply with CRAN NOTES
  . = NULL
  
  # DEPRECATION OF reference function
  if (is_present(gene_set) & !is.null(gene_set)) {
    
    # Signal the deprecation to the user
    deprecate_warn("1.3.1", "tidybulk::test_gene_rank(gene_set = )", details = "The argument gene_set is now deprecated please use gene_sets.")
    gene_sets = gene_set
    
  }
  
  # Get column names
  .arrange_desc = enquo(.arrange_desc)
  .entrez = enquo(.entrez)
  #
  # expr <- rlang::quo_get_expr(.do_test)
  # env <- quo_get_env(x)
  #
  
  # Check if entrez is set
  if(quo_is_missing(.entrez))
    stop("tidybulk says: the .entrez parameter appears to no be set")
  
  # Check packages msigdbr
  # Check if package is installed, otherwise install
  check_and_install_packages("msigdbr")
  
  
  # Check is correct species name
  if(species %in% msigdbr::msigdbr_species()$species_name |> not())
    stop(sprintf("tidybulk says: wrong species name. MSigDB uses the latin species names (e.g., %s)", paste(msigdbr::msigdbr_species()$species_name, collapse=", ")))
  
  .data |>
    pivot_transcript() |>
    arrange(desc(!!.arrange_desc)) |>
    select(!!.entrez, !!.arrange_desc) |>
    deframe() |>
    entrez_rank_to_gsea(species, gene_collections = gene_sets) |>
    
    # Add methods used. It is here and not in functions because I need the original .data
    memorise_methods_used(c("clusterProfiler", "enrichplot"), object_containing_methods = .data) |>
          when(
        gene_sets |> is("character") ~ (.) |> memorise_methods_used("msigdbr"),
      ~ (.)
    )
  
  
}

#' test_gene_rank
#'
#' @docType methods
#' @rdname test_gene_rank-methods
#'
#' @return A `SummarizedExperiment` object
setMethod("test_gene_rank",
          "SummarizedExperiment",
          .test_gene_rank_SE)

#' test_gene_rank
#'
#' @docType methods
#' @rdname test_gene_rank-methods
#'
#' @importFrom stringr str_replace
#'
#' @return A `RangedSummarizedExperiment` object
setMethod("test_gene_rank",
          "RangedSummarizedExperiment",
          .test_gene_rank_SE)





#' @details
#' This function plots the GSEA for gene overrepresentation
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom tibble rowid_to_column
#' @importFrom stats p.adjust
#' @importFrom purrr map
#'
entrez_rank_to_gsea = function(my_entrez_rank, species, gene_collections  = NULL){

	# From the page
	# https://yulab-smu.github.io/clusterProfiler-book/chapter5.html

	# Check if package is installed, otherwise install
  check_and_install_packages(c("fastmatch", "clusterProfiler", "enrichplot", "ggplot2"))

  # Get gene sets signatures
  if(is.null(gene_collections ) )
    my_gene_collection = msigdbr::msigdbr(species = species)
  else if(gene_collections |> is("character"))
    my_gene_collection = msigdbr::msigdbr(species = species) |>  filter( tolower(gs_collection) %in% tolower(gene_collections) )
  else if(gene_collections |> is("list"))
            my_gene_collection = tibble(gs_name=names(.), ncbi_gene = . ) |> unnest(ncbi_gene) |> mutate(gs_collection = "user_defined")
 else
   stop("tidybulk says: the gene sets should be either a character vector or a named list")


   my_gene_collection |>


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

					 	)) |>
			mutate(test =
					 	map(
					 		fit,
					 		~ .x |>
					 			# ggplot2::fortify(showCategory=Inf) %>%
					 			as_tibble() |>
					 			rowid_to_column(var = "idx_for_plotting")
					 			#%>%
					 			#	mutate(plot = future_imap(ID, ~ enrichplot::gseaplot2(fit, geneSetID = .y, title = .x)))

					 	)) |>
		select(-data)


}