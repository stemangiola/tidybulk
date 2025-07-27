#' analyse gene over-representation with GSEA
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description test_gene_overrepresentation() takes as input a `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a `tbl` with the GSEA statistics
#'
#' @importFrom rlang enquo quo_is_missing
#' @importFrom magrittr not
#' @importFrom dplyr filter distinct pull mutate
#' @importFrom SummarizedExperiment rowData
#'
#'
#' @name test_gene_overrepresentation
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .sample The name of the sample column
#' @param .entrez The ENTREZ ID of the transcripts/genes
#' @param .do_test A boolean column name symbol. It indicates the transcript to check
#' @param species A character. For example, human or mouse. MSigDB uses the latin species names (e.g., \"Mus musculus\", \"Homo sapiens\")
#' @param gene_sets  A character vector. The subset of MSigDB datasets you want to test against (e.g. \"C2\"). If NULL all gene sets are used (suggested). This argument was added to avoid time overflow of the examples.
#'
#' @param gene_set DEPRECATED. Use gene_sets instead.
#'
#' @details This wrapper execute gene enrichment analyses of the dataset using a list of transcripts and GSEA.
#' This wrapper uses clusterProfiler (DOI: doi.org/10.1089/omi.2011.0118) on the back-end.
#'
#' # Get MSigDB data
#' msigdb_data = msigdbr::msigdbr(species = species)
#' 
#' # Filter for specific gene collections if provided
#' if (!is.null(gene_collections)) {
#'   msigdb_data = filter(msigdb_data, gs_collection %in% gene_collections)
#' }
#' 
#' # Process the data
#' msigdb_data |>
#'   nest(data = -gs_collection) |>
#'   mutate(test =
#'            map(
#'              data,
#'              ~ clusterProfiler::enricher(
#'                my_entrez_rank,
#'                TERM2GENE=.x |> select(gs_name, ncbi_gene),
#'                pvalueCutoff = 1
#'              ) |>
#'                as_tibble()
#'            ))
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
#' # se_mini = tidybulk::se_mini[!rowData(tidybulk::se_mini)$entrez |> is.na(),] |> aggregate_duplicates(.transcript = entrez)
#' # df_entrez = mutate(df_entrez, do_test = feature %in% c("TNFRSF4", "PLCH2", "PADI4", "PAX7"))
#'
#' \dontrun{
#' 	test_gene_overrepresentation(
#' 		df_entrez,
#' 		.sample = sample,
#' 		.entrez = entrez,
#' 		.do_test = do_test,
#' 		species="Homo sapiens",
#'    gene_sets =c("C2")
#' 	)
#' }
#'
#' @references
#' Mangiola, S., Molania, R., Dong, R., Doyle, M. A., & Papenfuss, A. T. (2021). tidybulk: an R tidy framework for modular transcriptomic data analysis. Genome Biology, 22(1), 42. doi:10.1186/s13059-020-02233-7
#'
#' Yu, G., Wang, L. G., Han, Y., & He, Q. Y. (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology, 16(5), 284-287. doi:10.1089/omi.2011.0118
#'
#' @docType methods
#' @rdname test_gene_overrepresentation-methods
#' @export
#'
#'
setGeneric("test_gene_overrepresentation", function(
        .data,
        .entrez,
        .do_test,
        species,
        .sample = NULL,
        gene_sets = NULL,
        gene_set = NULL  # DEPRECATED
) standardGeneric("test_gene_overrepresentation"))




# Set internal
.test_gene_overrepresentation_SE = 		function(.data,
                                              .entrez,
                                              .do_test,
                                              species,
                                              .sample = NULL,
                                              gene_sets = NULL,
                                              gene_set = NULL  # DEPRECATED
)	{
    
    # Comply with CRAN NOTES
    . = NULL
    
    # DEPRECATION OF reference function
    if (is_present(gene_set) & !is.null(gene_set)) {
        
        # Signal the deprecation to the user
        deprecate_warn("1.3.1", "tidybulk::.test_gene_overrepresentation(gene_set = )", details = "The argument gene_set is now deprecated please use gene_sets.")
        gene_sets = gene_set
    }
    
    # Get column names
    .do_test = enquo(.do_test)
    .entrez = enquo(.entrez)
    #
    # expr <- rlang::quo_get_expr(.do_test)
    # env <- quo_get_env(x)
    #
    
    # Check if entrez is set
    if(quo_is_missing(.entrez))
        stop("tidybulk says: the .entrez parameter appears to no be set")
    
    # Check column type
    if (
        .data |> 
        rowData() |> 
        as_tibble(rownames = f_(.data)$name) |> 
        mutate(my_do_test = !!.do_test) |> 
        pull(my_do_test) |> 
        is("logical") |> 
        not()
    )
        stop("tidybulk says: .do_test column must be logical (i.e., TRUE or FALSE)")
    
    # Check packages msigdbr
    # Check if package is installed, otherwise install
    check_and_install_packages("msigdbr")
    
    
    # Check is correct species name
    if(species %in% msigdbr::msigdbr_species()$species_name |> not())
        stop(sprintf("tidybulk says: wrong species name. MSigDB uses the latin species names (e.g., %s)", paste(msigdbr::msigdbr_species()$species_name, collapse=", ")))
    
    # # Check if missing entrez
    # if(.data %>% filter(!!.entrez %>% is.na) %>% nrow() %>% gt(0) ){
    # 	warning("tidybulk says: there are .entrez that are NA. Those will be removed")
    # 	.data = .data %>%	filter(!!.entrez %>% is.na %>% not())
    # }
    
      .data |>
    pivot_transcript() |>
    filter(!!.do_test) |>
    distinct(!!.entrez) |>
    pull(!!.entrez) |>
    entrez_over_to_gsea(species, gene_collections = gene_sets) |>
        
        # Add methods used
        memorise_methods_used(c("clusterProfiler", "msigdbr", "msigdb"), object_containing_methods = .data)
    
    
}
#' test_gene_overrepresentation
#'
#' @docType methods
#' @rdname test_gene_overrepresentation-methods
#'
#' @return A `SummarizedExperiment` object
setMethod("test_gene_overrepresentation",
          "SummarizedExperiment",
          .test_gene_overrepresentation_SE)

#' test_gene_overrepresentation
#'
#' @docType methods
#' @rdname test_gene_overrepresentation-methods
#'
#' @return A `RangedSummarizedExperiment` object
setMethod("test_gene_overrepresentation",
          "RangedSummarizedExperiment",
          .test_gene_overrepresentation_SE)


#' @details
#' This function plots the GSEA for gene overrepresentation
#'
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom stats p.adjust
entrez_over_to_gsea = function(my_entrez_rank, species, gene_collections  = NULL){
    
    # From the page
    # https://yulab-smu.github.io/clusterProfiler-book/chapter5.html
    
    # Check if package is installed, otherwise install
    check_and_install_packages(c("fastmatch", "clusterProfiler"))
    
    # Check msigdbr version
    if(packageVersion("msigdbr") < "10.0.0") {
        stop("tidybulk says: msigdbr version must be >= 10.0.0")
    }
    
    
    # Get gene sets signatures
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
                   )) |>
        select(-data) |>
        unnest(test) |>
        
        # Order
        arrange(`p.adjust`) |>
        
        # format transcripts
        mutate(entrez = strsplit(geneID, "/")) |>
        select(-geneID)
    
}


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

