#' analyse gene over-representation with GSEA
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description test_gene_overrepresentation() takes as input a `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a `tbl` with the GSEA statistics
#'
#' @importFrom rlang enquo
#' @importFrom rlang quo_is_missing
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
#' @docType methods
#' @rdname test_gene_overrepresentation-methods
#' @export
#'
#'
setGeneric("test_gene_overrepresentation", function(
    .data,
    .formula,
    .entrez,
    .abundance = NULL,
    contrasts = NULL,
    methods = c("camera", "roast", "safe", "gage", "padog", "globaltest", "ora"),
    gene_sets = c("h", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "kegg_disease", "kegg_metabolism", "kegg_signaling"),
    species,
    cores = 10,
    method = NULL,
    .contrasts = NULL
) standardGeneric("test_gene_overrepresentation"))




# Set internal
.test_gene_overrepresentation_SE = 		function(.data,
                                              .formula,
                                              
                                              .entrez,
                                              .abundance = NULL,
                                              contrasts = NULL,
                                              methods = c("camera" ,    "roast" ,     "safe",       "gage"  ,     "padog" ,     "globaltest",  "ora" ),
                                              gene_sets = c("h", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "kegg_disease", "kegg_metabolism", "kegg_signaling"),
                                              species,
                                              cores = 10,
                                              
                                              # DEPRECATED
                                              method = NULL,
                                              .contrasts = NULL
)	{
    
    # Fix NOTEs
    . = NULL
    
    # DEPRECATION OF reference function
    if (is_present(method) & !is.null(method)) {
        
        # Signal the deprecation to the user
        deprecate_warn("1.3.2", "tidybulk::test_gene_enrichment(method = )", details = "The argument method is now deprecated please use methods")
        methods = method
    }
    
    # DEPRECATION OF .constrasts
    if (is_present(.contrasts) & !is.null(.contrasts)) {
        
        # Signal the deprecation to the user
        deprecate_warn("1.7.4", "tidybulk::test_differential_abundance(.contrasts = )", details = "The argument .contrasts is now deprecated please use contrasts (without the dot).")
        
        contrasts = .contrasts
    }
    
    .entrez = enquo(.entrez)
    
    # Check that there are no entrez missing
    .data =
        .data %>%
        when(
            filter(., !!.entrez %>% is.na) %>% nrow() %>% gt(0) ~ {
                warning("tidybulk says: There are NA entrez IDs. Those genes will be filtered")
                filter(., !!.entrez %>% is.na %>% not())
            },
            ~ (.)
        )
    
    # Check if duplicated entrez
    if(rowData(.data)[,quo_name(.entrez)] %>% duplicated() %>% any())
        stop("tidybulk says: There are duplicated .entrez IDs. Please use aggregate_duplicates(.transcript = entrez).")
    
    # For use within when
    .my_data = .data
    
    
    # Comply with CRAN NOTES
    . = NULL
    
    # Check if at least two samples for each group
    if (.data %>%
        pivot_sample() %>%
        count(!!as.symbol(parse_formula(.formula))) %>%
        distinct(n) %>%
        pull(n) %>%
        min %>%
        st(2))
        stop("tidybulk says: You need at least two replicates for each condition for EGSEA to work")
    
    
    # Create design matrix
    design =	model.matrix(	object = .formula,	data = .data %>% colData() 	)
    
    # Print the design column names in case I want contrasts
    message(
        sprintf(
            "tidybulk says: The design column names are \"%s\"",
            design %>% colnames %>% paste(collapse = ", ")
        )
    )
    
    my_contrasts =
        contrasts %>%
        when(
            length(.) > 0 ~ limma::makeContrasts(contrasts = ., levels = design),
            ~ NULL
        )
    
    # Check if package is installed, otherwise install
    check_and_install_packages("EGSEA")
    
    if (!"EGSEA" %in% (.packages())) {
        stop("EGSEA package not loaded. Please run library(\"EGSEA\"). With this setup, EGSEA require manual loading, for technical reasons.")
    }
    
    dge =
        .data %>%
        assays() %>%
        as.list() %>%
        .[[1]] %>%
        as.matrix %>%
        
        # Change rownames to entrez
        when(
            quo_is_null(.entrez) %>% `!` ~ {
                x = (.)
                rownames(x) =
                    .my_data %>%
                    pivot_transcript() %>%
                    pull(!!.entrez)
                x
            },
            ~ (.)
        ) %>%
        
        # Filter missing entrez
        .[rownames(.) %>% is.na %>% not, ] %>%
        
        # # Make sure transcript names are adjacent
        # arrange(!!.entrez) %>%
        
        # select(!!.sample, !!.entrez, !!.abundance) %>%
        # spread(!!.sample,!!.abundance) %>%
        # as_matrix(rownames = !!.entrez) %>%
        edgeR::DGEList(counts = .)
    
    # Add gene ids for Interpret Results tables in report
    dge$genes = rownames(dge$counts)
    
    if (is.list(gene_sets)) {
        
        idx =  buildCustomIdx(geneIDs = rownames(dge), species = species, gsets=gene_sets)
        nonkegg_genesets = idx
        kegg_genesets = NULL
        
    } else {
        
        # Specify gene sets to include
        msig_all <- c("h", "c1", "c2", "c3", "c4", "c5", "c6", "c7")
        kegg_all <- c("kegg_disease", "kegg_metabolism", "kegg_signaling")
        
        # Record which collections used (kegg, msigdb) for bibliography
        collections_bib = c()
        
        # Identify any msigdb sets to be included
        msigdb.gsets <- gene_sets[gene_sets %in% msig_all]
        if (length(msigdb.gsets) >= 1) {
            collections_bib = c(collections_bib, "msigdb")
        }
        
        # Have to identify kegg sets to exclude for EGSEA
        kegg_to_exclude = kegg_all[!(kegg_all %in% gene_sets)]
        
        # If all 3 kegg sets are excluded then set to "all" as specifying the 3 names gives empty kegg object
        if (length(kegg_to_exclude) == 3) {
            kegg.exclude = "all"
        } else {
            kegg.exclude = kegg_to_exclude %>% str_replace("kegg_", "")
            collections_bib = c(collections_bib, "kegg")
        }
        
        
        idx =  buildIdx(entrezIDs = rownames(dge), species = species,  msigdb.gsets = msigdb.gsets,
                        kegg.exclude = kegg.exclude)
        
        # Due to a bug with kegg pathview overlays, this collection is run without report
        # https://support.bioconductor.org/p/122172/#122218
        
        kegg_genesets = idx[which(names(idx)=="kegg")]
        nonkegg_genesets = idx[which(names(idx)!="kegg")]
    }
    
    # Specify column to use to sort results in output table
    # If only one method is specified there is no med.rank column
    if (length(methods) == 1) {
        sort_column = "p.value"
    } else {
        sort_column = "med.rank"
    }
    
    
    if (length(nonkegg_genesets) != 0) {
        res =
            dge %>%
            
            # Calculate weights
            limma::voom(design, plot = FALSE) %>%
            
            # Execute EGSEA
            egsea(
                contrasts = my_contrasts,
                gs.annots = nonkegg_genesets,
                baseGSEAs = methods,
                sort.by = sort_column,
                num.threads = cores
            )
        
        gsea_web_page = "https://www.gsea-msigdb.org/gsea/msigdb/cards/%s.html"
        
        res_formatted_nonkegg =
            res@results %>%
            map2_dfr(
                (.) %>% names,
                ~ .x[[1]][[1]] %>%
                    as_tibble(rownames = "pathway") %>%
                    mutate(data_base = .y)
            ) %>%
            arrange(sort_column) %>%
            
            # Add webpage
            mutate(web_page = sprintf(gsea_web_page, pathway)) %>%
            select(data_base, pathway, web_page, sort_column, everything())
    }
    
    if (length(kegg_genesets) != 0) {
        message("tidybulk says: due to a bug in the call to KEGG database (http://supportupgrade.bioconductor.org/p/122172/#122218), the analysis for this database is run without report production.")
        
        res_kegg =
            dge %>%
            
            # Calculate weights
            limma::voom(design, plot = FALSE) %>%
            
            # Execute EGSEA
            egsea(
                contrasts = my_contrasts,
                gs.annots = kegg_genesets,
                baseGSEAs = methods,
                sort.by = sort_column,
                num.threads = cores,
                report = FALSE
            )
        
        res_formatted_kegg =
            res_kegg@results %>%
            map2_dfr(
                (.) %>% names,
                ~ .x[[1]][[1]] %>%
                    as_tibble(rownames = "pathway") %>%
                    mutate(data_base = .y)
            ) %>%
            arrange(sort_column) %>%
            select(data_base, pathway, everything())
        
    }
    
    # output tibble
    if (exists("res_formatted_nonkegg") & exists("res_formatted_kegg")) {
        out = bind_rows(res_formatted_nonkegg, res_formatted_kegg)
    } else if (exists("res_formatted_nonkegg")) {
        out = res_formatted_nonkegg
    } else {
        out = res_formatted_kegg
    }
    
    # add to bibliography
    if (exists("collections_bib")) {
        out %>% memorise_methods_used(c("egsea", collections_bib, methods))
    }
    
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

