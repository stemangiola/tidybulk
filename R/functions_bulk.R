
#' Create tt object from tibble
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#' @param df A tibble
#' @param sample_column A character name of the sample column
#' @param transcript_column A character name of the gene/transcript name column
#' @param counts_column A character name of the read count column
#'
#' @return A tibble with an additional column
#'
#' @export
create_tt_from_tibble_bulk = function(input.df,
                                      sample_column,
                                      transcript_column,
                                      counts_column) {

  # Make col names
  sample_column = enquo(sample_column)
  transcript_column = enquo(transcript_column)
  counts_column = enquo(counts_column)

  input.df %>%

    # # Check input types
    # error_if_wrong_input(
    #   as.list(environment()),
    #   c("spec_tbl_df",  "quosure",  "quosure",  "quosure")
    # ) %>%

    # Add parameters attribute
    add_attr(
      list(
        sample_column = sample_column,
        transcript_column = transcript_column,
        counts_column = counts_column
      ),
      "parameters"
    ) %>%

    # Add class
    add_class("tt") %>%
    add_class("ttBulk")
}


#' Convert bam/sam files to a tidy gene transcript counts data frame
#'
#' @param file_names A character vector
#' @return A tibble of gene counts
#'
#' @export
create_tt_from_bam_sam_bulk <-
  function(file_names, genome = "hg38") {
    # This function uses Subread to count the gene features,
    # annotate gene features with symbols, and
    # convert the data frame to tibble format

    n_cores <- system("nproc", intern = TRUE) %>%
      as.integer() %>%
      `-`(2)

    file_names %>%

      # Run subread
      Rsubread::featureCounts(
        annot.inbuilt = genome,
        nthreads = n_cores,
        isPairedEnd = T,
        requireBothEndsMapped = T,
        checkFragLength = F,
        useMetaFeatures = T
      ) %>%

      # Anonymous function
      # input: Subread results
      # output edgeR::DGEList object
      {
        edgeR::DGEList(
          counts = (.)$counts,
          genes = (.)$annotation[, c("GeneID", "Length")],
          samples = (.) %$% stat %>% as_tibble() %>% gather(sample, temp,-Status) %>% spread(Status, temp)
        )
      } %>%

      # Anonymous function
      # input: edgeR::DGEList object
      # output: edgeR::DGEList object with added transcript symbol
      {
        dge <- (.)
        dge$genes$transcript <-
          AnnotationDbi:::mapIds(
            org.Hs.eg.db::org.Hs.eg.db,
            keys = as.character(dge$genes$GeneID),
            column = "transcript",
            keytype = "ENTREZID",
            multiVals = "first"
          )

        dge
      } %>%

      # Anonymous function
      # input: annotated edgeR::DGEList object
      # output: tibble
      {
        reduce(
          list(
            (.) %$% counts %>% as_tibble(rownames = "GeneID") %>% mutate(GeneID = GeneID %>% as.integer()) %>% gather(sample, `read count`,-GeneID),
            (.) %$% genes %>% select(GeneID, transcript) %>% as_tibble(),
            (.) %$% samples %>% as_tibble()
          ),
          left_join
        ) %>%
          rename(entrez = GeneID) %>%
          mutate(entrez = entrez %>% as.character())
      }
  }


#' Get count per million for TMM normalisation.
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#' @param df A tibble
#' @param cpm_threshold A real positive number
#' @return A tibble with an additional column
add_normalised_counts_bulk.get_cpm <- function(df,
                                               sample_column = `sample`,
                                               transcript_column = `transcript`,
                                               counts_column = `read count`,
                                               cpm_threshold = 0.5) {
  sample_column = enquo(sample_column)
  transcript_column = enquo(transcript_column)
  counts_column = enquo(counts_column)

  if (cpm_threshold < 0)
    stop("The parameter cpm_threshold must be > 0")

  # Adjust cpm_theshold based on the library size
  cpm_threshold <-
    cpm_threshold /
    (
      df %>%
        group_by(!!sample_column) %>%
        summarise(s = sum(!!counts_column)) %>%
        ungroup() %>%
        summarise(m = median(s)) %>%
        pull(m) /
        1e6
    )

  # Add cmp and cmp threshold to the data set, and return
  df %>%
    left_join(
      (.) %>%
        select(-contains("ct")) %>%
        select(!!transcript_column,!!sample_column,!!counts_column) %>%
        spread(!!sample_column,!!counts_column) %>%
        drop_na() %>%
        do(
          bind_cols(
            !!transcript_column := (.) %>% pull(!!transcript_column),
            tibble::as_tibble((.) %>% select(-!!transcript_column) %>% edgeR::cpm())
          )
        ) %>%
        gather(!!sample_column, cpm,-!!transcript_column) %>%
        mutate(cpm_threshold = cpm_threshold),
      by = c(quo_name(transcript_column), quo_name(sample_column))
    )

}


#' Drop lowly tanscribed genes for TMM normalization
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#'
#' @param df A tibble
#' @param cpm_threshold A real positive number
#' @param prop A number between 0 and 1
#' @param sample_column A character name of the sample column
#' @param transcript_column A character name of the gene/transcript name column
#' @param counts_column A character name of the read count column
#' @return A tibble filtered
add_normalised_counts_bulk.get_low_expressed <- function(df,
                                                         sample_column = `sample`,
                                                         transcript_column = `transcript`,
                                                         counts_column = `read count`,
                                                         cpm_threshold = 0.5,
                                                         prop = 3 / 4) {
  sample_column = enquo(sample_column)
  transcript_column = enquo(transcript_column)
  counts_column = enquo(counts_column)

  if (cpm_threshold < 0)
    stop("The parameter cpm_threshold must be > 0")
  if (prop < 0 |
      prop > 1)
    stop("The parameter prop must be between 0 and 1")

  df %>%

    # Prepare the data frame
    select(!!transcript_column,!!sample_column,!!counts_column) %>%
    tidyr::spread(!!sample_column,!!counts_column) %>%
    gather(!!sample_column,!!counts_column,-!!transcript_column) %>%
    group_by(!!transcript_column) %>%
    mutate(
      !!counts_column :=
        ifelse(
          !!counts_column %>% is.na(),!!counts_column %>% median(na.rm = T) %>% as.integer(),!!counts_column
        )
    ) %>%
    ungroup() %>%

    # Calculate cpm
    add_normalised_counts_bulk.get_cpm(!!sample_column,!!transcript_column,!!counts_column) %>%

    # Filter based on how many samples have a gene below the threshold
    mutate(`gene above threshold` = (cpm > cpm_threshold) %>% as.integer) %>%
    group_by(!!transcript_column) %>%
    summarise(n = `gene above threshold` %>% sum) %>%
    filter(n < (max(n) * !!prop) %>% floor) %>%

    # Pull information
    pull(!!transcript_column) %>%
    as.character()
}


#' Calculate the norm factor with calcNormFactor from limma
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#' @param df A tibble
#' @param reference A reference matrix, not sure if used anymore
#' @param cpm_threshold A real positive number
#' @param prop A number between 0 and 1
#' @return A list including the filtered data frame and the normalization factors
add_normalised_counts_bulk.calcNormFactor <- function(df,
                                                      reference = NULL,
                                                      cpm_threshold = 0.5,
                                                      prop = 3 / 4,
                                                      sample_column = `sample`,
                                                      transcript_column = `transcript`,
                                                      counts_column = `read count`) {
  sample_column = enquo(sample_column)
  transcript_column = enquo(transcript_column)
  counts_column = enquo(counts_column)

  error_if_log_transformed(df,!!counts_column)

  # Get list of low transcribed genes
  gene_to_exclude <-
    add_normalised_counts_bulk.get_low_expressed(
      df %>%
        filter(sample != "reference"),!!sample_column,!!transcript_column,!!counts_column,
      cpm_threshold = cpm_threshold,
      prop = prop
    )

  # Check if transcript after filtering is 0
  if (length(gene_to_exclude) == df %>%
      dplyr::distinct(!!transcript_column) %>%
      nrow()) {
    stop("The gene expression matrix has been filtered completely for lowly expressed genes")
  }

  # Get data frame for the higly transcribed transcripts
  df.filt <-
    df %>%
    dplyr::filter(!(!!transcript_column %in% gene_to_exclude)) %>%
    droplevels()



  # List of low abundant transcripts
  gene_to_exclude = gene_to_exclude

  # Normalised data set
  nf =
    tibble::tibble(
      # Sample factor
      sample = factor(levels(df.filt %>% pull(!!sample_column))),

      # normalised data frame
      nf = edgeR::calcNormFactors(
        df.filt %>%
          tidyr::spread(!!sample_column,!!counts_column) %>%
          tidyr::drop_na() %>%
          dplyr::select(-!!transcript_column),
        refColumn = which(reference == factor(levels(
          df.filt %>% pull(!!sample_column)
        ))),
        method = "TMM"
      )
    ) %>%

    # Add the statistics about the number of genes filtered
    dplyr::left_join(
      df.filt %>%
        dplyr::group_by(!!sample_column) %>%
        dplyr::summarise(tot_filt = sum(!!counts_column, na.rm = T)) %>%
        dplyr::mutate(!!sample_column := as.factor(as.character(!!sample_column))),
      by = quo_name(sample_column)
    )

  # Return
  list(gene_to_exclude= gene_to_exclude, nf = nf)
}

#' Get a tibble with normalised read counts using TMM
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr equals
#'
#'
#' @param input.df A tibble
#' @param reference A reference matrix, not sure if used anymore
#' @param cpm_threshold A real positive number
#' @param prop A number between 0 and 1
#' @param method A character string
#' @param sample_column A character name of the sample column
#' @param transcript_column A character name of the gene/transcript name column
#' @param counts_column A character name of the read count column
#' @return A tibble including additional columns
#'
#' @export
get_normalised_counts_bulk <- function(input.df,
                                       sample_column = NULL,
                                       transcript_column = NULL,
                                       counts_column = NULL,
                                       cpm_threshold = 0.5,
                                       prop = 3 / 4,
                                       method = "TMM",
                                       reference_selection_function = median) {

  # Get column names
  sample_column = enquo(sample_column)
  transcript_column = enquo(transcript_column)
  counts_column = enquo(counts_column)
  col_names = get_sample_transcript_counts(input.df, sample_column, transcript_column, counts_column)
  sample_column = col_names$sample_column
  transcript_column = col_names$transcript_column
  counts_column = col_names$counts_column

  # Check if package is installed, otherwise install
  if ("edgeR" %in% rownames(installed.packages()) == FALSE) {
    writeLines("Installing edgeR needed for differential transcript abundance analyses")
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("edgeR")
  }

  # Set column name for value normalised
  value_normalised = as.symbol(sprintf("%s normalised",  quo_name(counts_column)))

  # Reformat input data set
  df <-
    input.df %>%

    # # Check input types
    # error_if_wrong_input(
    #   as.list(environment())[-1],
    #   c("spec_tbl_df",  "quosure",  "quosure",  "quosure")
    # ) %>%

    # Stop if any counts is NA
    error_if_counts_is_na(!!counts_column) %>%

    # Stop if there are duplicated transcripts
    error_if_duplicated_genes(!!sample_column,!!transcript_column,!!counts_column) %>%

    # Rename
    dplyr::select(!!sample_column,!!transcript_column,!!counts_column) %>%
    #setNames(c("!!sample_column", "gene", "read count")) %>%

    # Set samples and genes as factors
    mutate(
      !!sample_column := factor(!!sample_column),!!transcript_column := factor(!!transcript_column)
    )

  # Get norm factor object
  reference <-
    df %>%
    group_by(!!sample_column) %>%
    summarise(sum = sum(!!counts_column)) %>%
    mutate(med = reference_selection_function(sum)) %>%
    mutate(diff = abs(sum - med)) %>%
    arrange(diff) %>%
    head(n = 1) %>%
    pull(!!sample_column) %>%
    as.character()

  nf_obj <-
    add_normalised_counts_bulk.calcNormFactor(
      df,
      reference,
      cpm_threshold,
      prop,
      sample_column = !!sample_column,
      transcript_column = !!transcript_column,
      counts_column = !!counts_column
    )

  # Calculate normalization factors
  nf <- nf_obj$nf %>%
    left_join(
      df %>%
        group_by(!!sample_column) %>%
        summarise(tot = sum(!!counts_column, na.rm = T)) %>%
        mutate(!!sample_column := as.factor(as.character(!!sample_column))),
      by = quo_name(sample_column)
    ) %>%
    mutate(multiplier =
             1 /
             (tot_filt * nf) *
             ((.) %>% filter(!!sample_column == reference) %>% pull(tot))) %>%

    # I have correct the strange behaviour of edgeR of reference
    # sample not being 1
    ifelse_pipe(
      "reference" %in% ((.) %>% pull(!!sample_column)),
      ~ .x %>%
        mutate(
          multiplier =
            multiplier /
            (.) %>%
            filter(!!sample_column == "reference") %>%
            pull(multiplier)
        )
    )

  # Return
  df %>%
    mutate(!!sample_column := as.factor(as.character(!!sample_column))) %>%
    left_join(nf, by = quo_name(sample_column)) %>%

    # Calculate normalised values
    mutate(!!value_normalised := !!counts_column * multiplier) %>%

    # Format df for join
    dplyr::select(!!sample_column,
                  !!transcript_column,
                  !!value_normalised,
                  everything()) %>%
    dplyr::mutate(filtered_out_low_counts = !!transcript_column %in% nf_obj$gene_to_exclude) %>%
    dplyr::select(-!!counts_column,-tot,-tot_filt) %>%
    dplyr::rename(TMM = nf) %>%
    #setNames(c("sample", "gene", sprintf("%s normalised", input.df %>% select(!!counts_column) %>% colnames), colnames(.)[4:ncol(.)])) %>%
    arrange(!!sample_column,!!transcript_column) %>%
    #dplyr::select(-!!sample_column,-!!transcript_column) %>%

    # Attach attributes
    add_attr(input.df %>% attr("parameters"), "parameters")

}

#' Add a tibble with normalised read counts using TMM
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#'
#' @param df A tibble
#' @param reference A reference matrix, not sure if used anymore
#' @param cpm_threshold A real positive number
#' @param prop A number between 0 and 1
#' @param method A character string
#' @param sample_column A character name of the sample column
#' @param transcript_column A character name of the gene/transcript name column
#' @param counts_column A character name of the read count column
#' @return A tibble including additional columns
#'
#' @export
add_normalised_counts_bulk <- function(input.df,
                                       sample_column = NULL,
                                       transcript_column = NULL,
                                       counts_column = NULL,
                                       cpm_threshold = 0.5,
                                       prop = 3 / 4,
                                       method = "TMM",
                                       reference_selection_function = median) {

  # Get column names
  sample_column = enquo(sample_column)
  transcript_column = enquo(transcript_column)
  counts_column = enquo(counts_column)
  col_names = get_sample_transcript_counts(input.df, sample_column, transcript_column, counts_column)
  sample_column = col_names$sample_column
  transcript_column = col_names$transcript_column
  counts_column = col_names$counts_column

  input.df %>%
    arrange(!!sample_column,!!transcript_column) %>%

    # Add normalised data set
    bind_cols(
      input.df %>%
        get_normalised_counts_bulk(
          sample_column = !!sample_column,
          transcript_column = !!transcript_column,
          counts_column = !!counts_column,
          cpm_threshold = cpm_threshold,
          prop = prop,
          method = method,
          reference_selection_function = reference_selection_function
        ) %>%
        select(-contains(quo_name(sample_column)), -contains(quo_name(transcript_column)))
    ) %>%

    # Attach attributes
    add_attr(input.df %>% attr("parameters"), "parameters")
}


#' Get differential transcription information to a tibble using edgeR. At the moment only one covariate is accepted
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr set_colnames
#'
#'
#'
#' @param df A tibble
#' @param formula a formula with no response variable, referring only to numeric variables
#' @param sample_column A character name of the sample column
#' @param transcript_column A character name of the gene/transcript name column
#' @param counts_column A character name of the read count column
#' @param design_column A character name of the covariate column
#'
#' @return A tibble with edgeR results
#'
#' @export
get_differential_transcript_abundance_bulk <- function(input.df,
                                                       formula,
                                                       sample_column = NULL,
                                                       transcript_column = NULL,
                                                       counts_column = NULL,
                                                       significance_threshold = 0.05) {

  # Get column names
  sample_column = enquo(sample_column)
  transcript_column = enquo(transcript_column)
  counts_column = enquo(counts_column)
  col_names = get_sample_transcript_counts(input.df, sample_column, transcript_column, counts_column)
  sample_column = col_names$sample_column
  transcript_column = col_names$transcript_column
  counts_column = col_names$counts_column

  # distinct_at is not released yet for dplyr, thus we have to use this trick
  df_for_edgeR <- input.df %>%

    # # Check input types
    # error_if_wrong_input(
    #   as.list(environment())[-1],
    #   c("spec_tbl_df","formula", "quosure",  "quosure",  "quosure", "numeric")
    # ) %>%

    # Stop if any counts is NA
    error_if_counts_is_na(!!counts_column) %>%

    # Stop if there are duplicated transcripts
    error_if_duplicated_genes(!!sample_column,!!transcript_column,!!counts_column) %>%

    # Prepare the data frame
    select(!!transcript_column,
           !!sample_column,
           !!counts_column,
           one_of(parse_formula(formula))) %>%
    distinct()

  # Check if at least two samples for each group
  if(
    df_for_edgeR %>%
    select(!!sample_column, one_of(parse_formula(formula))) %>%
    distinct %>%
    count(!!as.symbol(parse_formula(formula))) %>%
    distinct(n) %>%
    pull(1) %>%
    `==` (1)
  ) stop("You need at least two replicated for each condition for edgeR to work")

  # Create design matrix
  design =
    model.matrix(
      object = formula,
      data = df_for_edgeR %>% select(!!sample_column, one_of(parse_formula(formula))) %>% distinct %>% arrange(!!sample_column)
    ) %>%
    magrittr::set_colnames(
      c(
        "(Intercept)",
        (.) %>% colnames %>% `[` (-1)
      ))

  # Check if package is installed, otherwise install
  if ("edgeR" %in% rownames(installed.packages()) == FALSE) {
    writeLines("Installing edgeR needed for differential transcript abundance analyses")
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("edgeR")
  }

  df_for_edgeR.filt <-
    df_for_edgeR %>%
    select(!!transcript_column,!!sample_column,!!counts_column) %>%
    mutate(
      filtered_out_low_counts = !!transcript_column %in% add_normalised_counts_bulk.get_low_expressed(.,!!sample_column,!!transcript_column,!!counts_column)
    )

  df_for_edgeR.filt %>%
    filter(!filtered_out_low_counts) %>%
    select(!!transcript_column,!!sample_column,!!counts_column) %>%
    spread(!!sample_column,!!counts_column) %>%
    as_matrix(rownames = !!transcript_column) %>%
    edgeR::DGEList(counts = .) %>%
    edgeR::calcNormFactors(method = "TMM") %>%
    edgeR::estimateGLMCommonDisp(design) %>%
    edgeR::estimateGLMTagwiseDisp(design) %>%
    edgeR::glmFit(design) %>%
    edgeR::glmLRT(coef = 2) %>%
    edgeR::topTags(n = 999999) %$%
    table %>%
    as_tibble(rownames = quo_name(transcript_column)) %>%

    # Mark DE genes
    mutate(is_de = FDR < significance_threshold) %>%

    # Add filtering info
    full_join(
      df_for_edgeR.filt %>%
        select(!!transcript_column, filtered_out_low_counts) %>%
        distinct()
    ) %>%

    # Arrange
    arrange(FDR) %>%

    # Attach attributes
    add_attr(input.df %>% attr("parameters"), "parameters")
}

#' Add differential transcription information to a tibble using edgeR. At the moment only one covariate is accepted
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr set_colnames
#'
#'
#' @param input.df A tibble
#' @param formula a formula with no response variable, referring only to numeric variables
#' @param sample_column A character name of the sample column
#' @param transcript_column A character name of the gene/transcript name column
#' @param counts_column A character name of the read count column
#' @param design_column A character name of the covariate column
#'
#' @return A tibble with differential_transcript_abundance results
#'
#' @export
add_differential_transcript_abundance_bulk <- function(input.df,
                                                       formula,
                                                       sample_column = NULL,
                                                       transcript_column = NULL,
                                                       counts_column = NULL,
                                                       significance_threshold = 0.05) {

  # Get column names
  sample_column = enquo(sample_column)
  transcript_column = enquo(transcript_column)
  counts_column = enquo(counts_column)
  col_names = get_sample_transcript_counts(input.df, sample_column, transcript_column, counts_column)
  sample_column = col_names$sample_column
  transcript_column = col_names$transcript_column
  counts_column = col_names$counts_column

  input.df %>%
    left_join(
      (.) %>%
        get_differential_transcript_abundance_bulk(formula,
                                                   sample_column = !!sample_column,
                                                   transcript_column = !!transcript_column,
                                                   counts_column = !!counts_column,
                                                   significance_threshold = significance_threshold)
    ) %>%

    # Arrange
    arrange(FDR) %>%

    # Attach attributes
    add_attr(input.df %>% attr("parameters"), "parameters")
}

#' Get K-mean clusters to a tibble
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom stats kmeans
#'
#'
#' @param input.df A tibble
#' @param feature_column A character string. The column that is represents entities to cluster (i.e., normally samples)
#' @param elements_column A character string. The column that is used to calculate distance (i.e., normally genes)
#' @param counts_column   A character string. The column that contains the numeric value (i.e., normally read counts)
#' @param number_of_clusters A integer indicating how many clusters we are seeking
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#' @return A tibble with additional columns
#'
#' @export
get_clusters_kmeans_bulk <-
  function(input.df,
           value_column,
           number_of_clusters,
           elements_column = NULL,
           feature_column = NULL,
           of_samples = T,
           log_transform = T) {

    # Get column names
    elements_column = enquo(elements_column)
    feature_column = enquo(feature_column)
    col_names = get_elements_features(input.df, elements_column, feature_column, of_samples)
    elements_column = col_names$elements_column
    feature_column = col_names$feature_column

    value_column = enquo(value_column)

    input.df %>%

      # Through error if some read counts are NA
      error_if_counts_is_na(!!value_column) %>%

      # Prepare data frame
      distinct(!!feature_column,!!elements_column,!!value_column) %>%

      # Check if log tranfrom is needed
      ifelse_pipe(log_transform,
                  ~ .x %>% mutate(!!value_column := !!value_column %>%  `+`(1) %>%  log())) %>%

      # Prepare data frame for return
      spread(!!feature_column,!!value_column) %>%
      as_matrix(rownames = !!elements_column) %>%
      kmeans(centers = number_of_clusters, iter.max = 1000) %$%
      cluster %>%
      as.list() %>%
      as_tibble() %>%
      gather(!!elements_column, cluster) %>%
      mutate(cluster = cluster %>% as.factor()) %>%

      # Attach attributes
      add_attr(input.df %>% attr("parameters"), "parameters")
  }

#' Add K-mean clusters to a tibble
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom stats kmeans
#'
#'
#' @param input.df A tibble
#' @param feature_column A character string. The column that is represents entities to cluster (i.e., normally samples)
#' @param elements_column A character string. The column that is used to calculate distance (i.e., normally genes)
#' @param counts_column   A character string. The column that contains the numeric value (i.e., normally read counts)
#' @param number_of_clusters A integer indicating how many clusters we are seeking
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#' @return A tibble with additional columns
#'
#' @export
add_clusters_kmeans_bulk <-
  function(input.df,
           value_column,
           number_of_clusters,
           elements_column = NULL,
           feature_column = NULL,
           of_samples = T,
           log_transform = T) {

    # Get column names
    elements_column = enquo(elements_column)
    feature_column = enquo(feature_column)
    col_names = get_elements_features(input.df, elements_column, feature_column, of_samples)
    elements_column = col_names$elements_column
    feature_column = col_names$feature_column

    value_column = enquo(value_column)

    input.df %>%
      left_join(
        (.) %>%
          get_clusters_kmeans_bulk(
            value_column = !!value_column,
            number_of_clusters = number_of_clusters,
            elements_column = !!elements_column,
            feature_column = !!feature_column,
            log_transform = log_transform
          )
      ) %>%

      # Attach attributes
      add_attr(input.df %>% attr("parameters"), "parameters")
  }

#' Get dimensionality information to a tibble using MDS
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom purrr map_dfr
#'
#'
#' @param input.df A tibble
#' @param feature_column A character string. The column that is represents entities to cluster (i.e., normally genes)
#' @param elements_column A character string. The column that is used to calculate distance (i.e., normally samples)
#' @param counts_column   A character string. The column that contains the numeric value (i.e., normally read counts)
#' @param components A integer vector corresponding to principal components of interest (e.g., list(1:2, 3:4, 5:6))
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#' @return A tibble with additional columns
#'
#' @export
get_reduced_dimensions_MDS_bulk <-
  function(input.df,
           value_column,
           components = 1:2,
           elements_column = NULL,
           feature_column = NULL,
           of_samples = T,
           log_transform = T) {

    # Get column names
    elements_column = enquo(elements_column)
    feature_column = enquo(feature_column)
    col_names = get_elements_features(input.df, elements_column, feature_column, of_samples)
    elements_column = col_names$elements_column
    feature_column = col_names$feature_column

    value_column = enquo(value_column)

    # Convert components to components list
    if((length(components) %% 2) != 0 ) components = components %>% append(components[1])
    components_list = split(components, ceiling(seq_along(components)/2))

    # Loop over components list and calculate MDS. (I have to make this process more elegant)
    components_list %>%
      map_dfr(
        ~ input.df %>%

          # Through error if some read counts are NA
          error_if_counts_is_na(!!value_column) %>%

          # Filter lowly transcribed (I have to avoid the use of normalising function)
          add_normalised_counts_bulk(!!elements_column,!!feature_column,!!value_column) %>%
          filter(!filtered_out_low_counts) %>%
          distinct(!!feature_column,!!elements_column,!!value_column) %>%

          # Check if logtansform is needed
          ifelse_pipe(
            log_transform,
            ~ .x %>% mutate(!!value_column := !!value_column %>% `+`(1) %>%  log())
          ) %>%

          # Stop any column is not if not numeric or integer
          ifelse_pipe(
            (.) %>% select(!!value_column) %>% summarise_all(class) %>% `%in%`(c("numeric", "integer")) %>% `!`() %>% any(),
            ~ stop("value_column must be numerical or integer")
          ) %>%
          spread(!!elements_column,!!value_column) %>%
          as_matrix(rownames = !!feature_column, do_check = FALSE) %>%
          limma::plotMDS(dim.plot = .x, plot = F) %>%

          # Anonymous function
          # input: MDS object
          # output: tibble
          {
            tibble(names((.)$x), (.)$x, (.)$y) %>%
              rename(
                !!elements_column := `names((.)$x)`,
                !!as.symbol(.x[1]) := `(.)$x`,
                !!as.symbol(.x[2]) := `(.)$y`
              ) %>%
              gather(Component, `Component value`,-!!elements_column)
          }
      )  %>%
      distinct() %>%
      spread(Component, `Component value`) %>%
      setNames(c((.) %>% select(1) %>% colnames(),
                 paste("Dimension", (.) %>% select(-1) %>% colnames())
      )) %>%

      # Attach attributes
      add_attr(input.df %>% attr("parameters"), "parameters")
  }

#' Add dimensionality information to a tibble using MDS
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#'
#' @param input.df A tibble
#' @param feature_column A character string. The column that is represents entities to cluster (i.e., normally genes)
#' @param elements_column A character string. The column that is used to calculate distance (i.e., normally samples)
#' @param counts_column   A character string. The column that contains the numeric value (i.e., normally read counts)
#' @param components A integer vector corresponding to principal components of interest (e.g., list(1:2, 3:4, 5:6))
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#' @return A tibble with additional columns
#'
#' @export
add_reduced_dimensions_MDS_bulk <-
  function(input.df,
           value_column ,
           components = 1:2,
           elements_column = NULL,
           feature_column = NULL,
           of_samples = T,
           log_transform = T) {

    # Get column names
    elements_column = enquo(elements_column)
    feature_column = enquo(feature_column)
    col_names = get_elements_features(input.df, elements_column, feature_column, of_samples)
    elements_column = col_names$elements_column
    feature_column = col_names$feature_column

    value_column = enquo(value_column)

    input.df %>%
      left_join(
        (.) %>%
          get_reduced_dimensions_MDS_bulk(
            value_column = !!value_column,
            components = components,
            elements_column = !!elements_column,
            feature_column = !!feature_column,
            log_transform = log_transform
          ),
        by = quo_name(elements_column)
      ) %>%

      # Attach attributes
      add_attr(input.df %>% attr("parameters"), "parameters")
  }

#' Get principal component information to a tibble using PCA
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#'
#'
#' @param input.df A tibble
#' @param feature_column A character string. The column that is represents entities to cluster (i.e., normally genes)
#' @param elements_column A character string. The column that is used to calculate distance (i.e., normally samples)
#' @param counts_column   A character string. The column that contains the numeric value (i.e., normally read counts)
#' @param components An integer vector corresponding to principal components of interest (e.g., list(1:2, 3:4, 5:6))
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#' @return A tibble with additional columns
#'
#' @export
get_reduced_dimensions_PCA_bulk <-
  function(input.df,
           value_column ,
           components = 1:2,
           elements_column = NULL,
           feature_column = NULL,
           of_samples = T,
           log_transform = T) {

    # Get column names
    elements_column = enquo(elements_column)
    feature_column = enquo(feature_column)
    col_names = get_elements_features(input.df, elements_column, feature_column, of_samples)
    elements_column = col_names$elements_column
    feature_column = col_names$feature_column

    value_column = enquo(value_column)

    input.df %>%

      # Through error if some read counts are NA
      error_if_counts_is_na(!!value_column) %>%

      # Prepare data frame
      distinct(!!feature_column,!!elements_column,!!value_column) %>%

      # Check if logtansform is needed
      ifelse_pipe(log_transform,
                  ~ .x %>% mutate(!!value_column := !!value_column %>% `+`(1) %>%  log())) %>%

      # Stop any column is not if not numeric or integer
      ifelse_pipe(
        (.) %>% select(!!value_column) %>% summarise_all(class) %>% `%in%`(c("numeric", "integer")) %>% `!`() %>% any(),
        ~ stop("value_column must be numerical or integer")
      ) %>%

      spread(!!elements_column,!!value_column) %>%

      drop_na %>% # Is this necessary?

      # check that there are non-NA genes for enough samples
      ifelse2_pipe(# First condition
        (.) %>% nrow == 0,

        # Second condition
        (.) %>% nrow < 100,

        # First function
        ~ stop(
          "In calculating correlation there is no gene that have non NA values is all samples"
        ),

        # Second function
        ~ {
          warning(
            "
                  In calculating correlation there is < 100 genes that have non NA values is all samples.
                  The correlation calculation would not be reliable,
                  we suggest to partition the dataset for sample clusters.
                  "
          )
          .x
        }) %>%

      # Transform to matrix
      as_matrix(rownames = !!feature_column, do_check = FALSE) %>%

      # Calculate principal components
      prcomp(scale = TRUE) %>%

      # Anonymous function - Prints fraction of variance
      # input: PCA object
      # output: PCA object
      {
        writeLines("Fraction of variance explained by the selected principal components")

        (.) %$% sdev %>% `^` (2) %>% # Eigen value
          `/` (sum(.)) %>%
          `[` (components) %>%
          enframe() %>%
          select(-name) %>%
          rename(`Fraction of variance` = value) %>%
          mutate(PC = components) %>%
          print(n = 9999999)

        (.)

      } %$%

      # Parse the PCA results to a tibble
      rotation %>%
      as_tibble(rownames = quo_name(elements_column)) %>%
      select(!!elements_column, sprintf("PC%s", components))

  }

#' Add principal component information to a tibble using PCA
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#'
#' @param input.df A tibble
#' @param feature_column A character string. The column that is represents entities to cluster (i.e., normally genes)
#' @param elements_column A character string. The column that is used to calculate distance (i.e., normally samples)
#' @param counts_column   A character string. The column that contains the numeric value (i.e., normally read counts)
#' @param components An integer vector corresponding to principal components of interest (e.g., list(1:2, 3:4, 5:6))
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#' @return A tibble with additional columns
#'
#' @export
add_reduced_dimensions_PCA_bulk <-
  function(input.df,
           value_column ,
           components = 1:2,
           elements_column = NULL,
           feature_column = NULL,
           of_samples = T,
           log_transform = T) {

    # Get column names
    elements_column = enquo(elements_column)
    feature_column = enquo(feature_column)
    col_names = get_elements_features(input.df, elements_column, feature_column, of_samples)
    elements_column = col_names$elements_column
    feature_column = col_names$feature_column

    value_column = enquo(value_column)

    input.df %>%
      left_join(
        (.) %>%
          get_reduced_dimensions_PCA_bulk(
            value_column = !!value_column,
            components = components,
            elements_column = !!elements_column,
            feature_column = !!feature_column,
            log_transform = log_transform
          )
      )
  }

#' Get rotated dimensions of two principal components or MDS dimension of choice, of an angle
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang quo_is_null
#'
#'
#' @param input.df A tibble
#' @param rotation_degrees A real number between 0 and 360
#' @param dimension_1_column A character string. The column of the dimension 1
#' @param dimension_2_column   A character string. The column of the dimension 2
#' @return A tibble with additional rotated columns
#'
#' @export
get_rotated_dimensions =
  function(input.df,
           dimension_1_column,
           dimension_2_column,
           rotation_degrees,
           elements_column = NULL,
           of_samples = T,
           dimension_1_column_rotated = NULL,
           dimension_2_column_rotated = NULL) {

    # Get column names
    elements_column = enquo(elements_column)
    col_names = get_elements(input.df, elements_column)
    elements_column = col_names$elements_column

    # Parse other colnames
    dimension_1_column = enquo(dimension_1_column)
    dimension_2_column = enquo(dimension_2_column)
    dimension_1_column_rotated = enquo(dimension_1_column_rotated)
    dimension_2_column_rotated = enquo(dimension_2_column_rotated)

    if(
      input.df %>%
      distinct(
        !!elements_column,
        !!dimension_1_column,
        !!dimension_2_column) %>%
      count(
        !!elements_column,
        !!dimension_1_column,
        !!dimension_2_column) %>%
      pull(n) %>%
      max %>%
      `>` (1)
    )
    stop(sprintf("%s must be unique for each row for the calculation of rotation", quo_name(elements_column)))

    # Set default col names for rotated dimensions if not set
    if (quo_is_null(dimension_1_column_rotated))
      dimension_1_column_rotated = as.symbol(sprintf(
        "%s rotated %s",
        quo_name(dimension_1_column),
        rotation_degrees
      ))
    if (quo_is_null(dimension_2_column_rotated))
      dimension_2_column_rotated = as.symbol(sprintf(
        "%s rotated %s",
        quo_name(dimension_2_column),
        rotation_degrees
      ))

    # Function that rotates a 2D space of a arbitrary angle
    rotation = function(m, d) {
      r = d * pi / 180
      ((dplyr::bind_rows(
        c(`1` = cos(r), `2` = -sin(r)),
        c(`1` = sin(r), `2` = cos(r))
      ) %>% as_matrix) %*% m)
    }

    # Sanity check of the angle selected
    if (rotation_degrees %>% between(-360, 360) %>% `!`)
      stop("rotation_degrees must be between -360 and 360")

    # Return
    input.df %>%
      distinct(!!elements_column,
               !!dimension_1_column,
               !!dimension_2_column) %>%
      as_matrix(rownames = !!elements_column) %>% t %>%
      rotation(rotation_degrees) %>%
      as_tibble() %>%
      mutate(`rotated dimensions` =
               c(
                 quo_name(dimension_1_column_rotated),
                 quo_name(dimension_2_column_rotated)
               )) %>%
      gather(!!elements_column, value,-`rotated dimensions`) %>%
      spread(`rotated dimensions`, value)

  }

#' Add Rotated dimensions of two principal components or MDS dimensions, of an angle
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang quo_is_null
#'
#'
#' @param input.df A tibble
#' @param rotation_degrees A real number between 0 and 360
#' @param dimension_1_column A character string. The column of the dimension 1
#' @param dimension_2_column   A character string. The column of the dimension 2
#' @return A tibble with additional rotated columns
#'
#' @export
add_rotated_dimensions =
  function(input.df,
           dimension_1_column,
           dimension_2_column,
           rotation_degrees,
           elements_column = NULL,
           of_samples = T,
           dimension_1_column_rotated = NULL,
           dimension_2_column_rotated = NULL) {

    # Get column names
    elements_column = enquo(elements_column)
    col_names = get_elements(input.df, elements_column)
    elements_column = col_names$elements_column

    # Parse other colnames
    dimension_1_column = enquo(dimension_1_column)
    dimension_2_column = enquo(dimension_2_column)
    dimension_1_column_rotated = enquo(dimension_1_column_rotated)
    dimension_2_column_rotated = enquo(dimension_2_column_rotated)

    # Set default col names for rotated dimensions if not set
    if (quo_is_null(dimension_1_column_rotated))
      dimension_1_column_rotated = as.symbol(sprintf(
        "%s rotated %s",
        quo_name(dimension_1_column),
        rotation_degrees
      ))
    if (quo_is_null(dimension_2_column_rotated))
      dimension_2_column_rotated = as.symbol(sprintf(
        "%s rotated %s",
        quo_name(dimension_2_column),
        rotation_degrees
      ))

    input.df %>%
      left_join(
        (.) %>%
          get_rotated_dimensions(
            dimension_1_column = !!dimension_1_column,
            dimension_2_column = !!dimension_2_column,
            rotation_degrees = rotation_degrees,
            elements_column = !!elements_column,
            of_samples = of_samples,
            dimension_1_column_rotated = !!dimension_1_column_rotated,
            dimension_2_column_rotated = !!dimension_2_column_rotated
          ),
        by = quo_name(elements_column)
      )
  }

#' Aggregates multiple read counts from the same samples (e.g., from isoforms)
#' This function aggregates read counts over samples, concatenates other character columns, and averages other numeric columns
#'
#' @importFrom dplyr summarise_all
#' @importFrom dplyr bind_rows
#'
#' @param input.df A tibble
#' @param sample_column A character name of the sample column
#' @param transcript_column A character name of the gene/transcript name column
#' @param counts_column A character name of the read count column
#' @param aggregation_function A function for counts aggregation (e.g., sum)
#' @return A tibble with aggregated genes and annotation
#'
#' @export
aggregate_duplicated_transcripts_bulk =
  function(input.df,
           aggregation_function = sum,
           sample_column = NULL,
           transcript_column = NULL,
           counts_column = NULL,
           keep_integer = T) {

    # Get column names
    sample_column = enquo(sample_column)
    transcript_column = enquo(transcript_column)
    counts_column = enquo(counts_column)
    col_names = get_sample_transcript_counts(input.df, sample_column, transcript_column, counts_column)
    sample_column = col_names$sample_column
    transcript_column = col_names$transcript_column
    counts_column = col_names$counts_column

    # Robust paste function that preserves NAs
    paste3 <- function(..., sep = ", ") {
      L <- list(...)
      L <- lapply(L, function(x) {
        x[is.na(x)] <- ""
        x
      })
      ret <- gsub(paste0("(^", sep, "|", sep, "$)"),
                  "",
                  gsub(paste0(sep, sep), sep,
                       do.call(paste, c(
                         L, list(sep = sep)
                       ))))
      is.na(ret) <- ret == ""
      ret
    }

    # Through warning if there are logicals of factor in the data frame
    # because they cannot be merged if they are not unique
    if ((lapply(input.df, class) %>% unlist %in% c("logical", "factor")) %>% any) {
      warning("for aggregation fctors and logical columns were converted to character")
      writeLines("Converted to characters")
      lapply(input.df, class) %>% unlist %>% `[` (. %in% c("logical", "factor") %>% which) %>% print
    }

    # Select which are the numerical columns
    numerical_columns = input.df %>% ungroup() %>% select_if(is.numeric) %>% select(-!!counts_column) %>% colnames() %>% c("n_aggr")

    # ggregates read input.df over samples, concatenates other character columns, and averages other numeric columns
    input.df %>%

      # Through error if some read counts are NA
      error_if_counts_is_na(!!counts_column) %>%

      # transform logials and factors
      mutate_if(is.factor, as.character) %>%
      mutate_if(is.logical, as.character) %>%

      # Add the nuber of duplicates for each gene
      left_join(
        (.) %>% count(!!sample_column,!!transcript_column, name = "n_aggr"),
        by = c(quo_name(sample_column), quo_name(transcript_column))
      ) %>%

      # Anonymous function - binds the unique and the reduced genes,
      # in the way we have to reduce redundancy just for the duplicated genes
      # input: tibble
      # output tibble distinct
      {
        dplyr::bind_rows(
          # Unique symbols
          (.) %>%
            filter(n_aggr == 1),

          # Duplicated symbols
          (.) %>%
            filter(n_aggr > 1) %>%
            group_by(!!sample_column,!!transcript_column) %>%
            mutate(
              !!counts_column := !!counts_column %>% aggregation_function()
            ) %>%
            mutate_at(vars(numerical_columns), mean) %>%
            mutate_at(
              vars(-group_cols(),-!!counts_column,-!!numerical_columns),
              list( ~ paste3(unique(.), collapse = ", "))
            ) %>%
            distinct()
        )
      } %>%

      # Rename column of number of duplicates for each gene
      rename(`number of merged transcripts` = n_aggr)

  }

#' Drop redundant elements (e.g., samples) for which feature (e.g., genes) aboundances are correlated
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom widyr pairwise_cor
#'
#'
#' @param input.df A tibble
#' @param feature_column A character string. The column that is represents entities to cluster (i.e., normally genes)
#' @param elements_column A character string. The column that is used to calculate distance (i.e., normally samples)
#' @param counts_column   A character string. The column that contains the numeric value (i.e., normally read counts)
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#' @param correlation_threshold A real number between 0 and 1
#' @return A tibble with redundant elemens removed
#'
#' @export
drop_redundant_elements_through_correlation <- function(input.df,
                                                        value_column,
                                                        correlation_threshold = 0.9,
                                                        elements_column = NULL,
                                                        feature_column = NULL,
                                                        of_samples = T,
                                                        log_transform = F) {

  # Get column names
  elements_column = enquo(elements_column)
  feature_column = enquo(feature_column)
  col_names = get_elements_features(input.df, elements_column, feature_column, of_samples)
  elements_column = col_names$elements_column
  feature_column = col_names$feature_column

  value_column = enquo(value_column)

  # Get the redundant data frame
  input.df.correlated =
    input.df %>%

    # Stop if any counts is NA
    error_if_counts_is_na(!!value_column) %>%

    # Stop if there are duplicated transcripts
    error_if_duplicated_genes(!!elements_column,!!feature_column,!!value_column) %>%

    # Prepare the data frame
    select(!!feature_column,!!elements_column,!!value_column) %>%

    # Check if logtansform is needed
    ifelse_pipe(log_transform,
                ~ .x %>% mutate(!!value_column := !!value_column %>% `+`(1) %>%  log())) %>%
    distinct() %>%
    spread(!!elements_column,!!value_column) %>%
    drop_na() %>%

    # check that there are non-NA genes for enough samples
    ifelse2_pipe(# First condition
      (.) %>% nrow == 0,

      # Second condition
      (.) %>% nrow < 100,

      # First function
      ~ stop(
        "In calculating correlation there is no gene that have non NA values is all samples"
      ),

      # Second function
      ~ {
        warning(
          "
                  In calculating correlation there is < 100 genes that have non NA values is all samples.
                  The correlation calculation would not be reliable,
                  we suggest to partition the dataset for sample clusters.
                  "
        )
        .x
      }) %>%

    # Prepare the data frame
    gather(!!elements_column,!!value_column,-!!feature_column) %>%
    rename(rc := !!value_column,
           sample := !!elements_column,
           transcript := !!feature_column) %>% # Is rename necessary?
    mutate_if(is.factor, as.character) %>%

    # Run pairwise correlation and return a tibble
    widyr::pairwise_cor(
      sample,
      transcript,
      rc,
      sort = T,
      diag = FALSE,
      upper = F
    ) %>%
    filter(correlation > correlation_threshold) %>%
    distinct(item1) %>%
    rename(!!elements_column := item1)

  # Return non redudant data frame
  input.df %>% anti_join(input.df.correlated)
}

#' Identifies the closest pairs in a MDS contaxt and return one of them
#'
#' @param input.df A tibble
#' @param Dim_a_column A character string. The column of one principal component
#' @param Dim_b_column A character string. The column of another principal component
#' @param redundant_column A character string. The column that is represents entities to cluster (i.e., normally samples)
#' @param counts_column   A character string. The column that contains the numeric value (i.e., normally read counts)
#' @return A tibble with pairs dropped
#'
#' @export
drop_redundant_elements_though_reduced_dimensions <-
  function(input.df,
           Dim_a_column,
           Dim_b_column,
           elements_column = NULL,
           of_samples = T) {
    # This function identifies the closest pairs and return one of them

    # Get column names
    elements_column = enquo(elements_column)
    col_names = get_elements(input.df, elements_column)
    elements_column = col_names$elements_column

    Dim_a_column = enquo(Dim_a_column)
    Dim_b_column = enquo(Dim_b_column)

    # Find redundant samples
    input.df.redundant =

      # Calculate distances
      input.df %>%
      select(sample,!!Dim_a_column,!!Dim_b_column) %>%
      distinct() %>%
      as_matrix(rownames = !!redundant_column) %>%
      dist() %>%

      # Prepare matrix
      as.matrix() %>% as_tibble(rownames = "sample a") %>%
      gather(`sample b`, dist,-`sample a`) %>%
      filter(`sample a` != `sample b`) %>%

      # Sort the elements of the two columns to avoid eliminating all samples
      rowwise() %>%
      mutate(
        `sample 1` = c(`sample a`, `sample b`) %>% sort() %>% `[`(1),
        `sample 2` = c(`sample a`, `sample b`) %>% sort() %>% `[`(2)
      ) %>%
      ungroup() %>%
      select(`sample 1`, `sample 2`, dist) %>%
      distinct() %>%

      # Select closestpairs
      select_closest_pairs %>%

      # Select pair to keep
      select(1) %>%

      # Set the column names
      setNames(quo_name(elements_column))

    # Drop samples that are correlated with others and return
    input.df %>% anti_join(input.df.redundant)
  }

#' after wget, this function merges hg37 and hg38 mapping data bases - Do not execute!
#'
#' @return A tibble with ensembl-transcript mapping
#'
get_ensembl_symbol_mapping <- function() {
  # wget -O mapping_38.txt 'http://www.ensembl.org/biomart/martservice?query=  <Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6">  <Dataset name="hsapiens_gene_ensembl" interface="default"> <Attribute name="ensembl_transcript_id"/> <Attribute name="ensembl_gene_id"/><Attribute name="hgnc_symbol"/> </Dataset> </Query>'
  # wget -O mapping_37.txt 'http://grch37.ensembl.org/biomart/martservice?query=<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6"><Dataset name="hsapiens_gene_ensembl" interface="default"><Attribute name="ensembl_transcript_id"/><Attribute name="ensembl_gene_id"/><Attribute name="hgnc_symbol"/></Dataset></Query>'
  read_table2("~/third_party_sofware/ensembl_mapping/mapping_37.txt",
              col_names = F) %>%
    setNames(c("ensembl_transcript_id", "ensembl_gene_id", "hgnc_symbol")) %>%
    mutate(hg = "hg37") %>%
    bind_rows(
      read_table2(
        "~/third_party_sofware/ensembl_mapping/mapping_38.txt",
        col_names = F
      ) %>%
        setNames(
          c("ensembl_transcript_id", "ensembl_gene_id", "hgnc_symbol")
        ) %>%
        mutate(hg = "hg38")
    ) %>%
    drop_na() %>%
    select(-ensembl_transcript_id) %>%
    group_by(ensembl_gene_id) %>%
    arrange(hg %>% desc()) %>%
    slice(1) %>%
    ungroup() %>%
    {
      (.) %>% write_csv("~/third_party_sofware/ensembl_mapping/ensembl_symbol_mapping.csv")
      (.)
    }
}

#' Get transcript column from ensembl gene id
#'
#' @param input.df A tibble
#' @param ensembl_transcript_column A character string. The column that is represents ensembl gene id
#' @return A tibble with added annotation
#'
#' @export
get_symbol_from_ensembl <-
  function(input.df, ensembl_transcript_column) {
    ensembl_transcript_column = enquo(ensembl_transcript_column)

    input.df %>%
      select(!!ensembl_transcript_column) %>%
      distinct() %>%

      # Add name information
      left_join(
        ensembl_symbol_mapping %>%
          distinct(ensembl_gene_id, hgnc_symbol, hg) %>%
          rename(!!ensembl_transcript_column := ensembl_gene_id),
        by = "ens"
      )

  }

#' Add transcript column from ensembl gene id
#'
#' @param input.df A tibble
#' @param ensembl_transcript_column A character string. The column that is represents ensembl gene id
#' @return A tibble with added annotation
#'
#' @export
add_symbol_from_ensembl <-
  function(input.df, ensembl_transcript_column) {
    ensembl_transcript_column = enquo(ensembl_transcript_column)

    # Add new symbols column
    input.df %>%
      left_join((.) %>%
                  get_symbol_from_ensembl(!!ensembl_transcript_column)) %>%

      # Attach attributes
      add_attr(input.df %>% attr("parameters"), "parameters")
  }

#' Get cell type proportions from cibersort
#'
#' @import parallel
#'
#'
#' @param input.df A tibble
#' @param sample_column A character name of the sample column
#' @param transcript_column A character name of the gene/transcript name column
#' @param counts_column A character name of the read count column
#' @return A tibble including additional columns
#'
#' @export
get_cell_type_proportions = function(input.df,            sample_column = NULL,
                                     transcript_column = NULL,
                                     counts_column = NULL, ...) {

  # Get column names
  sample_column = enquo(sample_column)
  transcript_column = enquo(transcript_column)
  counts_column = enquo(counts_column)
  col_names = get_sample_transcript_counts(input.df, sample_column, transcript_column, counts_column)
  sample_column = col_names$sample_column
  transcript_column = col_names$transcript_column
  counts_column = col_names$counts_column

  # Check if package is installed, otherwise install
  if ("e1071" %in% rownames(installed.packages()) == FALSE) {
    writeLines("Installing e1071 needed for Cibersort")
    install.packages("e1071")
  }

  # Check if package is installed, otherwise install
  if ("preprocessCore" %in% rownames(installed.packages()) == FALSE) {
    writeLines("Installing preprocessCore needed for Cibersort")
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("preprocessCore")

  }

  # Load library which is optional for the whole package
  library(preprocessCore)

  # Check if there are enough genes for the signature
  if(
    (
      input.df %>%
      pull(!!transcript_column) %in% (X_cibersort %>% rownames)
    ) %>%
    which %>%
    length %>%
    `<` (50)
  ) stop("You have less than 50 genes that overlap the Cibersort signature. Please check again your input dataframe")

  input.df %>%

    # Check if some transcripts are duplicated
    error_if_duplicated_genes(!!sample_column,!!transcript_column,!!counts_column) %>%

    # Prepare data frame
    distinct(!!sample_column,!!transcript_column,!!counts_column) %>%
    spread(!!sample_column,!!counts_column) %>%
    data.frame(row.names = 1, check.names = F) %>%

    # Run Cibersort through custom function
    my_CIBERSORT(X_cibersort,	...) %$%
    proportions %>%

    # Parse results and return
    as_tibble(rownames = quo_name(sample_column)) %>%
    select(-`P-value`,-Correlation,-RMSE) %>%
    gather(`Cell type`, proportion,-!!sample_column) %>%

    # Attach attributes
    add_attr(input.df %>% attr("parameters"), "parameters")


}

#' Add cell type proportions from cibersort
#'
#' @import parallel
#'
#'
#' @param input.df A tibble
#' @param sample_column A character name of the sample column
#' @param transcript_column A character name of the gene/transcript name column
#' @param counts_column A character name of the read count column
#' @return A tibble including additional columns
#'
#' @export
add_cell_type_proportions = function(input.df,            sample_column = NULL,
                                     transcript_column = NULL,
                                     counts_column = NULL, ...) {

  # Get column names
  sample_column = enquo(sample_column)
  transcript_column = enquo(transcript_column)
  counts_column = enquo(counts_column)
  col_names = get_sample_transcript_counts(input.df, sample_column, transcript_column, counts_column)
  sample_column = col_names$sample_column
  transcript_column = col_names$transcript_column
  counts_column = col_names$counts_column

  input.df %>%

    # Add new annotation
    left_join(
      (.) %>%
        get_cell_type_proportions(
          sample_column = !!sample_column,
          transcript_column = !!transcript_column,
          counts_column = !!counts_column,
          ...
        ),
      by = "sample"
    )

}

#' Get adjusted read count for some batch effect
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr set_colnames
#'
#' @param df A tibble
#' @param formula a formula with no response variable, of the kind ~ factor_of_intrest + batch
#' @param sample_column A character name of the sample column
#' @param transcript_column A character name of the gene/transcript name column
#' @param counts_column A character name of the read count column
#' @param design_column A character name of the covariate column
#'
#' @return A tibble with adjusted counts
#'
#' @export
get_adjusted_counts_for_unwanted_variation_bulk <- function(input.df,
                                                            formula,
                                                            sample_column = NULL,
                                                            transcript_column = NULL,
                                                            counts_column = NULL,
                                                            log_transform = T) {

  # Get column names
  sample_column = enquo(sample_column)
  transcript_column = enquo(transcript_column)
  counts_column = enquo(counts_column)
  col_names = get_sample_transcript_counts(input.df, sample_column, transcript_column, counts_column)
  sample_column = col_names$sample_column
  transcript_column = col_names$transcript_column
  counts_column = col_names$counts_column

  # Check that formula includes at least two covariates
  if(parse_formula(formula) %>% length %>% `<` (2))
    stop("The formula must contain two covariates, the first being the factor of interest, the second being the factor of unwanted variation")

  # Check that formula includes no more than two covariates at the moment
  if(parse_formula(formula) %>% length %>% `>` (3))
    warning("Only the second covariate in the formula is adjusted for, at the moment")

  # New column name
  value_adjusted = as.symbol(sprintf("%s adjusted",  quo_name(counts_column)))

  # Stop is any counts are NAs
  input.df %>% error_if_counts_is_na(!!counts_column)

  df_for_combat <-
    input.df %>%
    select(!!transcript_column,
           !!sample_column,
           !!counts_column,
           one_of(parse_formula(formula))) %>%
    distinct() %>%

    # Check if logtansform is needed
    ifelse_pipe(log_transform,
                ~ .x %>% mutate(!!counts_column := !!counts_column %>% `+`(1) %>%  log())) %>%

    # Mark Filter low read counts
    mutate(
      filtered_out_low_counts =
        !!transcript_column %in%
        add_normalised_counts_bulk.get_low_expressed(.,!!sample_column,!!transcript_column,!!counts_column)
    )

  # Create design matrix
  design =
    model.matrix(
      object = as.formula("~" %>% paste0(parse_formula(formula)[1])),
      # get first argument of the formula
      data = df_for_combat %>% select(!!sample_column, one_of(parse_formula(formula))) %>% distinct %>% arrange(!!sample_column)
    ) %>%
    set_colnames(c("(Intercept)", parse_formula(formula)[1]))

  # Check if package is installed, otherwise install
  if ("sva" %in% rownames(installed.packages()) == FALSE) {
    writeLines("Installing sva - Combat needed for adjustment for unwanted variation")
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("sva")
  }

  df_for_combat %>%

    # Filter low read counts
    filter(!filtered_out_low_counts) %>%

    # Select relevant info
    distinct(!!transcript_column,!!sample_column,!!counts_column) %>%

    # Stop any column is not if not numeric or integer
    ifelse_pipe(
      (.) %>% select(!!counts_column) %>% summarise_all(class) %>% `%in%`(c("numeric", "integer")) %>% `!`() %>% any(),
      ~ stop("counts_column must be numerical or integer")
    ) %>%

    spread(!!sample_column,!!counts_column) %>%
    as_matrix(rownames = !!transcript_column,
              do_check = FALSE) %>%

    # Run combat
    sva::ComBat(
      batch =
        df_for_combat %>%
        distinct(!!sample_column,!!as.symbol(parse_formula(formula)[2])) %>%
        arrange(!!sample_column) %>%
        pull(2),
      mod = design
    ) %>%
    as_tibble(rownames = quo_name(transcript_column)) %>%
    gather(!!sample_column,!!counts_column,-!!transcript_column) %>%

    # ReverseLog transform if tranformed in the first place
    ifelse_pipe(
      log_transform,
      ~ .x %>%
        mutate(!!counts_column := !!counts_column %>% exp %>% `-`(1)) %>%
        mutate(
          !!counts_column := ifelse(!!counts_column < 0, 0,!!counts_column)
        ) %>%
        mutate(!!counts_column := !!counts_column %>% as.integer)
    ) %>%

    # Reset column names
    rename(!!value_adjusted := !!counts_column)  %>%


    # Add filtering info
    right_join(
      df_for_combat %>%
        distinct(
          !!transcript_column,
          !!sample_column,
          filtered_out_low_counts
        ),
      by = c(quo_name(transcript_column), quo_name(sample_column))
    )
}

#' Add adjusted read count for some batch effect
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr set_colnames
#'
#' @param df A tibble
#' @param formula a formula with no response variable, of the kind ~ factor_of_intrest + batch
#' @param sample_column A character name of the sample column
#' @param transcript_column A character name of the gene/transcript name column
#' @param counts_column A character name of the read count column
#' @param design_column A character name of the covariate column
#'
#' @return A tibble with adjusted counts
#'
#' @export
add_adjusted_counts_for_unwanted_variation_bulk <- function(input.df,
                                                            formula,
sample_column = NULL,
transcript_column = NULL,
counts_column = NULL,
                                                            log_transform = T) {

# Get column names
sample_column = enquo(sample_column)
transcript_column = enquo(transcript_column)
counts_column = enquo(counts_column)
col_names = get_sample_transcript_counts(input.df, sample_column, transcript_column, counts_column)
sample_column = col_names$sample_column
transcript_column = col_names$transcript_column
counts_column = col_names$counts_column

  input.df %>%

    # Add adjsted column
    left_join(
      (.) %>%
        get_adjusted_counts_for_unwanted_variation_bulk(
          formula,
          sample_column = !!sample_column,
          transcript_column = !!transcript_column,
          counts_column = !!counts_column,
          log_transform = log_transform
        ) ,
      by = c(quo_name(transcript_column), quo_name(sample_column))
    )

}
