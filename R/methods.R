

# CREATE


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
create_ttBulk <- function(input.df,
                          sample_column,
                          transcript_column,
                          counts_column) {
  UseMethod("create_ttBulk", input.df)
}

#' @export
create_ttBulk.default <- function(input.df,
                                  sample_column,
                                  transcript_column,
                                  counts_column,
                                  cell_column = NULL,
                                  type,
                                  ...)
{
  print("This function cannot be applied to this object")
}

#' @export
create_ttBulk.tbl_df <- function(input.df,
                                 sample_column,
                                 transcript_column,
                                 counts_column,
                                 type,
                                 ...)
{
  # Make col names
  sample_column = enquo(sample_column)
  transcript_column = enquo(transcript_column)
  counts_column = enquo(counts_column)

  create_tt_from_tibble_bulk(input.df, !!sample_column, !!transcript_column,  !!counts_column)

}


# NORMALISE


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
#' @param reference_selection_function A fucntion median, or max for the fererence sample
#' @param action A character string. Whether to join the new information to the input tibble (add), or just get the non-redundant tibble with the new information (get).
#'
#' @return A tibble including additional columns
#'
#' @export
normalise_counts <- function(input.df,
                             sample_column = NULL,
                             transcript_column = NULL,
                             counts_column = NULL,
                             cpm_threshold = 0.5,
                             prop = 3 / 4,
                             method = "TMM",
                             reference_selection_function = median,
                             action = "add",
                             ...) {
  UseMethod("normalise_counts", input.df)
}

#' @export
normalise_counts.default <-  function(input.df,
                                      sample_column = NULL,
                                      transcript_column = NULL,
                                      counts_column = NULL,
                                      cpm_threshold = 0.5,
                                      prop = 3 / 4,
                                      method = "TMM",
                                      reference_selection_function = median,
                                      action = "add",
                                      ...)
{
  print("This function cannot be applied to this object")
}

#' @export
normalise_counts.tbl_df = normalise_counts.ttBulk <-
  function(input.df,
           sample_column = NULL,
           transcript_column = NULL,
           counts_column = NULL,
           cpm_threshold = 0.5,
           prop = 3 / 4,
           method = "TMM",
           reference_selection_function = median,
           action = "add")
  {
    # Make col names
    sample_column = enquo(sample_column)
    transcript_column = enquo(transcript_column)
    counts_column = enquo(counts_column)

    if (action == "add")
      add_normalised_counts_bulk(input.df,
                                 !!sample_column,
                                 !!transcript_column,
                                 !!counts_column)
    else if (action == "get")
      get_normalised_counts_bulk(input.df,
                                 !!sample_column,
                                 !!transcript_column,
                                 !!counts_column)
    else
      stop(
        "action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
      )
  }

# CLUSTER

#' Get clusters of elements (e.g., samples or transcripts)
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
#' @param value_column   A character string. The column that contains the numeric value (i.e., normally read counts)
#' @param number_of_clusters A integer indicating how many clusters we are seeking
#' @param of_samples if the input is tt object, this will indicate whether the element column will be sample or transcript column
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#' @param action A character string. Whether to join the new information to the input tibble (add), or just get the non-redundant tibble with the new information (get).
#'
#' @return A tibble with additional columns
#'
#' @export
annotate_clusters <- function(input.df,
                          value_column,
                          number_of_clusters,
                          elements_column = NULL,
                          feature_column = NULL,
                          of_samples = T,
                          log_transform = T,
                          action = "add") {
  UseMethod("annotate_clusters", input.df)
}

#' @export
annotate_clusters.default <-  function(input.df,
                                   value_column,
                                   number_of_clusters,
                                   elements_column = NULL,
                                   feature_column = NULL,
                                   of_samples = T,
                                   log_transform = T,
                                   action = "add")
{
  print("This function cannot be applied to this object")
}

#' @export
annotate_clusters.tbl_df = annotate_clusters.ttBulk <-  function(input.df,
                                                         value_column,
                                                         number_of_clusters,
                                                         elements_column = NULL,
                                                         feature_column = NULL,
                                                         method = "kmeans",
                                                         of_samples = T,
                                                         log_transform = T,
                                                         action = "add")
{
  # Make col names
  value_column = enquo(value_column)
  elements_column = enquo(elements_column)
  feature_column = enquo(feature_column)

  if(method == "kmeans"){
    if (action == "add")
      add_clusters_kmeans_bulk(
        input.df,
        value_column = !!value_column,
        number_of_clusters = number_of_clusters,
        elements_column = !!elements_column,
        feature_column = !!feature_column,
        log_transform = log_transform
      )
    else if (action == "get")
      get_clusters_kmeans_bulk(
        input.df,
        value_column = !!value_column,
        number_of_clusters = number_of_clusters,
        elements_column = !!elements_column,
        feature_column = !!feature_column,
        log_transform = log_transform
      )
    else
      stop(
        "action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
      )
  }
  else
    stop("the only supported method is \"kmeans\" ")

}

# REDUCE DIMENTIONS

#' Find reduced dimensions
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
#' @param value_column   A character string. The column that contains the numeric value (i.e., normally read counts)
#' @param method A character string
#' @param components A list of integer vectors corresponding to principal components of interest (e.g., list(1:2, 3:4, 5:6))
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#' @param action A character string. Whether to join the new information to the input tibble (add), or just get the non-redundant tibble with the new information (get).
#'
#' @return A tibble with additional columns
#'
#' @export
reduce_dimensions <- function(input.df,
                              value_column,
                              method,
                              components = 1:2,
                              elements_column = NULL,
                              feature_column = NULL,
                              of_samples = T,
                              log_transform = T,
                              action = "add") {
  UseMethod("reduce_dimensions", input.df)
}

#' @export
reduce_dimensions.default <-  function(input.df,
                                       value_column,
                                       method,
                                       components = 1:2,
                                       elements_column = NULL,
                                       feature_column = NULL,
                                       of_samples = T,
                                       log_transform = T,
                                       action = "add")
{
  print("This function cannot be applied to this object")
}

#' @export
reduce_dimensions.tbl_df = reduce_dimensions.ttBulk <-
  function(input.df,
           value_column,
           method,
           components = 1:2,
           elements_column = NULL,
           feature_column = NULL,
           of_samples = T,
           log_transform = T,
           action = "add")
  {
    # Make col names
    value_column = enquo(value_column)
    elements_column = enquo(elements_column)
    feature_column = enquo(feature_column)

    if (method == "MDS") {
      if (action == "add")
        add_reduced_dimensions_MDS_bulk(
          input.df,
          value_column = !!value_column,
          components = components,
          elements_column = !!elements_column,
          feature_column = !!feature_column,
          log_transform = log_transform
        )
      else if (action == "get")
        get_reduced_dimensions_MDS_bulk(
          input.df,
          value_column = !!value_column,
          components = components,
          elements_column = !!elements_column,
          feature_column = !!feature_column,
          log_transform = log_transform
        )
      else
        stop(
          "action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
        )
    }
    else if (method == "PCA") {
      if (action == "add")
        add_reduced_dimensions_PCA_bulk(
          input.df,
          value_column = !!value_column,
          components = components,
          elements_column = !!elements_column,
          feature_column = !!feature_column,
          log_transform = log_transform
        )
      else if (action == "get")
        get_reduced_dimensions_PCA_bulk(
          input.df,
          value_column = !!value_column,
          components = components,
          elements_column = !!elements_column,
          feature_column = !!feature_column,
          log_transform = log_transform
        )
      else
        stop(
          "action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
        )

    }
    else
      stop("method must be either \"MDS\" or \"PCA\"")

  }

# ROTATE DIMENTIONS

#' Calculate the rotated dimensions, of an arbitrary angle
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
#' @param dimension_2_column  A character string. The column of the dimension 2
#' @param elements_column A character string. The column that is used to calculate distance (i.e., normally genes)
#' @param of_samples if the input is tt object, this will indicate whether the element column will be sample or transcript column
#' @param dimension_1_column_rotated A character string. The column of the rotated dimension 1
#' @param dimension_2_column_rotated A character string. The column of the rotated dimension 2
#' @param action A character string. Whether to join the new information to the input tibble (add), or just get the non-redundant tibble with the new information (get).
#'
#' @return A tibble with additional rotated columns
#'
#' @export
rotate_dimensions <- function(input.df,
                              dimension_1_column,
                              dimension_2_column,
                              rotation_degrees,
                              elements_column = NULL,
                              of_samples = T,
                              dimension_1_column_rotated = NULL,
                              dimension_2_column_rotated = NULL,
                              action = "add") {
  UseMethod("rotate_dimensions", input.df)
}

#' @export
rotate_dimensions.default <-  function(input.df,
                                       dimension_1_column,
                                       dimension_2_column,
                                       rotation_degrees,
                                       elements_column = NULL,
                                       of_samples = T,
                                       dimension_1_column_rotated = NULL,
                                       dimension_2_column_rotated = NULL,
                                       action = "add")
{
  print("This function cannot be applied to this object")
}

#' @export
rotate_dimensions.tbl_df = rotate_dimensions.ttBulk <-
  function(input.df,
           dimension_1_column,
           dimension_2_column,
           rotation_degrees,
           elements_column = NULL,
           of_samples = T,
           dimension_1_column_rotated = NULL,
           dimension_2_column_rotated = NULL,
           action =
             "add")
  {
    # Make col names
    elements_column = enquo(elements_column)
    dimension_1_column = enquo(dimension_1_column)
    dimension_2_column = enquo(dimension_2_column)
    dimension_1_column_rotated = enquo(dimension_1_column_rotated)
    dimension_2_column_rotated = enquo(dimension_2_column_rotated)


    if (action == "add")
      add_rotated_dimensions(
        input.df,
        dimension_1_column = !!dimension_1_column,
        dimension_2_column = !!dimension_2_column,
        rotation_degrees = rotation_degrees,
        elements_column = !!elements_column,
        of_samples = of_samples,
        dimension_1_column_rotated = !!dimension_1_column_rotated,
        dimension_2_column_rotated = !!dimension_2_column_rotated
      )
    else if (action == "get")
      get_rotated_dimensions(
        input.df,
        dimension_1_column = !!dimension_1_column,
        dimension_2_column = !!dimension_2_column,
        rotation_degrees = rotation_degrees,
        elements_column = !!elements_column,
        of_samples = of_samples,
        dimension_1_column_rotated = !!dimension_1_column_rotated,
        dimension_2_column_rotated = !!dimension_2_column_rotated
      )
    else
      stop(
        "action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
      )
  }

# DROP REDUNDANT

#' Drop redundant elements (e.g., samples) for which feature (e.g., genes) aboundances are correlated
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom widyr pairwise_cor
#'
#'
#' @param input.df A tibble
#' @param method A string character. Method used for calculating the redundancy (e.g., correlation, reduced_dimensions proximity)
#' @param elements_column A character string. The column that is used to calculate distance (i.e., normally samples)
#' @param of_samples if the input is tt object, this will indicate whether the element column will be sample or transcript column

#' @param value_column  A character string. For correlation based calculation. The column that contains the numeric value (i.e., normally read counts)
#' @param feature_column A character string. For correlation based calculation. The column that is represents entities to cluster (i.e., normally genes)
#' @param correlation_threshold A real number between 0 and 1. For correlation based calculation.
#' @param log_transform A boolean. For correlation based calculation. Whether the value should be log-transformed (e.g., TRUE for RNA sequencing data).
#' @param Dim_a_column A character string. For reduced_dimension based calculation. The column of one principal component
#' @param Dim_b_column A character string. For reduced_dimension based calculation. The column of another principal component
#' @return A tibble with redundant elemens removed
#'
#' @export
drop_redundant <- function(input.df,
                           method,
                           elements_column = NULL,
                           of_samples = T,

                           value_column,
                           feature_column = NULL,
                           correlation_threshold = 0.9,
                           log_transform = F,

                           Dim_a_column,
                           Dim_b_column) {
  UseMethod("drop_redundant", input.df)
}

#' @export
drop_redundant.default <-  function(input.df,
                                    method,
                                    elements_column = NULL,
                                    of_samples = T,

                                    value_column,
                                    feature_column = NULL,
                                    correlation_threshold = 0.9,
                                    log_transform = F,

                                    Dim_a_column,
                                    Dim_b_column)
{
  print("This function cannot be applied to this object")
}

#' @export
drop_redundant.tbl_df = drop_redundant.ttBulk <-  function(input.df,
                                                           method,
                                                           elements_column = NULL,
                                                           of_samples = T,

                                                           value_column,
                                                           feature_column = NULL,
                                                           correlation_threshold = 0.9,
                                                           log_transform = F,

                                                           Dim_a_column,
                                                           Dim_b_column)
{
  # Make col names
  value_column = enquo(value_column)
  elements_column = enquo(elements_column)
  feature_column = enquo(feature_column)


  if (method == "correlation")
    drop_redundant_elements_through_correlation(
      input.df,
      value_column = !!value_column,
      elements_column = !!elements_column,
      feature_column = !!feature_column,
      correlation_threshold = correlation_threshold,
      of_samples = of_samples,
      log_transform = log_transform
    )
  else if (method == "reduced_dimensions")
    drop_redundant_elements_though_reduced_dimensions(
      input.df,
      Dim_a_column = !!Dim_a_column,
      Dim_b_column = !!Dim_b_column,
      elements_column = !!elements_column,
      of_samples = of_samples,
      log_transform = log_transform
    )
  else
    stop(
      "method must be either \"correlation\" for dropping correlated elements or \"reduced_dimension\" to drop the closest pair according to two dimensions (e.g., PCA)"
    )

}


# ADJUST


#' Get adjusted read count for unwanted variation
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
#' @param log_transform A boolean. For correlation based calculation. Whether the value should be log-transformed (e.g., TRUE for RNA sequencing data).
#' @param action A character string. Whether to join the new information to the input tibble (add), or just get the non-redundant tibble with the new information (get).
#'
#' @return A tibble with edgeR results
#'
#' @export
adjust_counts <- function(input.df,
                          formula,
                          sample_column = NULL,
                          transcript_column = NULL,
                          counts_column = NULL,
                          log_transform = T,
                          action = "add",
                          ...) {
  UseMethod("adjust_counts", input.df)
}

#' @export
adjust_counts.default <-  function(input.df,
                                   formula,
                                   sample_column = NULL,
                                   transcript_column = NULL,
                                   counts_column = NULL,
                                   log_transform = T,
                                   action = "add",
                                   ...)
{
  print("This function cannot be applied to this object")
}

#' @export
adjust_counts.tbl_df = adjust_counts.ttBulk <-  function(input.df,
                                                         formula,
                                                         sample_column = NULL,
                                                         transcript_column = NULL,
                                                         counts_column = NULL,
                                                         log_transform = T,
                                                         action =
                                                           "add")
{
  # Make col names
  sample_column = enquo(sample_column)
  transcript_column = enquo(transcript_column)
  counts_column = enquo(counts_column)

  if (action == "add")
    add_adjusted_counts_for_unwanted_variation_bulk(
      input.df,
      formula,
      sample_column = !!sample_column,
      transcript_column = !!transcript_column,
      counts_column = !!counts_column,
      log_transform = log_transform
    )
  else if (action == "get")
    get_adjusted_counts_for_unwanted_variation_bulk(
      input.df,
      formula,
      sample_column = !!sample_column,
      transcript_column = !!transcript_column,
      counts_column = !!counts_column,
      log_transform = log_transform
    )
  else
    stop(
      "action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
    )
}

# AGGREGATE DUPLICATES

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
#' @param keep_integer A boolean. Whether to force the aggregate counts to integer
#'
#' @return A tibble with aggregated genes and annotation
#'
#' @export
aggregate_duplicates <- function(input.df,
                                 aggregation_function = sum,
                                 sample_column = NULL,
                                 transcript_column = NULL,
                                 counts_column = NULL,
                                 keep_integer = T) {
  UseMethod("aggregate_duplicates", input.df)
}

#' @export
aggregate_duplicates.default <-  function(input.df,
                                          aggregation_function = sum,
                                          sample_column = NULL,
                                          transcript_column = NULL,
                                          counts_column = NULL,
                                          keep_integer = T)
{
  print("This function cannot be applied to this object")
}

#' @export
aggregate_duplicates.tbl_df = aggregate_duplicates.ttBulk <-
  function(input.df,
           aggregation_function = sum,
           sample_column = NULL,
           transcript_column = NULL,
           counts_column = NULL,
           keep_integer = T)  {
    # Make col names
    sample_column = enquo(sample_column)
    transcript_column = enquo(transcript_column)
    counts_column = enquo(counts_column)


    aggregate_duplicated_transcripts_bulk(
      input.df,
      aggregation_function = aggregation_function,
      sample_column = !!sample_column,
      transcript_column = !!transcript_column,
      counts_column = !!counts_column,
      keep_integer = T
    )
  }



# ANNOTATE cell-type


#' Get cell type proportions from cibersort
#'
#' @import parallel
#'
#'
#' @param input.df A tibble
#' @param sample_column A character name of the sample column
#' @param transcript_column A character name of the gene/transcript name column
#' @param counts_column A character name of the read count column
#' @param action A character string. Whether to join the new information to the input tibble (add), or just get the non-redundant tibble with the new information (get).
#'
#' @return A tibble including additional columns
#'
#' @export
annotate_cell_type <- function(input.df,
                               sample_column = NULL,
                               transcript_column = NULL,
                               counts_column = NULL,
                               action = "add",
                               ...) {
  UseMethod("annotate_cell_type", input.df)
}

#' @export
annotate_cell_type.default <-  function(input.df,
                                        sample_column = NULL,
                                        transcript_column = NULL,
                                        counts_column = NULL,
                                        action = "add",
                                        ...)
{
  print("This function cannot be applied to this object")
}

#' @export
annotate_cell_type.tbl_df = annotate_cell_type.ttBulk <-
  function(input.df,
           sample_column = NULL,
           transcript_column = NULL,
           counts_column = NULL,
           action = "add", ...)  {
    # Make col names
    sample_column = enquo(sample_column)
    transcript_column = enquo(transcript_column)
    counts_column = enquo(counts_column)

    if (action == "add")
      add_cell_type_proportions(
        input.df,
        sample_column = !!sample_column,
        transcript_column = !!transcript_column,
        counts_column = !!counts_column,
        ...
      )
    else if (action == "get")
      get_cell_type_proportions(
        input.df,
        sample_column = !!sample_column,
        transcript_column = !!transcript_column,
        counts_column = !!counts_column,
        ...
      )
    else
      stop(
        "action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
      )
  }

# ANNOTATE GENE SYMBOL

#' Add transcript column from ensembl gene id
#'
#' @param input.df A tibble
#' @param ensembl_transcript_column A character string. The column that is represents ensembl gene id
#' @param action A character string. Whether to join the new information to the input tibble (add), or just get the non-redundant tibble with the new information (get).
#'
#' @return A tibble with added annotation
#'
#' @export
annotate_symbol <- function(input.df,
                            ensembl_transcript_column,
                            action = "add") {
  UseMethod("annotate_symbol", input.df)
}

#' @export
annotate_symbol.default <-  function(input.df,
                                     ensembl_transcript_column,
                                     action = "add")
{
  print("This function cannot be applied to this object")
}

#' @export
annotate_symbol.tbl_df = annotate_symbol.ttBulk <-
  function(input.df,
           ensembl_transcript_column,
           action = "add")
  {
    # Make col names
    ensembl_transcript_column = enquo(ensembl_transcript_column)


    if (action == "add")
      add_symbol_from_ensembl(input.df,!!ensembl_transcript_column)

    else if (action == "get")
      get_symbol_from_ensembl(input.df,!!ensembl_transcript_column)

    else
      stop(
        "action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
      )

  }

# ANNOTATE differential transcription

#' Add differential transcription information to a tibble using edgeR. At the moment only one covariate is accepted
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr set_colnames
#'
#'
#'
#' @param input.df A tibble
#' @param formula a formula with no response variable, referring only to numeric variables
#' @param sample_column A character name of the sample column
#' @param transcript_column A character name of the gene/transcript name column
#' @param counts_column A character name of the read count column
#' @param significance_threshold A real between 0 and 1.
#' @param action A character string. Whether to join the new information to the input tibble (add), or just get the non-redundant tibble with the new information (get).
#'
#' @return A tibble with differential_transcript_abundance results
#'
#' @export
annotate_differential_transcription <- function(input.df,
                                                formula,
                                                sample_column = NULL,
                                                transcript_column = NULL,
                                                counts_column = NULL,
                                                significance_threshold = 0.05,
                                                action = "add",
                                                ...) {
  UseMethod("annotate_differential_transcription", input.df)
}

#' @export
annotate_differential_transcription.default <-  function(input.df,
                                                         formula,
                                                         sample_column = NULL,
                                                         transcript_column = NULL,
                                                         counts_column = NULL,
                                                         significance_threshold = 0.05,
                                                         action = "add",
                                                         ...)
{
  print("This function cannot be applied to this object")
}

#' @export
annotate_differential_transcription.tbl_df = annotate_differential_transcription.ttBulk <-
  function(input.df,
           formula,
           sample_column = NULL,
           transcript_column = NULL,
           counts_column = NULL,
           significance_threshold = 0.05,
           action = "add")
  {
    # Make col names
    sample_column = enquo(sample_column)
    transcript_column = enquo(transcript_column)
    counts_column = enquo(counts_column)

    if (action == "add")
      add_differential_transcript_abundance_bulk(
        input.df,
        formula,
        sample_column = !!sample_column,
        transcript_column = !!transcript_column,
        counts_column = !!counts_column,
        significance_threshold = significance_threshold
      )
    else if (action == "get")
      get_differential_transcript_abundance_bulk(
        input.df,
        formula,
        sample_column = !!sample_column,
        transcript_column = !!transcript_column,
        counts_column = !!counts_column,
        significance_threshold = significance_threshold
      )
    else
      stop(
        "action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
      )
  }
