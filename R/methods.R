

# CREATE


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
  create_tt_from_tibble_bulk(input.df, sample_column, transcript_column,  counts_column)

}


# NORMALISE


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

#' @export
find_clusters <- function(input.df,
                          value_column,
                          number_of_clusters,
                          elements_column = NULL,
                          feature_column = NULL,
                          of_samples = T,
                          log_transform = T,
                          action = "add") {
  UseMethod("find_clusters", input.df)
}

#' @export
find_clusters.default <-  function(input.df,
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
find_clusters.tbl_df = find_clusters.ttBulk <-  function(input.df,
                                                         value_column,
                                                         number_of_clusters,
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

# REDUCE DIMENTIONS

#' @export
reduce_dimensions <- function(input.df,
                              value_column,
                              method,
                              components_list = list(c(1, 2)),
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
                                       components_list = list(c(1, 2)),
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
           components_list = list(c(1, 2)),
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
          components_list = components_list,
          elements_column = !!elements_column,
          feature_column = !!feature_column,
          log_transform = log_transform
        )
      else if (action == "get")
        get_reduced_dimensions_MDS_bulk(
          input.df,
          value_column = !!value_column,
          components_list = components_list,
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
          components_list = components_list,
          elements_column = !!elements_column,
          feature_column = !!feature_column,
          log_transform = log_transform
        )
      else if (action == "get")
        get_reduced_dimensions_PCA_bulk(
          input.df,
          value_column = !!value_column,
          components_list = components_list,
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


  if (action == "correlation")
    drop_redundant_elements_through_correlation(
      input.df,
      value_column = !!value_column,
      elements_column = !!elements_column,
      feature_column = !!feature_column,
      correlation_threshold = correlation_threshold,
      of_samples = of_samples,
      log_transform = log_transform
    )
  else if (action == "reduced_dimensions")
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
      "action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
    )

}


# ADJUST


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
           keep_integer = T)
  {
    # Make col names
    sample_column = enquo(sample_column)
    transcript_column = enquo(transcript_column)
    counts_column = enquo(counts_column)

    if (action == "add")
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
           action = "add")
  {
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
