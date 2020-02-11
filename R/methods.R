setOldClass("spec_tbl_df")
setOldClass("ttBulk")




#' Creates a `tt` object from a `tbl``
#'
#' \lifecycle{maturing}
#'
#' @description ttBulk() creates a `tt` object from a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#' @import lifecycle
#'
#' @name ttBulk
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param .abundance_normalised The name of the transcript/gene normalised abundance column
#'
#' @details This function created a ttBulk object and is useful if you want
#' to avoid to specify .sample, .transcript and .abundance arguments all the times.
#' The ttBulk object have an attribute called parameters where these three
#' arguments are stored as metadata. They can be extracted as attr(<object>, "parameters").
#'
#' @return A `ttBulk` object
#'
#'
#' @examples
#'
#'
#'
#'
#' my_tt =  ttBulk(ttBulk::counts_mini, sample, transcript, count)
#'
#'
#' @docType methods
#' @rdname ttBulk-methods
#' @export
#'
setGeneric("ttBulk", function(.data,
                              .sample,
                              .transcript,
                              .abundance,
															.abundance_normalised = NULL) standardGeneric("ttBulk"))

# Set internal
.ttBulk = function(.data,
                   .sample,
                   .transcript,
                   .abundance,
									 .abundance_normalised = NULL) {

  # Make col names
  .sample = enquo(.sample)
  .transcript = enquo(.transcript)
  .abundance = enquo(.abundance)
  .abundance_normalised = enquo(.abundance_normalised)

  # Validate data frame
  validation(.data,!!.sample,!!.transcript,!!.abundance, skip_dupli_check = TRUE)

  create_tt_from_tibble_bulk(.data,!!.sample,!!.transcript,!!.abundance, !!.abundance_normalised)
}
#' ttBulk
#' @inheritParams ttBulk
#' @return A `ttBulk` object
#'
setMethod("ttBulk", "spec_tbl_df", .ttBulk)

#' ttBulk
#'
#'
#' @inheritParams ttBulk
#' @return A `ttBulk` object
#'
setMethod("ttBulk", "tbl_df", .ttBulk)

.ttBulk_se = function(.data,
											 .sample,
											 .transcript,
											 .abundance,
											.abundance_normalised = NULL){

	# Make col names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	.abundance_normalised = enquo(.abundance_normalised)

	# Set normalised col names
	norm_col =
		assays(.data)[1] %>% names %>% paste("normalised") %>%
		ifelse_pipe(
			(.) %in% names(assays(.data)),
			~ as.symbol(.x),
			~ NULL
		)

	# Do donversion
	assays(.data) %>%
		as.list() %>%
		map2(
			assays(.data) %>%  names,
			~ .x %>%
				as_tibble(rownames = "feature") %>%
				gather(sample, !!.y, -feature)
		) %>%

		# Join the assays
		purrr::reduce(dplyr::left_join, by = c("sample", "feature")) %>%

		# Attach annotation
		left_join(
			rowData(.data) %>% as.data.frame() %>% as_tibble(rownames = "feature"),
			by = "feature"
		) %>%
		left_join(
			colData(.data) %>% as_tibble(rownames="sample"),
			by = "sample"
		) %>%
		mutate_if(is.character, as.factor) %>%
		ttBulk(
			sample,
			feature,
			!!as.symbol(assays(.data)[1] %>%  names),

			# Normalised counts if any
			!!norm_col
		)

}

#' ttBulk
#'
#' @importFrom tibble as_tibble
#' @importFrom purrr reduce
#' @import dplyr
#' @import tidyr
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors DataFrame
#' @importFrom purrr map2
#'
#'
#' @inheritParams ttBulk
#' @return A `ttBulk` object
#'
setMethod("ttBulk", "SummarizedExperiment", .ttBulk_se)

#' ttBulk
#' @inheritParams ttBulk
#' @return A `ttBulk` object
#'
setMethod("ttBulk", "RangedSummarizedExperiment", .ttBulk_se)


#' Creates a `tt` object from a list of file names of BAM/SAM
#'
#' \lifecycle{maturing}
#'
#' @description ttBulk_SAM_BAM() creates a `tt` object from a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name ttBulk_SAM_BAM
#'
#' @param file_names A character vector
#' @param genome A character string
#' @param ... Further parameters passed to the function Rsubread::featureCounts
#'
#' @details This function is based on FeatureCounts package. This function created a ttBulk object and is useful if you want
#' to avoid to specify .sample, .transcript and .abundance arguments all the times.
#' The ttBulk object have an attribute called parameters where these three
#' arguments are stored as metadata. They can be extracted as attr(<object>, "parameters").
#'
#' @return A `ttBulk` object
#'
#'
#'
#'
#'
#' @docType methods
#' @rdname ttBulk_SAM_BAM-methods
#' @export
#'
setGeneric("ttBulk_SAM_BAM", function(file_names, genome = "hg38", ...) standardGeneric("ttBulk_SAM_BAM"))

#' ttBulk_SAM_BAM
#' @inheritParams ttBulk_SAM_BAM
#' @return A `ttBulk` object
#'
setMethod("ttBulk_SAM_BAM", c(file_names = "character", genome = "character"), 	function(file_names, genome = "hg38", ...) create_tt_from_bam_sam_bulk(file_names = file_names, genome = genome, ...))

#' Normalise the counts of transcripts/genes
#'
#' \lifecycle{maturing}
#'
#' @description scale_abundance() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and normalises the data for the library size (e.g., with TMM algorithm, Robinson and Oshlack doi.org/10.1186/gb-2010-11-3-r25).
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#' @importFrom stats median
#'
#' @name scale_abundance
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#'
#' @param cpm_threshold A real positive number. It is the threshold of count per million that is used to filter transcripts/genes out from the normalisation procedure. The normalisation inference is then applied back to all unfiltered data.
#' @param prop A real positive number between 0 and 1. It is the threshold of proportion of samples for each transcripts/genes that have to be characterised by a cmp bigger than the threshold to be included for normalisation procedure.
#' @param method A character string. The methods used by the function. The method is passed to the function `calcNormFactors` from edgeR package.
#' @param reference_selection_function A fucntion that is used to selecting the reference sample for normalisation. It could be max (default), which choose the sample with maximum library size; or median, which chooses the sample with median library size.
#' @param action A character string between "add" (default) and "get". "add" joins the new information to the input tbl (default), "get" return a non-redundant tbl with the just new information.
#'
#' @details normalises the data for the library size
#' (e.g., with TMM algorithm, Robinson and Oshlack doi.org/10.1186/gb-2010-11-3-r25).
#' Lowly transcribed transcripts/genes (defined with cpm_threshold and prop parameters)
#' are filtered out from the normalisation procedure.
#' The normalisation inference is then applied back to all unfiltered data.
#'
#' @return A tbl object with additional columns with normalised data as `<NAME OF COUNT COLUMN> normalised`
#'
#'
#' @examples
#'
#'
#'  scale_abundance(ttBulk::counts_mini,  sample, transcript, `count`)
#'
#'
#'
#' @docType methods
#' @rdname scale_abundance-methods
#' @export

setGeneric("scale_abundance", function(.data,
                              .sample = NULL,
                              .transcript = NULL,
                              .abundance = NULL,
                              cpm_threshold = 0.5,
                              prop = 3 / 4,
                              method = "TMM",
                              reference_selection_function = median,
                              action = "add") standardGeneric("scale_abundance"))

# Set internal
.scale_abundance = 	function(.data,
                             .sample = NULL,
                             .transcript = NULL,
                             .abundance = NULL,
                             cpm_threshold = 0.5,
                             prop = 3 / 4,
                             method = "TMM",
                             reference_selection_function = median,
                             action = "add")
{
  # Make col names
  .sample = enquo(.sample)
  .transcript = enquo(.transcript)
  .abundance = enquo(.abundance)

  # Validate data frame
  validation(.data,!!.sample,!!.transcript,!!.abundance)

  if (action == "add")
    add_normalised_counts_bulk(
      .data,!!.sample,!!.transcript,!!.abundance,
      cpm_threshold = cpm_threshold,
      prop = prop,
      method = method,
      reference_selection_function = reference_selection_function
    )
  else if (action == "get")
    get_normalised_counts_bulk(
      .data,!!.sample,!!.transcript,!!.abundance,
      cpm_threshold = cpm_threshold,
      prop = prop,
      method = method,
      reference_selection_function = reference_selection_function
    )
  else
    stop(
      "action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
    )
}

#' scale_abundance
#' @inheritParams scale_abundance
#' @return A tbl object with additional columns with normalised data as `<NAME OF COUNT COLUMN> normalised`
#'
setMethod("scale_abundance", "spec_tbl_df", .scale_abundance)

#' scale_abundance
#' @inheritParams scale_abundance
#' @return A tbl object with additional columns with normalised data as `<NAME OF COUNT COLUMN> normalised`
#'
setMethod("scale_abundance", "tbl_df", .scale_abundance)

#' scale_abundance
#' @inheritParams scale_abundance
#' @return A tbl object with additional columns with normalised data as `<NAME OF COUNT COLUMN> normalised`
#'
setMethod("scale_abundance", "ttBulk", .scale_abundance)

.scale_abundance_se = function(.data,
															 .sample = NULL,
															 .transcript = NULL,
															 .abundance = NULL,
															 cpm_threshold = 0.5,
															 prop = 3 / 4,
															 method = "TMM",
															 reference_selection_function = median,
															 action = "add"){

	# Make col names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	.data %>%

		# Convert to ttBulk
		ttBulk() %>%

		# Apply scale method
		scale_abundance(!!.sample,
										!!.transcript,
										!!.abundance,
										cpm_threshold = cpm_threshold,
										prop = prop,
										method = method,
										reference_selection_function = reference_selection_function,
										action = action) %>%

		# Convert to SummaizedExperiment
		ttBulk_to_SummarizedExperiment()

}

#' scale_abundance
#' @inheritParams scale_abundance
#' @return A `ttBulk` object
#'
setMethod("scale_abundance", "SummarizedExperiment", .scale_abundance_se)

#' scale_abundance
#' @inheritParams scale_abundance
#' @return A `ttBulk` object
#'
setMethod("scale_abundance", "RangedSummarizedExperiment", .scale_abundance_se)



#' Get clusters of elements (e.g., samples or transcripts)
#'
#' \lifecycle{maturing}
#'
#' @description cluster_elements() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and identify clusters in the data.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name cluster_elements
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .element The name of the element column (normally samples).
#' @param .feature The name of the feature column (normally transcripts/genes)
#' @param .abundance The name of the column including the numerical value the clustering is based on (normally transcript abundance)
#'
#' @param method A character string. The cluster algorithm to use, ay the moment k-means is the only algorithm included.
#' @param of_samples A boolean. In case the input is a ttBulk object, it indicates Whether the element column will be sample or transcript column
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#' @param ... Further parameters passed to the function kmeans

#' @details identifies clusters in the data, normally of samples.
#' This function returns a tibble with additional columns for the cluster annotation.
#' At the moment only k-means clustering is supported, the plan is to introduce more clustering methods.
#'
#' @return A tbl object with additional columns with cluster labels
#'
#'
#' @examples
#'
#'
#'     cluster_elements(ttBulk::counts_mini, sample, transcript, count,	centers = 2, method="kmeans")
#'
#' @docType methods
#' @rdname cluster_elements-methods
#' @export
#'
setGeneric("cluster_elements", function(.data,
                                        .element = NULL,
                                        .feature = NULL,
                                        .abundance = NULL,
                                        method,
                                        of_samples = TRUE,
                                        log_transform = TRUE,
                                        action = "add",
                                        ...) standardGeneric("cluster_elements"))

# Set internal
.cluster_elements = 		function(.data,
                               .element = NULL,
                               .feature = NULL,
                               .abundance = NULL,
                               method ,
                               of_samples = TRUE,
                               log_transform = TRUE,
                               action = "add",
                               ...)
{
  # Make col names
  .abundance = enquo(.abundance)
  .element = enquo(.element)
  .feature = enquo(.feature)

  # Validate data frame
  validation(.data,!!.element,!!.feature,!!.abundance)

  if (method == "kmeans") {
    if (action == "add")
      add_clusters_kmeans_bulk(
        .data,
        .abundance = !!.abundance,
        .element = !!.element,
        .feature = !!.feature,
        of_samples = of_samples,
        log_transform = log_transform,
        ...
      )
    else if (action == "get")
      get_clusters_kmeans_bulk(
        .data,
        .abundance = !!.abundance,
        .element = !!.element,
        .feature = !!.feature,
        of_samples = of_samples,
        log_transform = log_transform,
        ...
      )
    else
      stop(
        "action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
      )
  }
  else if (method == "SNN") {
    if (action == "add")
      add_clusters_SNN_bulk(
        .data,
        .abundance = !!.abundance,
        .element = !!.element,
        .feature = !!.feature,
        of_samples = of_samples,
        log_transform = log_transform,
        ...
      )
    else if (action == "get")
      get_clusters_SNN_bulk(
        .data,
        .abundance = !!.abundance,
        .element = !!.element,
        .feature = !!.feature,
        of_samples = of_samples,
        log_transform = log_transform,
        ...
      )
    else
      stop(
        "action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
      )
  }
  else
    stop("the only supported methods are \"kmeans\" or \"SNN\" ")

}

#' cluster_elements
#' @inheritParams cluster_elements
#' @return A tbl object with additional columns with cluster labels
#'
setMethod("cluster_elements", "spec_tbl_df", .cluster_elements)

#' cluster_elements
#' @inheritParams cluster_elements
#' @return A tbl object with additional columns with cluster labels
#'
setMethod("cluster_elements", "tbl_df", .cluster_elements)

#' cluster_elements
#' @inheritParams cluster_elements
#' @return A tbl object with additional columns with cluster labels
#'
setMethod("cluster_elements", "ttBulk", .cluster_elements)

.cluster_elements_se = function(.data,
																.element = NULL,
																.feature = NULL,
																.abundance = NULL,
																method ,
																of_samples = TRUE,
																log_transform = TRUE,
																action = "add",
																...){

	# Make col names
	.abundance = enquo(.abundance)
	.element = enquo(.element)
	.feature = enquo(.feature)

	.data %>%

		# Convert to ttBulk
		ttBulk() %>%

		# Apply scale method
		cluster_elements(.element = !!.element ,
										 .feature = !!.feature ,
										 .abundance = !!.abundance,
										 method = method ,
										 of_samples = of_samples,
										 log_transform = log_transform,
										 action = action,
										 ...) %>%

		# Convert to SummaizedExperiment
		ttBulk_to_SummarizedExperiment()

}

#' cluster_elements
#' @inheritParams cluster_elements
#' @return A `ttBulk` object
#'
setMethod("cluster_elements", "SummarizedExperiment", .cluster_elements_se)

#' cluster_elements
#' @inheritParams cluster_elements
#' @return A `ttBulk` object
#'
setMethod("cluster_elements", "RangedSummarizedExperiment", .cluster_elements_se)


#' Dimension reduction of the transcript abundance data
#'
#' \lifecycle{maturing}
#'
#' @description reduce_dimensions() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and calculates the reduced dimensional space of the transcript abundance.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name reduce_dimensions
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .element The name of the element column (normally samples).
#' @param .feature The name of the feature column (normally transcripts/genes)
#' @param .abundance The name of the column including the numerical value the clustering is based on (normally transcript abundance)
#'
#' @param method A character string. The dimension reduction algorithm to use (PCA, MDS, tSNE).
#' @param top An integer. How many top genes to select for dimensionality reduction
#' @param of_samples A boolean. In case the input is a ttBulk object, it indicates Whether the element column will be sample or transcript column
#' @param .dims A list of integer vectors corresponding to principal components of interest (e.g., list(1:2, 3:4, 5:6))

#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#' @param scale A boolean for method="PCA", this will be passed to the `prcomp` function. It is not included in the ... argument because although the default for `prcomp` if FALSE, it is advisable to set it as TRUE.
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#' @param ... Further parameters passed to the function prcomp if you choose method="PCA" or Rtsne if you choose method="tSNE"
#'
#' @details This function reduces the dimensions of the transcript abundances.
#' It can use multi-dimensional scaling (MDS) of principal component analysis (PCA).
#'
#' @return A tbl object with additional columns for the reduced dimensions
#'
#'
#' @examples
#'
#'
#'
#' counts.MDS =  reduce_dimensions(ttBulk::counts_mini, sample, transcript, count, method="MDS", .dims = 3)
#'
#'
#' counts.PCA =  reduce_dimensions(ttBulk::counts_mini, sample, transcript, count, method="PCA", .dims = 3)
#'
#'
#'
#' @docType methods
#' @rdname reduce_dimensions-methods
#' @export
#'
#'
setGeneric("reduce_dimensions", function(.data,
															.element = NULL,
															.feature = NULL,
															.abundance = NULL,
															method,
															.dims = 2,

															top = 500,
															of_samples = TRUE,
															log_transform = TRUE,
															scale = TRUE,
															action = "add",
															...) standardGeneric("reduce_dimensions"))

# Set internal
.reduce_dimensions = 		function(.data,
                                .element = NULL,
                                .feature = NULL,
                                .abundance = NULL,
                                method,
                                .dims = 2,

                                top = 500,
                                of_samples = TRUE,
                                log_transform = TRUE,
                                scale = TRUE,
                                action = "add",
                                ...)
{
  # Make col names
  .abundance = enquo(.abundance)
  .element = enquo(.element)
  .feature = enquo(.feature)

  # Validate data frame
  validation(.data,!!.element,!!.feature,!!.abundance)

  if (method == "MDS") {
    if (action == "add")
      add_reduced_dimensions_MDS_bulk(
        .data,
        .abundance = !!.abundance,
        .dims = .dims,
        .element = !!.element,
        .feature = !!.feature,
        top = top,
        of_samples = of_samples,
        log_transform = log_transform,
        ...
      )
    else if (action == "get")
      get_reduced_dimensions_MDS_bulk(
        .data,
        .abundance = !!.abundance,
        .dims = .dims,
        .element = !!.element,
        .feature = !!.feature,
        top = top,
        of_samples = of_samples,
        log_transform = log_transform,
        ...
      )
    else
      stop(
        "action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
      )
  }
  else if (method == "PCA") {
    if (action == "add")
      add_reduced_dimensions_PCA_bulk(
        .data,
        .abundance = !!.abundance,
        .dims = .dims,
        .element = !!.element,
        .feature = !!.feature,
        top = top,
        of_samples = of_samples,
        log_transform = log_transform,
        scale = scale,
        ...
      )
    else if (action == "get")
      get_reduced_dimensions_PCA_bulk(
        .data,
        .abundance = !!.abundance,
        .dims = .dims,
        .element = !!.element,
        .feature = !!.feature,
        top = top,
        of_samples = of_samples,
        log_transform = log_transform,
        scale = scale,
        ...
      )
    else
      stop(
        "action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
      )

  }
  else if (method == "tSNE") {
    if (action == "add")
      add_reduced_dimensions_TSNE_bulk(
        .data,
        .abundance = !!.abundance,
        .dims = .dims,
        .element = !!.element,
        .feature = !!.feature,
        top = top,
        of_samples = of_samples,
        log_transform = log_transform,
        scale = scale,
        ...
      )
    else if (action == "get")
      get_reduced_dimensions_TSNE_bulk(
        .data,
        .abundance = !!.abundance,
        .dims = .dims,
        .element = !!.element,
        .feature = !!.feature,
        top = top,
        of_samples = of_samples,
        log_transform = log_transform,
        scale = scale,
        ...
      )
    else
      stop(
        "action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
      )

  }
  else
    stop("method must be either \"MDS\" or \"PCA\"")

}

#' reduce_dimensions
#' @inheritParams reduce_dimensions
#' @return A tbl object with additional columns for the reduced dimensions
setMethod("reduce_dimensions", "spec_tbl_df", .reduce_dimensions)

#' reduce_dimensions
#' @inheritParams reduce_dimensions
#' @return A tbl object with additional columns for the reduced dimensions
setMethod("reduce_dimensions", "tbl_df", .reduce_dimensions)

#' reduce_dimensions
#' @inheritParams reduce_dimensions
#' @return A tbl object with additional columns for the reduced dimensions
setMethod("reduce_dimensions", "ttBulk", .reduce_dimensions)

.reduce_dimensions_se = function(.data,
																 .element = NULL,
																 .feature = NULL,
																 .abundance = NULL,
																 method,
																 .dims = 2,

																 top = 500,
																 of_samples = TRUE,
																 log_transform = TRUE,
																 scale = TRUE,
																 action = "add",
																 ...){

	# Make col names
	.abundance = enquo(.abundance)
	.element = enquo(.element)
	.feature = enquo(.feature)

	.data %>%

		# Convert to ttBulk
		ttBulk() %>%

		# Apply scale method
		reduce_dimensions(.element = !!.element,
											.feature  = !!.feature,
											.abundance  = !!.abundance,
											method = method,
											.dims = .dims,

											top = top,
											of_samples = of_samples,
											log_transform = log_transform,
											scale = scale,
											action = action,
											...) %>%

		# Convert to SummaizedExperiment
		ttBulk_to_SummarizedExperiment()

}

#' reduce_dimensions
#' @inheritParams reduce_dimensions
#' @return A `ttBulk` object
#'
setMethod("reduce_dimensions", "SummarizedExperiment", .reduce_dimensions_se)

#' reduce_dimensions
#' @inheritParams reduce_dimensions
#' @return A `ttBulk` object
#'
setMethod("reduce_dimensions", "RangedSummarizedExperiment", .reduce_dimensions_se)

#' Rotate two dimensions (e.g., principal components) of an arbitrary angle
#'
#' \lifecycle{maturing}
#'
#' @description rotate_dimensions() takes as imput a `tbl` formatted as | <DIMENSION 1> | <DIMENSION 2> | <...> | and calculates the rotated dimensional space of the transcript abundance.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name rotate_dimensions
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .element The name of the element column (normally samples).
#'
#' @param dimension_1_column A character string. The column of the dimension 1
#' @param dimension_2_column  A character string. The column of the dimension 2
#' @param rotation_degrees A real number between 0 and 360
#' @param of_samples A boolean. In case the input is a ttBulk object, it indicates Whether the element column will be sample or transcript column
#' @param dimension_1_column_rotated A character string. The column of the rotated dimension 1 (optional)
#' @param dimension_2_column_rotated A character string. The column of the rotated dimension 2 (optional)
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#'
#' @details This function to rotate two dimensions such as the reduced dimensions.
#'
#' @return A tbl object with additional columns for the reduced dimensions. additional columns for the rotated dimensions. The rotated dimensions will be added to the original data set as `<NAME OF DIMENSION> rotated <ANGLE>` by default, or as specified in the input arguments.
#'
#'
#' @examples
#'
#' counts.MDS =  reduce_dimensions(ttBulk::counts_mini, sample, transcript, count, method="MDS", .dims = 3)
#'
#' counts.MDS.rotated =  rotate_dimensions(counts.MDS, `Dim1`, `Dim2`, rotation_degrees = 45, .element = sample)
#'
#'
#' @docType methods
#' @rdname rotate_dimensions-methods
#' @export
#'
setGeneric("rotate_dimensions", function(.data,
															dimension_1_column,
															dimension_2_column,
															rotation_degrees,
															.element = NULL,
															of_samples = TRUE,
															dimension_1_column_rotated = NULL,
															dimension_2_column_rotated = NULL,
															action = "add") standardGeneric("rotate_dimensions"))

# Set internal
.rotate_dimensions = 		function(.data,
					 dimension_1_column,
					 dimension_2_column,
					 rotation_degrees,
					 .element = NULL,
					 of_samples = TRUE,
					 dimension_1_column_rotated = NULL,
					 dimension_2_column_rotated = NULL,
					 action =
					 	"add")
	{
		# Make col names
		.element = enquo(.element)
		dimension_1_column = enquo(dimension_1_column)
		dimension_2_column = enquo(dimension_2_column)
		dimension_1_column_rotated = enquo(dimension_1_column_rotated)
		dimension_2_column_rotated = enquo(dimension_2_column_rotated)


		if (action == "add")
			add_rotated_dimensions(
				.data,
				dimension_1_column = !!dimension_1_column,
				dimension_2_column = !!dimension_2_column,
				rotation_degrees = rotation_degrees,
				.element = !!.element,
				of_samples = of_samples,
				dimension_1_column_rotated = !!dimension_1_column_rotated,
				dimension_2_column_rotated = !!dimension_2_column_rotated
			)
		else if (action == "get")
			get_rotated_dimensions(
				.data,
				dimension_1_column = !!dimension_1_column,
				dimension_2_column = !!dimension_2_column,
				rotation_degrees = rotation_degrees,
				.element = !!.element,
				of_samples = of_samples,
				dimension_1_column_rotated = !!dimension_1_column_rotated,
				dimension_2_column_rotated = !!dimension_2_column_rotated
			)
		else
			stop(
				"action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
			)
	}

#' rotate_dimensions
#' @inheritParams rotate_dimensions
#' @return A tbl object with additional columns for the reduced dimensions. additional columns for the rotated dimensions. The rotated dimensions will be added to the original data set as `<NAME OF DIMENSION> rotated <ANGLE>` by default, or as specified in the input arguments.
setMethod("rotate_dimensions", "spec_tbl_df", .rotate_dimensions)

#' rotate_dimensions
#' @inheritParams rotate_dimensions
#' @return A tbl object with additional columns for the reduced dimensions. additional columns for the rotated dimensions. The rotated dimensions will be added to the original data set as `<NAME OF DIMENSION> rotated <ANGLE>` by default, or as specified in the input arguments.
setMethod("rotate_dimensions", "tbl_df", .rotate_dimensions)

#' rotate_dimensions
#' @inheritParams rotate_dimensions
#' @return A tbl object with additional columns for the reduced dimensions. additional columns for the rotated dimensions. The rotated dimensions will be added to the original data set as `<NAME OF DIMENSION> rotated <ANGLE>` by default, or as specified in the input arguments.
setMethod("rotate_dimensions", "ttBulk", .rotate_dimensions)


.rotate_dimensions_se = function(.data,
																 dimension_1_column,
																 dimension_2_column,
																 rotation_degrees,
																 .element = NULL,
																 of_samples = TRUE,
																 dimension_1_column_rotated = NULL,
																 dimension_2_column_rotated = NULL,
																 action =
																 	"add"){

	# Make col names
	.element = enquo(.element)
	dimension_1_column = enquo(dimension_1_column)
	dimension_2_column = enquo(dimension_2_column)
	dimension_1_column_rotated = enquo(dimension_1_column_rotated)
	dimension_2_column_rotated = enquo(dimension_2_column_rotated)

	.data %>%

		# Convert to ttBulk
		ttBulk() %>%

		# Apply scale method
		rotate_dimensions(dimension_1_column = !!dimension_1_column,
											dimension_2_column = !!dimension_2_column,
											rotation_degrees = rotation_degrees,
											.element = !!.element,
											of_samples = of_samples,
											dimension_1_column_rotated = !!dimension_1_column_rotated,
											dimension_2_column_rotated = !!dimension_2_column_rotated,
											action = action) %>%

		# Convert to SummaizedExperiment
		ttBulk_to_SummarizedExperiment()

}

#' rotate_dimensions
#' @inheritParams rotate_dimensions
#' @return A `ttBulk` object
#'
setMethod("rotate_dimensions", "SummarizedExperiment", .rotate_dimensions_se)

#' rotate_dimensions
#' @inheritParams rotate_dimensions
#' @return A `ttBulk` object
#'
setMethod("rotate_dimensions", "RangedSummarizedExperiment", .rotate_dimensions_se)


#' Drop redundant elements (e.g., samples) for which feature (e.g., transcript/gene) aboundances are correlated
#'
#' \lifecycle{maturing}
#'
#' @description remove_redundancy() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | for correlation method or | <DIMENSION 1> | <DIMENSION 2> | <...> | for reduced_dimensions method, and returns a `tbl` with dropped elements (e.g., samples).
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name remove_redundancy
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .element The name of the element column (normally samples).
#' @param .feature The name of the feature column (normally transcripts/genes)
#' @param .abundance The name of the column including the numerical value the clustering is based on (normally transcript abundance)
#'
#' @param method A character string. The cluster algorithm to use, ay the moment k-means is the only algorithm included.
#' @param of_samples A boolean. In case the input is a ttBulk object, it indicates Whether the element column will be sample or transcript column
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#' @param correlation_threshold A real number between 0 and 1. For correlation based calculation.
#' @param top An integer. How many top genes to select for correlation based method
#' @param Dim_a_column A character string. For reduced_dimension based calculation. The column of one principal component
#' @param Dim_b_column A character string. For reduced_dimension based calculation. The column of another principal component
#'
#'
#' @details This function removes redundant elements from the original data set (e.g., samples or transcripts). For example, if we want to define cell-type specific signatures with low sample redundancy. This function returns a tibble with dropped recundant elements (e.g., samples). Two redundancy estimation approaches are supported: (i) removal of highly correlated clusters of elements (keeping a representative) with method="correlation"; (ii) removal of most proximal element pairs in a reduced dimensional space.
#'
#' @return A tbl object with with dropped recundant elements (e.g., samples).
#'
#' @examples
#'
#'
#'
#'    remove_redundancy(
#'     ttBulk::counts_mini,
#' 	   .element = sample,
#' 	   .feature = transcript,
#' 	   	.abundance =  count,
#' 	   	method = "correlation"
#' 	   	)
#'
#' counts.MDS =  reduce_dimensions(ttBulk::counts_mini, sample, transcript, count, method="MDS", .dims = 3)
#'
#' remove_redundancy(
#' 	counts.MDS,
#' 	Dim_a_column = `Dim1`,
#' 	Dim_b_column = `Dim2`,
#' 	.element = sample,
#'   method = "reduced_dimensions"
#' )
#'
#' @docType methods
#' @rdname remove_redundancy-methods
#' @export
#'
#'
setGeneric("remove_redundancy", function(.data,
													 .element = NULL,
													 .feature = NULL,
													 .abundance = NULL,
													 method,

													 of_samples = TRUE,



													 correlation_threshold = 0.9,
													 top = Inf,
													 log_transform = FALSE,

													 Dim_a_column,
													 Dim_b_column) standardGeneric("remove_redundancy"))

# Set internal
.remove_redundancy = 	 function(.data,
																													 .element = NULL,
																													 .feature = NULL,
																													 .abundance = NULL,
																													 method,

																													 of_samples = TRUE,



																													 correlation_threshold = 0.9,
																													 top = Inf,
																													 log_transform = FALSE,

																													 Dim_a_column = NULL,
																													 Dim_b_column = NULL)
{
	# Make col names
	.abundance = enquo(.abundance)
	.element = enquo(.element)
	.feature = enquo(.feature)

	Dim_a_column = enquo(Dim_a_column)
	Dim_b_column = enquo(Dim_b_column)

	if (method == "correlation"){

		# Validate data frame
		validation(.data,!!.element,!!.feature,!!.abundance)

		remove_redundancy_elements_through_correlation(
			.data,
			.abundance = !!.abundance,
			.element = !!.element,
			.feature = !!.feature,
			correlation_threshold = correlation_threshold,
			top = top,
			of_samples = of_samples,
			log_transform = log_transform
		)
	}
	else if (method == "reduced_dimensions"){

		# Validate data frame
		# MISSING because feature not needed. I should build a custom funtion.

		remove_redundancy_elements_though_reduced_dimensions(
			.data,
			Dim_a_column = !!Dim_a_column,
			Dim_b_column = !!Dim_b_column,
			.element = !!.element,
			of_samples = of_samples
		)
	}
	else
		stop(
			"method must be either \"correlation\" for dropping correlated elements or \"reduced_dimension\" to drop the closest pair according to two dimensions (e.g., PCA)"
		)

}

#' remove_redundancy
#' @inheritParams remove_redundancy
#' @return A tbl object with with dropped recundant elements (e.g., samples).
setMethod("remove_redundancy", "spec_tbl_df", .remove_redundancy)

#' remove_redundancy
#' @inheritParams remove_redundancy
#' @return A tbl object with with dropped recundant elements (e.g., samples).
setMethod("remove_redundancy", "tbl_df", .remove_redundancy)

#' remove_redundancy
#' @inheritParams remove_redundancy
#' @return A tbl object with with dropped recundant elements (e.g., samples).
setMethod("remove_redundancy", "ttBulk", .remove_redundancy)

.remove_redundancy_se = function(.data,
																 .element = NULL,
																 .feature = NULL,
																 .abundance = NULL,
																 method,

																 of_samples = TRUE,



																 correlation_threshold = 0.9,
																 top = Inf,
																 log_transform = FALSE,

																 Dim_a_column = NULL,
																 Dim_b_column = NULL){

	# Make col names
	.abundance = enquo(.abundance)
	.element = enquo(.element)
	.feature = enquo(.feature)

	Dim_a_column = enquo(Dim_a_column)
	Dim_b_column = enquo(Dim_b_column)

	.data %>%

		# Convert to ttBulk
		ttBulk() %>%

		# Apply scale method
		remove_redundancy( .element = !!.element,
											 .feature = !!.feature,
											 .abundance = !!.abundance,
											 method = method,

											 of_samples = of_samples,



											 correlation_threshold = correlation_threshold,
											 top = top,
											 log_transform = log_transform,

											 Dim_a_column = !!Dim_a_column,
											 Dim_b_column = !!Dim_b_column) %>%

		# Convert to SummaizedExperiment
		ttBulk_to_SummarizedExperiment()

}

#' remove_redundancy
#' @inheritParams remove_redundancy
#' @return A `ttBulk` object
#'
setMethod("remove_redundancy", "SummarizedExperiment", .remove_redundancy_se)

#' remove_redundancy
#' @inheritParams remove_redundancy
#' @return A `ttBulk` object
#'
setMethod("remove_redundancy", "RangedSummarizedExperiment", .remove_redundancy_se)



#' Adjust transcript abundance for unwanted variation
#'
#' \lifecycle{maturing}
#'
#' @description adjust_abundance() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and returns a `tbl` with an edditional adjusted abundance column..
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name adjust_abundance
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .formula A formula with no response variable, representing the desired linear model where the first covariate is the factor of interest and the second covariate is the unwanted variation (of the kind ~ factor_of_intrest + batch)
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#'
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#' @param ... Further parameters passed to the function sva::ComBat
#'
#' @details This function adjusts the abundance for (known) unwanted variation. At the moment just an unwanted covariated is allowed at a time.
#'
#' @return A `tbl` with additional columns for the adjusted counts as `<COUNT COLUMN> adjusted`
#'
#'
#'
#'
#' @examples
#'
#'
#'
#' cm = ttBulk::counts_mini
#' cm$batch = 0
#' cm$batch[cm$sample %in% c("SRR1740035", "SRR1740043")] = 1
#'
#' res =
#' 	adjust_abundance(
#' 		cm,
#'		~ condition + batch,
#'		.sample = sample,
#'		.transcript = transcript,
#'		.abundance = count
#'	)
#'
#'
#' @docType methods
#' @rdname adjust_abundance-methods
#' @export
#'
#'
setGeneric("adjust_abundance", function(.data,
														 .formula,
														 .sample = NULL,
														 .transcript = NULL,
														 .abundance = NULL,
														 log_transform = TRUE,
														 action = "add",
														 ...) standardGeneric("adjust_abundance"))

# Set internal
.adjust_abundance = 	function(.data,
					 .formula,
					 .sample = NULL,
					 .transcript = NULL,
					 .abundance = NULL,
					 log_transform = TRUE,
					 action = "add",
					 ...)
	{
		# Make col names
		.sample = enquo(.sample)
		.transcript = enquo(.transcript)
		.abundance = enquo(.abundance)

		# Validate data frame
		validation(.data,!!.sample,!!.transcript,!!.abundance)

		if (action == "add")
			add_adjusted_counts_for_unwanted_variation_bulk(
				.data,
				.formula,
				.sample = !!.sample,
				.transcript = !!.transcript,
				.abundance = !!.abundance,
				log_transform = log_transform,
				...
			)
		else if (action == "get")
			get_adjusted_counts_for_unwanted_variation_bulk(
				.data,
				.formula,
				.sample = !!.sample,
				.transcript = !!.transcript,
				.abundance = !!.abundance,
				log_transform = log_transform,
				...
			)
		else
			stop(
				"action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
			)
	}

#' adjust_abundance
#' @inheritParams adjust_abundance
#' @return A `tbl` with additional columns for the adjusted counts as `<COUNT COLUMN> adjusted`
setMethod("adjust_abundance", "spec_tbl_df", .adjust_abundance)

#' adjust_abundance
#' @inheritParams adjust_abundance
#' @return A `tbl` with additional columns for the adjusted counts as `<COUNT COLUMN> adjusted`
setMethod("adjust_abundance", "tbl_df", .adjust_abundance)

#' adjust_abundance
#' @inheritParams adjust_abundance
#' @return A `tbl` with additional columns for the adjusted counts as `<COUNT COLUMN> adjusted`
setMethod("adjust_abundance", "ttBulk", .adjust_abundance)

.adjust_abundance_se = function(.data,
																.formula,
																.sample = NULL,
																.transcript = NULL,
																.abundance = NULL,
																log_transform = TRUE,
																action = "add",
																...){

	# Make col names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	.data %>%

		# Convert to ttBulk
		ttBulk() %>%

		# Apply scale method
		adjust_abundance( .formula = .formula,
											.sample = !!.sample,
											.transcript = !!.transcript,
											.abundance = !!.abundance,
											log_transform = log_transform,
											action = action,
											...) %>%

		# Convert to SummaizedExperiment
		ttBulk_to_SummarizedExperiment()

}

#' adjust_abundance
#' @inheritParams adjust_abundance
#' @return A `ttBulk` object
#'
setMethod("adjust_abundance", "SummarizedExperiment", .adjust_abundance_se)

#' adjust_abundance
#' @inheritParams adjust_abundance
#' @return A `ttBulk` object
#'
setMethod("adjust_abundance", "RangedSummarizedExperiment", .adjust_abundance_se)



#' Aggregates multiple counts from the same samples (e.g., from isoforms), concatenates other character columns, and averages other numeric columns
#'
#' \lifecycle{maturing}
#'
#' @description aggregate_duplicates() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and returns a `tbl` with aggregated transcripts that were duplicated.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name aggregate_duplicates
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#'
#' @param aggregation_function A function for counts aggregation (e.g., sum,  median, or mean)
#' @param keep_integer A boolean. Whether to force the aggregated counts to integer
#'
#' @details This function aggregates duplicated transcripts (e.g., isoforms, ensembl).
#' For example, we often have to convert ensembl symbols to gene/transcript symbol,
#'  but in doing so we have to deal with duplicates. `aggregate_duplicates` takes a tibble
#'  and column names (as symbols; for `sample`, `transcript` and `count`) as arguments and
#'  returns a tibble with aggregate transcript with the same name. All the rest of the column
#'  are appended, and factors and boolean are appended as characters.
#'
#' @return A `tbl` object with aggregated transcript abundance and annotation
#'
#'
#'
#'
#' @examples
#'
#'     aggregate_duplicates(
#'     ttBulk::counts_mini,
#'     sample,
#'     transcript,
#'     `count`,
#'     aggregation_function = sum
#'     )
#'
#'
#' @docType methods
#' @rdname aggregate_duplicates-methods
#' @export
#'
#'
setGeneric("aggregate_duplicates", function(.data,

																 .sample = NULL,
																 .transcript = NULL,
																 .abundance = NULL,
																 aggregation_function = sum,
																 keep_integer = TRUE) standardGeneric("aggregate_duplicates"))

# Set internal
.aggregate_duplicates = 	function(.data,

					 .sample = NULL,
					 .transcript = NULL,
					 .abundance = NULL,
					 aggregation_function = sum,
					 keep_integer = TRUE)  {
		# Make col names
		.sample = enquo(.sample)
		.transcript = enquo(.transcript)
		.abundance = enquo(.abundance)

		# Validate data frame
		validation(.data,!!.sample,!!.transcript,!!.abundance, skip_dupli_check = TRUE)

		aggregate_duplicated_transcripts_bulk(
			.data,

			.sample = !!.sample,
			.transcript = !!.transcript,
			.abundance = !!.abundance,
			aggregation_function = aggregation_function,
			keep_integer = TRUE
		)
	}

#' aggregate_duplicates
#' @inheritParams aggregate_duplicates
#' @return A `tbl` object with aggregated transcript abundance and annotation
setMethod("aggregate_duplicates", "spec_tbl_df", .aggregate_duplicates)

#' aggregate_duplicates
#' @inheritParams aggregate_duplicates
#' @return A `tbl` object with aggregated transcript abundance and annotation
setMethod("aggregate_duplicates", "tbl_df", .aggregate_duplicates)

#' aggregate_duplicates
#' @inheritParams aggregate_duplicates
#' @return A `tbl` object with aggregated transcript abundance and annotation
setMethod("aggregate_duplicates", "ttBulk", .aggregate_duplicates)

.aggregate_duplicates_se = function(.data,

																		.sample = NULL,
																		.transcript = NULL,
																		.abundance = NULL,
																		aggregation_function = sum,
																		keep_integer = TRUE){

	# Make col names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	.data %>%

		# Convert to ttBulk
		ttBulk() %>%

		# Apply scale method
		aggregate_duplicates(
													.sample = !!.sample,
													.transcript = !!.transcript,
													.abundance = !!.abundance,
													aggregation_function = aggregation_function,
													keep_integer = keep_integer) %>%

		# Convert to SummaizedExperiment
		ttBulk_to_SummarizedExperiment()

}

#' aggregate_duplicates
#' @inheritParams aggregate_duplicates
#' @return A `ttBulk` object
#'
setMethod("aggregate_duplicates", "SummarizedExperiment", .aggregate_duplicates_se)

#' aggregate_duplicates
#' @inheritParams aggregate_duplicates
#' @return A `ttBulk` object
#'
setMethod("aggregate_duplicates", "RangedSummarizedExperiment", .aggregate_duplicates_se)



#' Get cell type proportions from samples
#'
#' \lifecycle{maturing}
#'
#' @description deconvolve_cellularity() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and returns a `tbl` with the estimated cell type abundance for each sample
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name deconvolve_cellularity
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param reference A data frame. The transcript/cell_type data frame of integer transcript abundance
#' @param method A character string. The method to be used. At the moment Cibersort (default) and llsr (linear least squares regression) are available.
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#' @param ... Further parameters passed to the function Cibersort
#'
#' @details This function infers the cell type composition of our samples (with the algorithm Cibersort; Newman et al., 10.1038/nmeth.3337).
#'
#' @return A `tbl` object including additional columns for each cell type estimated
#'
#'
#'
#'
#' @examples
#'
#'
#' deconvolve_cellularity(ttBulk::counts, sample, transcript, `count`, cores = 2)
#'
#'
#' @docType methods
#' @rdname deconvolve_cellularity-methods
#' @export
#'
setGeneric("deconvolve_cellularity", function(.data,
															 .sample = NULL,
															 .transcript = NULL,
															 .abundance = NULL,
															 reference = X_cibersort,
															 method = "cibersort",
															 action = "add",
															 ...) standardGeneric("deconvolve_cellularity"))

# Set internal
.deconvolve_cellularity = 		function(.data,
					 .sample = NULL,
					 .transcript = NULL,
					 .abundance = NULL,
					 reference = X_cibersort,
					 method = "cibersort",
					 action = "add",
					 ...)  {
		# Make col names
		.sample = enquo(.sample)
		.transcript = enquo(.transcript)
		.abundance = enquo(.abundance)

		# Validate data frame
		validation(.data,!!.sample,!!.transcript,!!.abundance)

		if (action == "add")
			add_cell_type_proportions(
				.data,
				.sample = !!.sample,
				.transcript = !!.transcript,
				.abundance = !!.abundance,
				reference = reference,
				method = method,
				...
			)
		else if (action == "get")
			get_cell_type_proportions(
				.data,
				.sample = !!.sample,
				.transcript = !!.transcript,
				.abundance = !!.abundance,
				reference = reference,
				method = method,
				...
			)
		else
			stop(
				"action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
			)
	}

#' deconvolve_cellularity
#' @inheritParams deconvolve_cellularity
#' @return A `tbl` object including additional columns for each cell type estimated
setMethod("deconvolve_cellularity", "spec_tbl_df", .deconvolve_cellularity)

#' deconvolve_cellularity
#' @inheritParams deconvolve_cellularity
#' @return A `tbl` object including additional columns for each cell type estimated
setMethod("deconvolve_cellularity", "tbl_df", .deconvolve_cellularity)

#' deconvolve_cellularity
#' @inheritParams deconvolve_cellularity
#' @return A `tbl` object including additional columns for each cell type estimated
setMethod("deconvolve_cellularity", "ttBulk", .deconvolve_cellularity)



.deconvolve_cellularity_se = function(.data,
																			.sample = NULL,
																			.transcript = NULL,
																			.abundance = NULL,
																			reference = X_cibersort,
																			method = "cibersort",
																			action = "add",
																			...){

	# Make col names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	.data %>%

		# Convert to ttBulk
		ttBulk() %>%

		# Apply scale method
		deconvolve_cellularity(	 .sample = !!.sample,
													 .transcript = !!.transcript,
													 .abundance = !!.abundance,
													 reference = reference,
													 method = method,
													 action = action,
													 ...) %>%

		# Convert to SummaizedExperiment
		ttBulk_to_SummarizedExperiment()

}

#' deconvolve_cellularity
#' @inheritParams deconvolve_cellularity
#' @return A `ttBulk` object
#'
setMethod("deconvolve_cellularity", "SummarizedExperiment", .deconvolve_cellularity_se)

#' deconvolve_cellularity
#' @inheritParams deconvolve_cellularity
#' @return A `ttBulk` object
#'
setMethod("deconvolve_cellularity", "RangedSummarizedExperiment", .deconvolve_cellularity_se)




#' Add transcript symbol column from ensembl id
#'
#' \lifecycle{maturing}
#'
#' @description annotate_symbol() takes as imput a `tbl` formatted as | <SAMPLE> | <ENSEMBL_ID> | <COUNT> | <...> | and returns a `tbl` with the the additional transcript symbol column
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name annotate_symbol
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <ENSEMBL_ID> | <COUNT> | <...> |
#' @param .ensembl A character string. The column that is represents ensembl gene id
#'
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#'
#' @details This is useful since different resources use ensembl IDs while others use gene symbol IDs.
#'
#' @return A `tbl` object including additional columns for transcript symbol
#'
#'
#'
#'
#' @examples
#'
#'
#' 	annotate_symbol(ttBulk::counts_ensembl, ens)
#'
#'
#' @docType methods
#' @rdname annotate_symbol-methods
#' @export
#'
#'
setGeneric("annotate_symbol", function(.data,
														.ensembl,
														action = "add") standardGeneric("annotate_symbol"))

# Set internal
.annotate_symbol = 		function(.data,
					 .ensembl,
					 action = "add")
	{
		# Make col names
		.ensembl = enquo(.ensembl)


		if (action == "add")
			add_symbol_from_ensembl(.data, !!.ensembl)

		else if (action == "get")
			get_symbol_from_ensembl(.data, !!.ensembl)

		else
			stop(
				"action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
			)

	}

#' annotate_symbol
#' @inheritParams annotate_symbol
#' @return A `tbl` object including additional columns for transcript symbol
setMethod("annotate_symbol", "spec_tbl_df", .annotate_symbol)

#' annotate_symbol
#' @inheritParams annotate_symbol
#' @return A `tbl` object including additional columns for transcript symbol
setMethod("annotate_symbol", "tbl_df", .annotate_symbol)

#' annotate_symbol
#' @inheritParams annotate_symbol
#' @return A `tbl` object including additional columns for transcript symbol
setMethod("annotate_symbol", "ttBulk", .annotate_symbol)


#' Add differential transcription information to a tbl using edgeR.
#'
#' \lifecycle{maturing}
#'
#' @description test_differential_abundance() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and returns a `tbl` with additional columns for the statistics from the hypothesis test.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name test_differential_abundance
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .formula A formula with no response variable, representing the desired linear model
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param .coef An integer. See edgeR specifications
#' @param .contrasts A character vector. See edgeR makeContrasts specification for the parameter `contrasts`
#' @param significance_threshold A real between 0 and 1 (usually 0.05).
#' @param fill_missing_values A boolean. Whether to fill missing sample/transcript values with the median of the transcript. This is rarely needed.
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#'
#' @details At the moment this function uses edgeR only, but other inference algorithms will be added in the near future.
#'
#' @return A `tbl` with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
#'
#'
#'
#'
#' @examples
#'
#'
#' 	test_differential_abundance(
#' 	 ttBulk::counts_mini,
#' 	    ~ condition,
#' 	    sample,
#' 	    transcript,
#' 	    `count`
#' 	)
#'
#'
#' @docType methods
#' @rdname test_differential_abundance-methods
#' @export
#'
setGeneric("test_differential_abundance", function(.data,
																						.formula,
																						.sample = NULL,
																						.transcript = NULL,
																						.abundance = NULL,
																				.coef = 2,
																				.contrasts = NULL,
																						significance_threshold = 0.05,
																				fill_missing_values = FALSE,

																						action = "add") standardGeneric("test_differential_abundance"))

# Set internal
.test_differential_abundance = 		function(.data,
					 .formula,
					 .sample = NULL,
					 .transcript = NULL,
					 .abundance = NULL,
					 .coef = 2,
					 .contrasts = NULL,
					 significance_threshold = 0.05,
					 fill_missing_values = FALSE,
					 action = "add")
	{
		# Make col names
		.sample = enquo(.sample)
		.transcript = enquo(.transcript)
		.abundance = enquo(.abundance)

		# Validate data frame
		validation(.data,!!.sample,!!.transcript,!!.abundance)

		if (action == "add")
			add_differential_transcript_abundance_bulk(
				.data,
				.formula,
				.sample = !!.sample,
				.transcript = !!.transcript,
				.abundance = !!.abundance,
				.coef = .coef,
				.contrasts = .contrasts,
				significance_threshold = significance_threshold,
				fill_missing_values = fill_missing_values
			)
		else if (action == "get")
			get_differential_transcript_abundance_bulk(
				.data,
				.formula,
				.sample = !!.sample,
				.transcript = !!.transcript,
				.abundance = !!.abundance,
				.coef = .coef,
				.contrasts = .contrasts,
				significance_threshold = significance_threshold,
				fill_missing_values = fill_missing_values
			)
		else
			stop(
				"action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
			)
	}

#' test_differential_abundance
#' @inheritParams test_differential_abundance
#' @return A `tbl` with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
setMethod("test_differential_abundance", "spec_tbl_df", .test_differential_abundance)

#' test_differential_abundance
#' @inheritParams test_differential_abundance
#' @return A `tbl` with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
setMethod("test_differential_abundance", "tbl_df", .test_differential_abundance)

#' test_differential_abundance
#' @inheritParams test_differential_abundance
#' @return A `tbl` with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
setMethod("test_differential_abundance", "ttBulk", .test_differential_abundance)



.test_differential_abundance_se = function(.data,
																					 .formula,
																					 .sample = NULL,
																					 .transcript = NULL,
																					 .abundance = NULL,
																					 .coef = 2,
																					 .contrasts = NULL,
																					 significance_threshold = 0.05,
																					 fill_missing_values = FALSE,
																					 action = "add")
{
	# Make col names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	.data %>%

		# Convert to ttBulk
		ttBulk() %>%

		# Apply scale method
		test_differential_abundance(.formula,
																.sample = !!.sample,
																.transcript = !!.transcript,
																.abundance = !!.abundance,
																.coef = .coef,
																.contrasts = .contrasts,
																significance_threshold = significance_threshold,
																fill_missing_values = fill_missing_values,
																action = action) %>%

		# Convert to SummaizedExperiment
		ttBulk_to_SummarizedExperiment()

}

#' test_differential_abundance
#' @inheritParams test_differential_abundance
#' @return A `ttBulk` object
#'
setMethod("test_differential_abundance", "SummarizedExperiment", .test_differential_abundance_se)

#' test_differential_abundance
#' @inheritParams test_differential_abundance
#' @return A `ttBulk` object
#'
setMethod("test_differential_abundance", "RangedSummarizedExperiment", .test_differential_abundance_se)



#' Filter variable transcripts
#'
#' \lifecycle{maturing}
#'
#' @description filter_variable() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and returns a `tbl` with additional columns for the statistics from the hypothesis test.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name filter_variable
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param top Integer. Number of top transcript to consider
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#'
#' @details At the moment this function uses edgeR only, but other inference algorithms will be added in the near future.
#'
#' @return A `tbl` with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
#'
#'
#'
#'
#' @examples
#'
#'
#'
#' 	filter_variable(
#' 	ttBulk::counts_mini,
#' 	    sample,
#' 	    transcript,
#' 	    `count`,
#' 	    top = 500
#' 	)
#'
#'
#' @docType methods
#' @rdname filter_variable-methods
#' @export
#'
setGeneric("filter_variable", function(.data,
																				.sample = NULL,
																				.transcript = NULL,
																				.abundance = NULL,
																				top = 500,
														log_transform = TRUE) standardGeneric("filter_variable"))

# Set internal
.filter_variable = 		function(.data,
					 .sample = NULL,
					 .transcript = NULL,
					 .abundance = NULL,
					 top = 500,
					 log_transform = TRUE)
	{
		# Make col names
		.sample = enquo(.sample)
		.transcript = enquo(.transcript)
		.abundance = enquo(.abundance)

		# Validate data frame
		validation(.data,!!.sample,!!.transcript,!!.abundance)

		filter_variable_transcripts(
				.data,
				.sample = !!.sample,
				.transcript = !!.transcript,
				.abundance = !!.abundance,
				top = top,
				log_transform = log_transform
			)
	}

#' filter_variable
#' @inheritParams filter_variable
#' @return A `tbl` with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
setMethod("filter_variable", "spec_tbl_df", .filter_variable)

#' filter_variable
#' @inheritParams filter_variable
#' @return A `tbl` with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
setMethod("filter_variable", "tbl_df", .filter_variable)

#' filter_variable
#' @inheritParams filter_variable
#' @return A `tbl` with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
setMethod("filter_variable", "ttBulk", .filter_variable)

.filter_variable_se = function(.data,
															 .sample = NULL,
															 .transcript = NULL,
															 .abundance = NULL,
															 top = 500,
															 log_transform = TRUE)
{
	# Make col names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	.data %>%

		# Convert to ttBulk
		ttBulk() %>%

		# Apply scale method
		filter_variable(
										.sample = !!.sample,
										.transcript = !!.transcript,
										.abundance = !!.abundance,
										top = top,
										log_transform = log_transform) %>%

		# Convert to SummaizedExperiment
		ttBulk_to_SummarizedExperiment()

}

#' filter_variable
#' @inheritParams filter_variable
#' @return A `ttBulk` object
#'
setMethod("filter_variable", "SummarizedExperiment", .filter_variable_se)

#' filter_variable
#' @inheritParams filter_variable
#' @return A `ttBulk` object
#'
setMethod("filter_variable", "RangedSummarizedExperiment", .filter_variable_se)



#' Filter abundant transcripts
#'
#' \lifecycle{maturing}
#'
#' @description filter_abundant() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and returns a `tbl` with additional columns for the statistics from the hypothesis test.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#' @importFrom dplyr filter
#'
#' @name filter_abundant
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param cpm_threshold A real positive number
#' @param prop A number between 0 and 1
#'
#' @details At the moment this function uses edgeR only, but other inference algorithms will be added in the near future.
#'
#' @return A `tbl` with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
#'
#'
#'
#'
#' @examples
#'
#'
#'
#' 	filter_abundant(
#' 	ttBulk::counts_mini,
#' 	    sample,
#' 	    transcript,
#' 	    `count`
#' 	)
#'
#'
#' @docType methods
#' @rdname filter_abundant-methods
#' @export
#'
setGeneric("filter_abundant", function(.data,
														.sample = NULL,
														.transcript = NULL,
														.abundance = NULL,
														cpm_threshold = 0.5,
														prop = 3 / 4) standardGeneric("filter_abundant"))

# Set internal
.filter_abundant = 		function(.data,
					 .sample = NULL,
					 .transcript = NULL,
					 .abundance = NULL,
					 cpm_threshold = 0.5,
					 prop = 3 / 4)
	{

	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
	.sample = col_names$.sample
	.transcript = col_names$.transcript
	.abundance = col_names$.abundance

		# Validate data frame
		validation(.data,!!.sample,!!.transcript,!!.abundance)

		gene_to_exclude =
			add_normalised_counts_bulk.get_low_expressed(.data,
																								 .sample = !!.sample,
																								 .transcript = !!.transcript,
																								 .abundance = !!.abundance,
																														 cpm_threshold = cpm_threshold,
																														 prop = prop)

		.data %>%

			# Filter
			dplyr::filter(!!.transcript %in% gene_to_exclude %>% `!`) %>%

			# Attach attributes
			add_attr(.data %>% attr("parameters"), "parameters")
	}

#' filter_abundant
#' @inheritParams filter_abundant
#' @return A `tbl` with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
setMethod("filter_abundant", "spec_tbl_df", .filter_abundant)

#' filter_abundant
#' @inheritParams filter_abundant
#' @return A `tbl` with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
setMethod("filter_abundant", "tbl_df", .filter_abundant)

#' filter_abundant
#' @inheritParams filter_abundant
#' @return A `tbl` with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
setMethod("filter_abundant", "ttBulk", .filter_abundant)

.filter_abundant_se = function(.data,
															 .sample = NULL,
															 .transcript = NULL,
															 .abundance = NULL,
															 cpm_threshold = 0.5,
															 prop = 3 / 4)
{
	# Make col names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	.data %>%

		# Convert to ttBulk
		ttBulk() %>%

		# Apply scale method
		filter_abundant(
			.sample = !!.sample,
			.transcript = !!.transcript,
			.abundance = !!.abundance,
			cpm_threshold = cpm_threshold,
			prop = prop) %>%

		# Convert to SummaizedExperiment
		ttBulk_to_SummarizedExperiment()

}

#' filter_abundant
#' @inheritParams filter_abundant
#' @return A `ttBulk` object
#'
setMethod("filter_abundant", "SummarizedExperiment", .filter_abundant_se)

#' filter_abundant
#' @inheritParams filter_abundant
#' @return A `ttBulk` object
#'
setMethod("filter_abundant", "RangedSummarizedExperiment", .filter_abundant_se)



#' Analise gene enrichment with EGSEA
#'
#' \lifecycle{maturing}
#'
#' @description analise_gene_enrichment() takes as imput a `tbl` formatted as | <SAMPLE> | <ENSEMBL_ID> | <COUNT> | <...> | and returns a `tbl` with the the additional transcript symbol column
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name analise_gene_enrichment
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .formula A formula with no response variable, representing the desired linear model
#' @param .sample The name of the sample column
#' @param .entrez The ENTREZ doce of the transcripts/genes
#' @param .abundance The name of the transcript/gene abundance column
#' @param .contrasts = NULL,
#' @param species A character. For example, human or mouse
#' @param cores An integer. The number of cores available
#'
#'
#' @details This wrapper execute gene enrichment analyses of the dataset
#'
#' @return A `tbl` object
#'
#'
#'
#'
#' @examples
#' \donttest{
#'
#' df_entrez = symbol_to_entrez(ttBulk::counts_mini, .transcript = transcript, .sample = sample)
#' df_entrez = aggregate_duplicates(df_entrez, aggregation_function = sum, .sample = sample, .transcript = entrez, .abundance = count)
#'
#' library("EGSEA")
#'
#' 	analise_gene_enrichment(
#'			df_entrez,
#'			~ condition,
#'			.sample = sample,
#'			.entrez = entrez,
#'			.abundance = count,
#'			species="human"
#'		)
#'
#'}
#'
#' @docType methods
#' @rdname analise_gene_enrichment-methods
#' @export
#'
#'
setGeneric("analise_gene_enrichment", function(.data,
																		.formula,
																		.sample = NULL,
																		.entrez,
																		.abundance = NULL,
																		.contrasts = NULL,
																		species,
																		cores = 10) standardGeneric("analise_gene_enrichment"))

# Set internal
.analise_gene_enrichment = 		function(.data,
					 .formula,
					 .sample = NULL,
					 .entrez,
					 .abundance = NULL,
					 .contrasts = NULL,
					 species,
					 cores = 10)	{
		# Make col names
		.sample = enquo(.sample)
		.abundance = enquo(.abundance)
		.entrez = enquo(.entrez)

		# Validate data frame
		validation(.data,!!.sample,!!.entrez)

		analyse_gene_enrichment_bulk_EGSEA(.data,
																			 .formula,
																			 .sample = !!.sample,
																			 .entrez = !!.entrez,
																			 .abundance = !!.abundance,
																			 .contrasts = .contrasts,
																			 species = species,
																			 cores = cores)



	}

#' analise_gene_enrichment
#' @inheritParams analise_gene_enrichment
#' @return A `tbl` object
setMethod("analise_gene_enrichment", "spec_tbl_df", .analise_gene_enrichment)

#' analise_gene_enrichment
#' @inheritParams analise_gene_enrichment
#' @return A `tbl` object
setMethod("analise_gene_enrichment", "tbl_df", .analise_gene_enrichment)

#' analise_gene_enrichment
#' @inheritParams analise_gene_enrichment
#' @return A `tbl` object
setMethod("analise_gene_enrichment", "ttBulk", .analise_gene_enrichment)
