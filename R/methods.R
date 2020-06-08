# setOldClass("spec_tbl_df")
setOldClass("tidybulk")

#' Creates a `tt` object from a `tbl``
#'
#' \lifecycle{maturing}
#'
#' @description tidybulk() creates a `tt` object from a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#' @import readr
#'
#' @name tidybulk
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param .abundance_scaled The name of the transcript/gene scaled abundance column
#'
#' @details This function created a tidybulk object and is useful if you want
#' to avoid to specify .sample, .transcript and .abundance arguments all the times.
#' The tidybulk object have an attribute called internals where these three
#' arguments are stored as metadata. They can be extracted as attr(<object>, "internals").
#'
#' @return A `tidybulk` object
#'
#'
#' @examples
#'
#'
#'
#'
#' my_tt =  tidybulk(tidybulk::counts_mini, sample, transcript, count)
#'
#'
#' @docType methods
#' @rdname tidybulk-methods
#' @export
#'
setGeneric("tidybulk", function(.data,
																.sample,
																.transcript,
																.abundance,
																.abundance_scaled = NULL)
	standardGeneric("tidybulk"))

# Set internal
.tidybulk = function(.data,
										 .sample,
										 .transcript,
										 .abundance,
										 .abundance_scaled = NULL) {
	# Make col names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	.abundance_scaled = enquo(.abundance_scaled)

	# Validate data frame
	validation(.data,
						 !!.sample,
						 !!.transcript,
						 !!.abundance,
						 skip_dupli_check = TRUE)

	create_tt_from_tibble_bulk(.data,
														 !!.sample,
														 !!.transcript,
														 !!.abundance,
														 !!.abundance_scaled)
}
#' tidybulk
#' @inheritParams tidybulk
#' 
#' @docType methods
#' @rdname tidybulk-methods
#' 
#' @return A `tidybulk` object
#'
setMethod("tidybulk", "spec_tbl_df", .tidybulk)

#' tidybulk
#'
#' @importFrom purrr map2
#'
#' @inheritParams tidybulk
#' 
#' @docType methods
#' @rdname tidybulk-methods
#' 
#' @return A `tidybulk` object
#'
setMethod("tidybulk", "tbl_df", .tidybulk)

.tidybulk_se = function(.data,
											.sample,
											.transcript,
											.abundance,
											.abundance_scaled = NULL) {
	# Check if package is installed, otherwise install
	if (find.package("SummarizedExperiment", quiet = TRUE) %>% length %>% equals(0)) {
		message("Installing SummarizedExperiment")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("SummarizedExperiment", ask = FALSE)
	}

	# Make col names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	.abundance_scaled = enquo(.abundance_scaled)

	# Set scaled col names
	norm_col =
		SummarizedExperiment::assays(.data)[1] %>% names %>% paste0(scaled_string) %>%
		ifelse_pipe((.) %in% names(SummarizedExperiment::assays(.data)),
								~ as.symbol(.x),
								~ NULL)

	# Do conversion
	SummarizedExperiment::assays(.data) %>%
		as.list() %>%
		map2(
			SummarizedExperiment::assays(.data) %>%  names,
			~ .x %>%
				as_tibble(rownames = "feature") %>%
				gather(sample,!!.y,-feature)
		) %>%

		# Join the assays
		purrr::reduce(dplyr::left_join, by = c("sample", "feature")) %>%

		# Attach annotation
		left_join(
			SummarizedExperiment::rowData(.data) %>% as.data.frame() %>% as_tibble(rownames = "feature"),
			by = "feature"
		) %>%
		left_join(SummarizedExperiment::colData(.data) %>% as_tibble(rownames =
																																 	"sample"),
							by = "sample") %>%
		mutate_if(is.character, as.factor) %>%
		tidybulk(
			sample,
			feature,
			!!as.symbol(SummarizedExperiment::assays(.data)[1] %>%  names	),
			!!norm_col # scaled counts if any
		)

}

#' tidybulk
#'
#' @importFrom tibble as_tibble
#' @importFrom purrr reduce
#' @import dplyr
#' @import tidyr
#'
#'
#' @inheritParams tidybulk
#' 
#' @docType methods
#' @rdname tidybulk-methods
#' 
#' @return A `tidybulk` object
#'
setMethod("tidybulk", "SummarizedExperiment", .tidybulk_se)

#' tidybulk
#' @inheritParams tidybulk
#' 
#' @docType methods
#' @rdname tidybulk-methods
#' 
#' @return A `tidybulk` object
#'
setMethod("tidybulk", "RangedSummarizedExperiment", .tidybulk_se)


#' Creates a `tt` object from a list of file names of BAM/SAM
#'
#' \lifecycle{maturing}
#'
#' @description tidybulk_SAM_BAM() creates a `tt` object from a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name tidybulk_SAM_BAM
#'
#' @param file_names A character vector
#' @param genome A character string
#' @param ... Further parameters passed to the function Rsubread::featureCounts
#'
#' @details This function is based on FeatureCounts package. This function created a tidybulk object and is useful if you want
#' to avoid to specify .sample, .transcript and .abundance arguments all the times.
#' The tidybulk object have an attribute called internals where these three
#' arguments are stored as metadata. They can be extracted as attr(<object>, "internals").
#'
#' @return A `tidybulk` object
#'
#'
#'
#'
#'
#' @docType methods
#' @rdname tidybulk_SAM_BAM-methods
#' @export
#'
setGeneric("tidybulk_SAM_BAM", function(file_names, genome = "hg38", ...)
	standardGeneric("tidybulk_SAM_BAM"))

#' tidybulk_SAM_BAM
#' @inheritParams tidybulk_SAM_BAM-methods
#' 
#' @docType methods
#' @rdname tidybulk_SAM_BAM-methods
#' 
#' @return A `tidybulk` object
#'
setMethod("tidybulk_SAM_BAM", c(file_names = "character", genome = "character"), 	function(file_names, genome = "hg38", ...)
	create_tt_from_bam_sam_bulk(file_names = file_names, genome = genome, ...))

#' Scale the counts of transcripts/genes
#'
#' \lifecycle{maturing}
#'
#' @description scale_abundance() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and Scales transcript abundance compansating for sequencing depth (e.g., with TMM algorithm, Robinson and Oshlack doi.org/10.1186/gb-2010-11-3-r25).
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
#' @param factor_of_interest The name of the column of the factor of interest. This is used for identifying lowly abundant transcript, to be ignored for calculating scaling fators.
#' @param minimum_counts A real positive number. It is the threshold of count per million that is used to filter transcripts/genes out from the scaling procedure. The scaling inference is then applied back to all unfiltered data.
#' @param minimum_proportion A real positive number between 0 and 1. It is the threshold of proportion of samples for each transcripts/genes that have to be characterised by a cmp bigger than the threshold to be included for scaling procedure.
#' @param method A character string. The scaling method passed to the backend function (i.e., edgeR::calcNormFactors; "TMM","TMMwsp","RLE","upperquartile")
#' @param reference_selection_function A fucntion that is used to selecting the reference sample for scaling. It could be max (default), which choose the sample with maximum library size; or median, which chooses the sample with median library size.
#' @param action A character string between "add" (default) and "only". "add" joins the new information to the input tbl (default), "only" return a non-redundant tbl with the just new information.
#'
#' @details Scales transcript abundance compansating for sequencing depth
#' (e.g., with TMM algorithm, Robinson and Oshlack doi.org/10.1186/gb-2010-11-3-r25).
#' Lowly transcribed transcripts/genes (defined with minimum_counts and minimum_proportion parameters)
#' are filtered out from the scaling procedure.
#' The scaling inference is then applied back to all unfiltered data.
#'
#' @return A tbl object with additional columns with scaled data as `<NAME OF COUNT COLUMN>_scaled`
#'
#'
#' @examples
#'
#'
#'  scale_abundance(tidybulk::counts_mini,  sample, transcript, `count`)
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
																			 factor_of_interest = NULL,
																			 minimum_counts = 10,
																			 minimum_proportion = 0.7,
																			 method = "TMM",
																			 reference_selection_function = median,
																			 action = "add")
	standardGeneric("scale_abundance"))

# Set internal
.scale_abundance = 	function(.data,
														 .sample = NULL,
														 .transcript = NULL,
														 .abundance = NULL,
														 factor_of_interest = NULL,
														 minimum_counts = 10,
														 minimum_proportion = 0.7,
														 method = "TMM",
														 reference_selection_function = median,
														 action = "add")
{
	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
	.sample = col_names$.sample
	.transcript = col_names$.transcript
	.abundance = col_names$.abundance

	factor_of_interest = enquo(factor_of_interest)

	# Validate data frame
	validation(.data, !!.sample, !!.transcript, !!.abundance)

	.data_norm =
		.data %>%
		get_scaled_counts_bulk(
			.sample = !!.sample,
			.transcript = !!.transcript,
			.abundance = !!.abundance,
			factor_of_interest = !!factor_of_interest,
			minimum_counts = minimum_counts,
			minimum_proportion = minimum_proportion,
			method = method,
			reference_selection_function = reference_selection_function
		) %>%
		arrange(!!.sample,!!.transcript)


	if (action == "add"){

		.data %>%
			arrange(!!.sample,!!.transcript) %>%

			# Add scaled data set
			bind_cols(.data_norm %>%
									select(-one_of(quo_name(.sample)), -one_of(quo_name(.transcript))))		%>%

			# Attach attributes
			reattach_internals(.data_norm)

	}
	else if (action == "get"){

		.data %>%

			# Selecting the right columns
			select(
				!!.sample,
				get_x_y_annotation_columns(.data, !!.sample,!!.transcript, !!.abundance, NULL)$horizontal_cols
			) %>%
			distinct() %>%
			mutate_if(is.character, as.factor) %>%

			# Join result
			left_join(.data_norm, by=quo_name(.sample)) %>%

			# Attach attributes
			reattach_internals(.data_norm)
	}
	else if (action == "only") .data_norm
	else
		stop(
			"tidybulk says: action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
		)
}

#' scale_abundance
#' @inheritParams scale_abundance
#' 
#' @docType methods
#' @rdname scale_abundance-methods
#' 
#' @return A tbl object with additional columns with scaled data as `<NAME OF COUNT COLUMN>_scaled`
#'
setMethod("scale_abundance", "spec_tbl_df", .scale_abundance)

#' scale_abundance
#' @inheritParams scale_abundance
#' 
#' @docType methods
#' @rdname scale_abundance-methods
#' 
#' @return A tbl object with additional columns with scaled data as `<NAME OF COUNT COLUMN>_scaled`
#'
setMethod("scale_abundance", "tbl_df", .scale_abundance)

#' scale_abundance
#' @inheritParams scale_abundance
#' 
#' @docType methods
#' @rdname scale_abundance-methods
#' 
#' @return A tbl object with additional columns with scaled data as `<NAME OF COUNT COLUMN>_scaled`
#'
setMethod("scale_abundance", "tidybulk", .scale_abundance)

.scale_abundance_se = function(.data,
															 .sample = NULL,
															 .transcript = NULL,
															 .abundance = NULL,
															 factor_of_interest = NULL,
															 minimum_counts = 10,
															 minimum_proportion = 0.7,
															 method = "TMM",
															 reference_selection_function = median,
															 action = "add") {
	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	factor_of_interest = enquo(factor_of_interest)

	.data %>%

		# Convert to tidybulk
		tidybulk() %>%

		# Apply scale method
		scale_abundance(
			!!.sample,
			!!.transcript,
			!!.abundance,
			factor_of_interest = !!factor_of_interest,
			minimum_counts = minimum_counts,
			minimum_proportion = minimum_proportion,
			method = method,
			reference_selection_function = reference_selection_function,
			action = action
		) %>%

		# Convert to SummaizedExperiment
		tidybulk_to_SummarizedExperiment()

}

#' scale_abundance
#' @inheritParams scale_abundance
#' 
#' @docType methods
#' @rdname scale_abundance-methods
#' 
#' @return A `SummarizedExperiment` object
#'
setMethod("scale_abundance",
					"SummarizedExperiment",
					.scale_abundance_se)

#' scale_abundance
#' @inheritParams scale_abundance
#' 
#' @docType methods
#' @rdname scale_abundance-methods
#' 
#' @return A `SummarizedExperiment` object
#'
setMethod("scale_abundance",
					"RangedSummarizedExperiment",
					.scale_abundance_se)



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
#' @param of_samples A boolean. In case the input is a tidybulk object, it indicates Whether the element column will be sample or transcript column
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#' @param ... Further parameters passed to the function kmeans
#'
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
#'     cluster_elements(tidybulk::counts_mini, sample, transcript, count,	centers = 2, method="kmeans")
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
																				...)
	standardGeneric("cluster_elements"))

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
	# Get column names
	.element = enquo(.element)
	.feature = enquo(.feature)
	col_names = get_elements_features(.data, .element, .feature, of_samples)
	.element = col_names$.element
	.feature = col_names$.feature

	# Get scaled abundance if present, otherwise get abundance
	.abundance = enquo(.abundance)
	col_names = get_abundance_norm_if_exists(.data, .abundance)
	.abundance = col_names$.abundance

	# Validate data frame
	validation(.data, !!.element, !!.feature, !!.abundance)

	if (method == "kmeans") {
		if (action == "add"){

			.data %>%
				dplyr::left_join(
					(.) %>%
						get_clusters_kmeans_bulk(
							.abundance = !!.abundance,
							.element = !!.element,
							.feature = !!.feature,
							of_samples = of_samples,
							log_transform = log_transform,
							...
						),
					by=quo_name(.element)
				) %>%

				# Attach attributes
				reattach_internals(.data)

		}
		else if (action == "get"){

			.data %>%

				# Selecting the right columns
				select(
					!!.element,
					get_x_y_annotation_columns(.data, !!.element,!!.feature, !!.abundance, NULL)$horizontal_cols
				) %>%
				distinct() %>%

				dplyr::left_join(
					.data %>%
						get_clusters_kmeans_bulk(
							.abundance = !!.abundance,
							.element = !!.element,
							.feature = !!.feature,
							of_samples = of_samples,
							log_transform = log_transform,
							...
						),
					by=quo_name(.element)
				) %>%

				# Attach attributes
				reattach_internals(.data)

		}
		else if (action == "only")
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
				"tidybulk says: action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
			)
	}
	else if (method == "SNN") {
		if (action == "add"){

			.data %>%
				dplyr::left_join(
					(.) %>%
						get_clusters_SNN_bulk(
							.abundance = !!.abundance,
							.element = !!.element,
							.feature = !!.feature,
							of_samples = of_samples,
							log_transform = log_transform,
							...
						),
					by=quo_name(.element)
				) %>%

				# Attach attributes
				reattach_internals(.data)

		}
		else if (action == "get"){

			.data %>%

				# Selecting the right columns
				select(
					!!.element,
					get_x_y_annotation_columns(.data, !!.element,!!.feature, !!.abundance, NULL)$horizontal_cols
				) %>%
				distinct() %>%

				dplyr::left_join(
					.data %>%
						get_clusters_SNN_bulk(
							.abundance = !!.abundance,
							.element = !!.element,
							.feature = !!.feature,
							of_samples = of_samples,
							log_transform = log_transform,
							...
						),
					by=quo_name(.element)
				) %>%

				# Attach attributes
				reattach_internals(.data)

		}

		else if (action == "only")
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
				"tidybulk says: action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
			)
	}
	else
		stop("tidybulk says: the only supported methods are \"kmeans\" or \"SNN\" ")

}

#' cluster_elements
#' @inheritParams cluster_elements
#' 
#' @docType methods
#' @rdname cluster_elements-methods
#' 
#' @return A tbl object with additional columns with cluster labels
#'
setMethod("cluster_elements", "spec_tbl_df", .cluster_elements)

#' cluster_elements
#' @inheritParams cluster_elements
#' 
#' @docType methods
#' @rdname cluster_elements-methods
#' 
#' @return A tbl object with additional columns with cluster labels
#'
setMethod("cluster_elements", "tbl_df", .cluster_elements)

#' cluster_elements
#' @inheritParams cluster_elements
#' 
#' @docType methods
#' @rdname cluster_elements-methods
#' 
#' @return A tbl object with additional columns with cluster labels
#'
setMethod("cluster_elements", "tidybulk", .cluster_elements)

.cluster_elements_se = function(.data,
																.element = NULL,
																.feature = NULL,
																.abundance = NULL,
																method ,
																of_samples = TRUE,
																log_transform = TRUE,
																action = "add",
																...) {
	# Get column names
	.element = enquo(.element)
	.feature = enquo(.feature)
	.abundance = enquo(.abundance)


	.data %>%

		# Convert to tidybulk
		tidybulk() %>%

		# Apply scale method
		cluster_elements(
			.element = !!.element ,
			.feature = !!.feature ,
			.abundance = !!.abundance,
			method = method ,
			of_samples = of_samples,
			log_transform = log_transform,
			action = action,
			...
		) %>%

		# Convert to SummaizedExperiment
		tidybulk_to_SummarizedExperiment()

}

#' cluster_elements
#' @inheritParams cluster_elements
#' 
#' @docType methods
#' @rdname cluster_elements-methods
#' 
#' @return A `SummarizedExperiment` object
#'
setMethod("cluster_elements",
					"SummarizedExperiment",
					.cluster_elements_se)

#' cluster_elements
#' @inheritParams cluster_elements
#' 
#' @docType methods
#' @rdname cluster_elements-methods
#' 
#' @return A `SummarizedExperiment` object
#'
setMethod("cluster_elements",
					"RangedSummarizedExperiment",
					.cluster_elements_se)


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
#' @param of_samples A boolean. In case the input is a tidybulk object, it indicates Whether the element column will be sample or transcript column
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
#' counts.MDS =  reduce_dimensions(tidybulk::counts_mini, sample, transcript, count, method="MDS", .dims = 3)
#'
#'
#' counts.PCA =  reduce_dimensions(tidybulk::counts_mini, sample, transcript, count, method="PCA", .dims = 3)
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
																				 ...)
					 standardGeneric("reduce_dimensions"))

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
	# Get column names
	.element = enquo(.element)
	.feature = enquo(.feature)
	col_names = get_elements_features(.data, .element, .feature, of_samples)
	.element = col_names$.element
	.feature = col_names$.feature

	# Get scaled abundance if present, otherwise get abundance
	.abundance = enquo(.abundance)
	col_names = get_abundance_norm_if_exists(.data, .abundance)
	.abundance = col_names$.abundance

	# Validate data frame
	validation(.data, !!.element, !!.feature, !!.abundance)

	if (method == "MDS") {

		.data_processed =
			.data %>%
			get_reduced_dimensions_MDS_bulk(
				.abundance = !!.abundance,
				.dims = .dims,
				.element = !!.element,
				.feature = !!.feature,
				top = top,
				of_samples = of_samples,
				log_transform = log_transform,
				...
			)

		if (action == "add"){

			.data %>%	dplyr::left_join(.data_processed,	by = quo_name(.element)) %>%

				# Attach attributes
				reattach_internals(.data_processed)

		}
		else if (action == "get"){

			.data %>%

				# Selecting the right columns
				select(
					!!.element,
					get_x_y_annotation_columns(.data, !!.element,!!.feature, !!.abundance, NULL)$horizontal_cols
				) %>%
				distinct() %>%

				dplyr::left_join(.data_processed,	by = quo_name(.element)) %>%

				# Attach attributes
				reattach_internals(.data_processed)

		}

		else if (action == "only") .data_processed
		else
			stop(
				"tidybulk says: action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
			)
	}
	else if (method == "PCA") {

		.data_processed =
			.data %>%
			get_reduced_dimensions_PCA_bulk(
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

		if (action == "add"){

			.data %>%
				dplyr::left_join(.data_processed,	by = quo_name(.element)) %>%

				# Attach attributes
				reattach_internals(.data_processed)

		}

		else if (action == "get"){

			.data %>%

				# Selecting the right columns
				select(
					!!.element,
					get_x_y_annotation_columns(.data, !!.element,!!.feature, !!.abundance, NULL)$horizontal_cols
				) %>%
				distinct() %>%

				dplyr::left_join(.data_processed,	by = quo_name(.element)) %>%

				# Attach attributes
				reattach_internals(.data_processed)

		}

		else if (action == "only")	.data_processed
		else
			stop(
				"tidybulk says: action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
			)

	}
	else if (method == "tSNE") {

		.data_processed =
			.data %>%
			get_reduced_dimensions_TSNE_bulk(
				.abundance = !!.abundance,
				.dims = .dims,
				.element = !!.element,
				.feature = !!.feature,
				top = top,
				of_samples = of_samples,
				log_transform = log_transform,
				...
			)

		if (action == "add"){

			.data %>%
				dplyr::left_join(.data_processed,	by = quo_name(.element)	) %>%

				# Attach attributes
				reattach_internals(.data)

		}
		else if (action == "get"){

			.data %>%

				# Selecting the right columns
				select(
					!!.element,
					get_x_y_annotation_columns(.data, !!.element,!!.feature, !!.abundance, NULL)$horizontal_cols
				) %>%
				distinct() %>%

				dplyr::left_join(.data_processed,	by = quo_name(.element)	) %>%

				# Attach attributes
				reattach_internals(.data)

		}
		else if (action == "only") .data_processed
		else
			stop(
				"tidybulk says: action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
			)

	}
	else
		stop("tidybulk says: method must be either \"MDS\" or \"PCA\"")

}

#' reduce_dimensions
#' @inheritParams reduce_dimensions
#' 
#' @docType methods
#' @rdname reduce_dimensions-methods
#' 
#' @return A tbl object with additional columns for the reduced dimensions
setMethod("reduce_dimensions", "spec_tbl_df", .reduce_dimensions)

#' reduce_dimensions
#' @inheritParams reduce_dimensions
#' 
#' @docType methods
#' @rdname reduce_dimensions-methods
#' 
#' @return A tbl object with additional columns for the reduced dimensions
setMethod("reduce_dimensions", "tbl_df", .reduce_dimensions)

#' reduce_dimensions
#' @inheritParams reduce_dimensions
#' 
#' @docType methods
#' @rdname reduce_dimensions-methods
#' 
#' @return A tbl object with additional columns for the reduced dimensions
setMethod("reduce_dimensions", "tidybulk", .reduce_dimensions)

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
																 ...) {
	# Get column names
	.element = enquo(.element)
	.feature = enquo(.feature)
	.abundance = enquo(.abundance)

	.data %>%

		# Convert to tidybulk
		tidybulk() %>%

		# Apply scale method
		reduce_dimensions(
			.element = !!.element,
			.feature  = !!.feature,
			.abundance  = !!.abundance,
			method = method,
			.dims = .dims,

			top = top,
			of_samples = of_samples,
			log_transform = log_transform,
			scale = scale,
			action = action,
			...
		) %>%

		# Convert to SummaizedExperiment
		tidybulk_to_SummarizedExperiment()

}

#' reduce_dimensions
#' @inheritParams reduce_dimensions
#' 
#' @docType methods
#' @rdname reduce_dimensions-methods
#' 
#' @return A `SummarizedExperiment` object
#'
setMethod("reduce_dimensions",
					"SummarizedExperiment",
					.reduce_dimensions_se)

#' reduce_dimensions
#' @inheritParams reduce_dimensions
#' 
#' @docType methods
#' @rdname reduce_dimensions-methods
#' 
#' @return A `SummarizedExperiment` object
#'
setMethod("reduce_dimensions",
					"RangedSummarizedExperiment",
					.reduce_dimensions_se)

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
#' @param of_samples A boolean. In case the input is a tidybulk object, it indicates Whether the element column will be sample or transcript column
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
#' counts.MDS =  reduce_dimensions(tidybulk::counts_mini, sample, transcript, count, method="MDS", .dims = 3)
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
																				 action = "add")
	standardGeneric("rotate_dimensions"))

# Set internal
.rotate_dimensions = 		function(.data,
																dimension_1_column,
																dimension_2_column,
																rotation_degrees,
																.element = NULL,
																of_samples = TRUE,
																dimension_1_column_rotated = NULL,
																dimension_2_column_rotated = NULL,
																action =	"add")
{
	# Get column names
	.element = enquo(.element)
	col_names = get_elements(.data, .element)
	.element = col_names$.element

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

	.data_processed =
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

	if (action == "add"){

		.data %>%
			dplyr::left_join(	.data_processed,	by = quo_name(.element)	) %>%

			# Attach attributes
			reattach_internals(.data)

	}
	else if (action == "get"){

		.data %>%

			# Selecting the right columns
			select(
				!!.element,
				get_specific_annotation_columns(.data, !!.element)
			) %>%
			distinct() %>%

			dplyr::left_join(	.data_processed,	by = quo_name(.element)	) %>%

			# Attach attributes
			reattach_internals(.data)

	}
	else if (action == "only") .data_processed
	else
		stop(
			"tidybulk says: action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
		)
}

#' rotate_dimensions
#' @inheritParams rotate_dimensions
#' 
#' @docType methods
#' @rdname rotate_dimensions-methods
#' 
#' @return A tbl object with additional columns for the reduced dimensions. additional columns for the rotated dimensions. The rotated dimensions will be added to the original data set as `<NAME OF DIMENSION> rotated <ANGLE>` by default, or as specified in the input arguments.
setMethod("rotate_dimensions", "spec_tbl_df", .rotate_dimensions)

#' rotate_dimensions
#' @inheritParams rotate_dimensions
#' 
#' @docType methods
#' @rdname rotate_dimensions-methods
#' 
#' @return A tbl object with additional columns for the reduced dimensions. additional columns for the rotated dimensions. The rotated dimensions will be added to the original data set as `<NAME OF DIMENSION> rotated <ANGLE>` by default, or as specified in the input arguments.
setMethod("rotate_dimensions", "tbl_df", .rotate_dimensions)

#' rotate_dimensions
#' @inheritParams rotate_dimensions
#' 
#' @docType methods
#' @rdname rotate_dimensions-methods
#' 
#' @return A tbl object with additional columns for the reduced dimensions. additional columns for the rotated dimensions. The rotated dimensions will be added to the original data set as `<NAME OF DIMENSION> rotated <ANGLE>` by default, or as specified in the input arguments.
setMethod("rotate_dimensions", "tidybulk", .rotate_dimensions)


.rotate_dimensions_se = function(.data,
																 dimension_1_column,
																 dimension_2_column,
																 rotation_degrees,
																 .element = NULL,
																 of_samples = TRUE,
																 dimension_1_column_rotated = NULL,
																 dimension_2_column_rotated = NULL,
																 action =
																 	"add") {
	# Get column names
	.element = enquo(.element)

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

	.data %>%

		# Convert to tidybulk
		tidybulk() %>%

		# Apply scale method
		rotate_dimensions(
			dimension_1_column = !!dimension_1_column,
			dimension_2_column = !!dimension_2_column,
			rotation_degrees = rotation_degrees,
			.element = !!.element,
			of_samples = of_samples,
			dimension_1_column_rotated = !!dimension_1_column_rotated,
			dimension_2_column_rotated = !!dimension_2_column_rotated,
			action = action
		) %>%

		# Convert to SummaizedExperiment
		tidybulk_to_SummarizedExperiment()

}

#' rotate_dimensions
#' @inheritParams rotate_dimensions
#' 
#' @docType methods
#' @rdname rotate_dimensions-methods
#' 
#' @return A `SummarizedExperiment` object
#'
setMethod("rotate_dimensions",
					"SummarizedExperiment",
					.rotate_dimensions_se)

#' rotate_dimensions
#' @inheritParams rotate_dimensions
#' 
#' @docType methods
#' @rdname rotate_dimensions-methods
#' 
#' @return A `SummarizedExperiment` object
#'
setMethod("rotate_dimensions",
					"RangedSummarizedExperiment",
					.rotate_dimensions_se)


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
#' @param of_samples A boolean. In case the input is a tidybulk object, it indicates Whether the element column will be sample or transcript column
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
#'     tidybulk::counts_mini,
#' 	   .element = sample,
#' 	   .feature = transcript,
#' 	   	.abundance =  count,
#' 	   	method = "correlation"
#' 	   	)
#'
#' counts.MDS =  reduce_dimensions(tidybulk::counts_mini, sample, transcript, count, method="MDS", .dims = 3)
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
																				 Dim_b_column)
					 standardGeneric("remove_redundancy"))

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

	if (method == "correlation") {
		# Validate data frame
		validation(.data, !!.element, !!.feature, !!.abundance)

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
	else if (method == "reduced_dimensions") {
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
			"tidybulk says: method must be either \"correlation\" for dropping correlated elements or \"reduced_dimension\" to drop the closest pair according to two dimensions (e.g., PCA)"
		)

}

#' remove_redundancy
#' @inheritParams remove_redundancy
#' 
#' @docType methods
#' @rdname remove_redundancy-methods
#' 
#' @return A tbl object with with dropped recundant elements (e.g., samples).
setMethod("remove_redundancy", "spec_tbl_df", .remove_redundancy)

#' remove_redundancy
#' @inheritParams remove_redundancy
#' 
#' @docType methods
#' @rdname remove_redundancy-methods
#' 
#' @return A tbl object with with dropped recundant elements (e.g., samples).
setMethod("remove_redundancy", "tbl_df", .remove_redundancy)

#' remove_redundancy
#' @inheritParams remove_redundancy
#' 
#' @docType methods
#' @rdname remove_redundancy-methods
#' 
#' @return A tbl object with with dropped recundant elements (e.g., samples).
setMethod("remove_redundancy", "tidybulk", .remove_redundancy)

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
																 Dim_b_column = NULL) {
	# Make col names
	.abundance = enquo(.abundance)
	.element = enquo(.element)
	.feature = enquo(.feature)

	Dim_a_column = enquo(Dim_a_column)
	Dim_b_column = enquo(Dim_b_column)

	.data %>%

		# Convert to tidybulk
		tidybulk() %>%

		# Apply scale method
		remove_redundancy(
			.element = !!.element,
			.feature = !!.feature,
			.abundance = !!.abundance,
			method = method,

			of_samples = of_samples,



			correlation_threshold = correlation_threshold,
			top = top,
			log_transform = log_transform,

			Dim_a_column = !!Dim_a_column,
			Dim_b_column = !!Dim_b_column
		) %>%

		# Convert to SummaizedExperiment
		tidybulk_to_SummarizedExperiment()

}

#' remove_redundancy
#' @inheritParams remove_redundancy
#' 
#' @docType methods
#' @rdname remove_redundancy-methods
#' 
#' @return A `SummarizedExperiment` object
#'
setMethod("remove_redundancy",
					"SummarizedExperiment",
					.remove_redundancy_se)

#' remove_redundancy
#' @inheritParams remove_redundancy
#' 
#' @docType methods
#' @rdname remove_redundancy-methods
#' 
#' @return A `SummarizedExperiment` object
#'
setMethod("remove_redundancy",
					"RangedSummarizedExperiment",
					.remove_redundancy_se)

#' Adjust transcript abundance for unwanted variation
#'
#' \lifecycle{maturing}
#'
#' @description adjust_abundance() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and returns a `tbl` with an edditional adjusted abundance column. This method uses scaled counts if present.
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
#' @return A `tbl` with additional columns for the adjusted counts as `<COUNT COLUMN>_adjusted`
#'
#'
#'
#'
#' @examples
#'
#'
#'
#' cm = tidybulk::counts_mini
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
																				...)
	standardGeneric("adjust_abundance"))

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
	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	col_names = get_sample_transcript(.data, .sample, .transcript)
	.sample = col_names$.sample
	.transcript = col_names$.transcript

	# Get scaled abundance if present, otherwise get abundance
	.abundance = enquo(.abundance)
	col_names = get_abundance_norm_if_exists(.data, .abundance)
	.abundance = col_names$.abundance

	# Validate data frame
	validation(.data, !!.sample, !!.transcript, !!.abundance)

	.data_processed =
		get_adjusted_counts_for_unwanted_variation_bulk(
			.data,
			.formula,
			.sample = !!.sample,
			.transcript = !!.transcript,
			.abundance = !!.abundance,
			log_transform = log_transform,
			...
		)

	if (action == "add"){

		.data %>%

			# Add adjsted column
			dplyr::left_join(.data_processed,	by = c(quo_name(.transcript), quo_name(.sample))) %>%

			# Attach attributes
			reattach_internals(.data)

	}
	else if (action == "get"){

		.data %>%

			# Selecting the right columns
			select(
				!!.sample,
				get_x_y_annotation_columns(.data, !!.sample,!!.transcript, !!.abundance, NULL)$horizontal_cols
			) %>%
			distinct() %>%

			# Add adjsted column
			dplyr::left_join(.data_processed,	by = quo_name(.sample)) %>%

			# Attach attributes
			reattach_internals(.data)

	}
	else if (action == "only") .data_processed
	else
		stop(
			"tidybulk says: action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
		)
}

#' adjust_abundance
#' @inheritParams adjust_abundance
#' 
#' @docType methods
#' @rdname adjust_abundance-methods
#' 
#' @return A `tbl` with additional columns for the adjusted counts as `<COUNT COLUMN>_adjusted`
setMethod("adjust_abundance", "spec_tbl_df", .adjust_abundance)

#' adjust_abundance
#' @inheritParams adjust_abundance
#' 
#' @docType methods
#' @rdname adjust_abundance-methods
#' 
#' @return A `tbl` with additional columns for the adjusted counts as `<COUNT COLUMN>_adjusted`
setMethod("adjust_abundance", "tbl_df", .adjust_abundance)

#' adjust_abundance
#' @inheritParams adjust_abundance
#' 
#' @docType methods
#' @rdname adjust_abundance-methods
#' 
#' @return A `tbl` with additional columns for the adjusted counts as `<COUNT COLUMN>_adjusted`
setMethod("adjust_abundance", "tidybulk", .adjust_abundance)

.adjust_abundance_se = function(.data,
																.formula,
																.sample = NULL,
																.transcript = NULL,
																.abundance = NULL,
																log_transform = TRUE,
																action = "add",
																...) {
	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	.data %>%

		# Convert to tidybulk
		tidybulk() %>%

		# Apply scale method
		adjust_abundance(
			.formula = .formula,
			.sample = !!.sample,
			.transcript = !!.transcript,
			.abundance = !!.abundance,
			log_transform = log_transform,
			action = action,
			...
		) %>%

		# Convert to SummaizedExperiment
		tidybulk_to_SummarizedExperiment()

}

#' adjust_abundance
#' @inheritParams adjust_abundance
#' 
#' @docType methods
#' @rdname adjust_abundance-methods
#' 
#' @return A `SummarizedExperiment` object
#'
setMethod("adjust_abundance",
					"SummarizedExperiment",
					.adjust_abundance_se)

#' adjust_abundance
#' @inheritParams adjust_abundance
#' 
#' @docType methods
#' @rdname adjust_abundance-methods
#' 
#' @return A `SummarizedExperiment` object
#'
setMethod("adjust_abundance",
					"RangedSummarizedExperiment",
					.adjust_abundance_se)



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
#'     tidybulk::counts_mini,
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
																						keep_integer = TRUE)
	standardGeneric("aggregate_duplicates"))

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
	validation(.data,
						 !!.sample,
						 !!.transcript,
						 !!.abundance,
						 skip_dupli_check = TRUE)

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
#' 
#' @docType methods
#' @rdname aggregate_duplicates-methods
#' 
#' @return A `tbl` object with aggregated transcript abundance and annotation
setMethod("aggregate_duplicates", "spec_tbl_df", .aggregate_duplicates)

#' aggregate_duplicates
#' @inheritParams aggregate_duplicates
#' 
#' @docType methods
#' @rdname aggregate_duplicates-methods
#' 
#' @return A `tbl` object with aggregated transcript abundance and annotation
setMethod("aggregate_duplicates", "tbl_df", .aggregate_duplicates)

#' aggregate_duplicates
#' @inheritParams aggregate_duplicates
#' 
#' @docType methods
#' @rdname aggregate_duplicates-methods
#' 
#' @return A `tbl` object with aggregated transcript abundance and annotation
setMethod("aggregate_duplicates", "tidybulk", .aggregate_duplicates)

.aggregate_duplicates_se = function(.data,

																		.sample = NULL,
																		.transcript = NULL,
																		.abundance = NULL,
																		aggregation_function = sum,
																		keep_integer = TRUE) {
	# Make col names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	.data %>%

		# Convert to tidybulk
		tidybulk() %>%

		# Apply scale method
		aggregate_duplicates(
			.sample = !!.sample,
			.transcript = !!.transcript,
			.abundance = !!.abundance,
			aggregation_function = aggregation_function,
			keep_integer = keep_integer
		) %>%

		# Convert to SummaizedExperiment
		tidybulk_to_SummarizedExperiment()

}

#' aggregate_duplicates
#' @inheritParams aggregate_duplicates
#' 
#' @docType methods
#' @rdname aggregate_duplicates-methods
#' 
#' @return A `SummarizedExperiment` object
#'
setMethod("aggregate_duplicates",
					"SummarizedExperiment",
					.aggregate_duplicates_se)

#' aggregate_duplicates
#' @inheritParams aggregate_duplicates
#' 
#' @docType methods
#' @rdname aggregate_duplicates-methods
#' 
#' @return A `SummarizedExperiment` object
#'
setMethod("aggregate_duplicates",
					"RangedSummarizedExperiment",
					.aggregate_duplicates_se)



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
#' deconvolve_cellularity(tidybulk::counts, sample, transcript, `count`, cores = 2)
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
																							...)
	standardGeneric("deconvolve_cellularity"))

# Set internal
.deconvolve_cellularity = 		function(.data,
																		 .sample = NULL,
																		 .transcript = NULL,
																		 .abundance = NULL,
																		 reference = X_cibersort,
																		 method = "cibersort",
																		 action = "add",
																		 ...)  {
	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
	.sample = col_names$.sample
	.transcript = col_names$.transcript
	.abundance = col_names$.abundance

	# Validate data frame
	validation(.data, !!.sample, !!.transcript, !!.abundance)

	.data_processed =
		get_cell_type_proportions(
		.data,
		.sample = !!.sample,
		.transcript = !!.transcript,
		.abundance = !!.abundance,
		reference = reference,
		method = method,
		...
	)

	if (action == "add"){
		.data %>%

			# Add new annotation
			dplyr::left_join(.data_processed,				by = quo_name(.sample)			) %>%

			# Attach attributes
			reattach_internals(.data)
	}

	else if (action == "get"){
		.data %>%


			# Selecting the right columns
			select(
				!!.sample,
				get_x_y_annotation_columns(.data, !!.sample,!!.transcript, !!.abundance, NULL)$horizontal_cols
			) %>%
			distinct() %>%

			# Add new annotation
			dplyr::left_join(.data_processed,				by = quo_name(.sample)			) %>%

			# Attach attributes
			reattach_internals(.data)
	}

	else if (action == "only") .data_processed
	else
		stop(
			"tidybulk says: action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
		)
}

#' deconvolve_cellularity
#' @inheritParams deconvolve_cellularity
#' 
#' @docType methods
#' @rdname deconvolve_cellularity-methods
#' 
#' @return A `tbl` object including additional columns for each cell type estimated
setMethod("deconvolve_cellularity",
					"spec_tbl_df",
					.deconvolve_cellularity)

#' deconvolve_cellularity
#' @inheritParams deconvolve_cellularity
#' 
#' @docType methods
#' @rdname deconvolve_cellularity-methods
#' 
#' @return A `tbl` object including additional columns for each cell type estimated
setMethod("deconvolve_cellularity", "tbl_df", .deconvolve_cellularity)

#' deconvolve_cellularity
#' @inheritParams deconvolve_cellularity
#' 
#' @docType methods
#' @rdname deconvolve_cellularity-methods
#' 
#' @return A `tbl` object including additional columns for each cell type estimated
setMethod("deconvolve_cellularity",
					"tidybulk",
					.deconvolve_cellularity)



.deconvolve_cellularity_se = function(.data,
																			.sample = NULL,
																			.transcript = NULL,
																			.abundance = NULL,
																			reference = X_cibersort,
																			method = "cibersort",
																			action = "add",
																			...) {
	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	.data %>%

		# Convert to tidybulk
		tidybulk() %>%

		# Apply scale method
		deconvolve_cellularity(
			.sample = !!.sample,
			.transcript = !!.transcript,
			.abundance = !!.abundance,
			reference = reference,
			method = method,
			action = action,
			...
		) %>%

		# Convert to SummaizedExperiment
		tidybulk_to_SummarizedExperiment()

}

#' deconvolve_cellularity
#' @inheritParams deconvolve_cellularity
#' 
#' @docType methods
#' @rdname deconvolve_cellularity-methods
#' 
#' @return A `SummarizedExperiment` object
#'
setMethod("deconvolve_cellularity",
					"SummarizedExperiment",
					.deconvolve_cellularity_se)

#' deconvolve_cellularity
#' @inheritParams deconvolve_cellularity
#' 
#' @docType methods
#' @rdname deconvolve_cellularity-methods
#' 
#' @return A `SummarizedExperiment` object
#'
setMethod(
	"deconvolve_cellularity",
	"RangedSummarizedExperiment",
	.deconvolve_cellularity_se
)

#' Get ENTREZ id from gene SYMBOL
#'
#' @param .data A tt or tbl object.
#' @param .transcript A character. The name of the ene symbol column.
#' @param .sample The name of the sample column
#'
#' @return A tbl
#'
#' @examples
#'
#' symbol_to_entrez(tidybulk::counts_mini, .transcript = transcript, .sample = sample)
#'
#' @export
#'
symbol_to_entrez = function(.data,
														.transcript = NULL,
														.sample = NULL) {
	# Get column names
	.transcript = enquo(.transcript)
	.sample = enquo(.sample)
	col_names = get_sample_transcript(.data, .sample, .transcript)
	.transcript = col_names$.transcript
	
	# Check if package is installed, otherwise install
	if (find.package("org.Hs.eg.db", quiet = TRUE) %>% length %>% equals(0)) {
		message("Installing org.Hs.eg.db needed for differential transcript abundance analyses")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("org.Hs.eg.db", ask = FALSE)
	}
	
	.data %>%
		dplyr::left_join(
			# Get entrez mapping 1:1
			AnnotationDbi::mapIds(
				org.Hs.eg.db::org.Hs.eg.db,
				.data %>% distinct(!!.transcript) %>% pull(!!.transcript) %>% as.character,
				'ENTREZID',
				'SYMBOL'
			) %>%
				enframe(name = quo_name(.transcript), value = "entrez") %>%
				filter(entrez %>% is.na %>% `!`) %>%
				group_by(!!.transcript) %>%
				slice(1) %>%
				ungroup(),
			by = quo_name(.transcript)
		)
	
}


#' Add transcript symbol column from ensembl id for human and mouse data
#'
#' \lifecycle{maturing}
#'
#' @description ensembl_to_symbol() takes as imput a `tbl` formatted as | <SAMPLE> | <ENSEMBL_ID> | <COUNT> | <...> | and returns a `tbl` with the additional transcript symbol column
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name ensembl_to_symbol
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <ENSEMBL_ID> | <COUNT> | <...> |
#' @param .ensembl A character string. The column that is represents ensembl gene id
#'
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#'
#' @details This is useful since different resources use ensembl IDs while others use gene symbol IDs. At the moment this work for human (genes and transcripts) and mouse (genes) data.
#'
#' @return A `tbl` object including additional columns for transcript symbol
#'
#'
#'
#'
#' @examples
#'
#'
#' 	ensembl_to_symbol(tidybulk::counts_ensembl, ens)
#'
#'
#' @docType methods
#' @rdname ensembl_to_symbol-methods
#' @export
#'
#'
setGeneric("ensembl_to_symbol", function(.data,
																			 .ensembl,
																			 action = "add")
	standardGeneric("ensembl_to_symbol"))

# Set internal
.ensembl_to_symbol = 		function(.data,
															.ensembl,
															action = "add")
{
	# Make col names
	.ensembl = enquo(.ensembl)

	.data_processed = get_symbol_from_ensembl(.data,!!.ensembl)

	if (action == "add"){

		# Add new symbols column
		.data %>%
			dplyr::left_join(.data_processed, by=quo_name(.ensembl)) %>%

			# Attach attributes
			reattach_internals(.data)

	}
	# else if (action == "get"){
	#
	# 	# Add new symbols column
	# 	.data %>%
	#
	#
	# 		dplyr::left_join(.data_processed) %>%
	#
	# 		# Attach attributes
	# 		reattach_internals(.data)
	#
	# }

	else if (action == "only") .data_processed

	else
		stop(
			"tidybulk says: action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
		)

}

#' ensembl_to_symbol
#' @inheritParams ensembl_to_symbol
#' 
#' @docType methods
#' @rdname ensembl_to_symbol-methods
#' 
#' @return A `tbl` object including additional columns for transcript symbol
setMethod("ensembl_to_symbol", "spec_tbl_df", .ensembl_to_symbol)

#' ensembl_to_symbol
#' @inheritParams ensembl_to_symbol
#' 
#' @docType methods
#' @rdname ensembl_to_symbol-methods
#' 
#' @return A `tbl` object including additional columns for transcript symbol
setMethod("ensembl_to_symbol", "tbl_df", .ensembl_to_symbol)

#' ensembl_to_symbol
#' @inheritParams ensembl_to_symbol
#' 
#' @docType methods
#' @rdname ensembl_to_symbol-methods
#' 
#' @return A `tbl` object including additional columns for transcript symbol
setMethod("ensembl_to_symbol", "tidybulk", .ensembl_to_symbol)


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
#' @param .contrasts A character vector. See edgeR makeContrasts specification for the parameter `contrasts`. If contrasts are not present the first covariate is the one the model is tested against (e.g., ~ factor_of_interest)
#' @param method A string character. Either "edgeR_quasi_likelihood" (i.e., QLF), "edgeR_likelihood_ratio" (i.e., LRT)
#' @param significance_threshold A real between 0 and 1 (usually 0.05).
#' @param minimum_counts A real positive number. It is the threshold of count per million that is used to filter transcripts/genes out from the scaling procedure.
#' @param minimum_proportion A real positive number between 0 and 1. It is the threshold of proportion of samples for each transcripts/genes that have to be characterised by a cmp bigger than the threshold to be included for scaling procedure.
#' @param fill_missing_values A boolean. Whether to fill missing sample/transcript values with the median of the transcript. This is rarely needed.
#' @param scaling_method A character string. The scaling method passed to the backend function (i.e., edgeR::calcNormFactors; "TMM","TMMwsp","RLE","upperquartile")
#' @param omit_contrast_in_colnames If just one contrast is specified you can choose to omit the contrast label in the colnames.
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
#' 	 tidybulk::counts_mini,
#' 	    ~ condition,
#' 	    sample,
#' 	    transcript,
#' 	    `count`
#' 	)
#'
#' 	# The functon `test_differential_abundance` operated with contrasts too
#'
#'  test_differential_abundance(
#' 	    tidybulk::counts_mini,
#' 	    ~ 0 + condition,
#' 	    sample,
#' 	    transcript,
#' 	    `count`,
#' 	    .contrasts = c( "conditionTRUE - conditionFALSE")
#'  )
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
																									 .contrasts = NULL,
																									 method = "edgeR_quasi_likelihood",
																									 significance_threshold = 0.05,
																									 minimum_counts = 10,
																									 minimum_proportion = 0.7,
																									 fill_missing_values = FALSE,
																									 scaling_method = "TMM",
																									 omit_contrast_in_colnames = FALSE,

																									 action = "add")
					 standardGeneric("test_differential_abundance"))

# Set internal
.test_differential_abundance = 		function(.data,
																					.formula,
																					.sample = NULL,
																					.transcript = NULL,
																					.abundance = NULL,
																					.contrasts = NULL,
																					method = "edgeR_quasi_likelihood",
																					significance_threshold = 0.05,
																					minimum_counts = 10,
																					minimum_proportion = 0.7,
																					fill_missing_values = FALSE,
																					scaling_method = "TMM",
																					omit_contrast_in_colnames = FALSE,

																					action = "add")
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
	validation(.data, !!.sample, !!.transcript, !!.abundance)

	if(grepl("edgeR", method)){
		.data_processed =
			get_differential_transcript_abundance_bulk(
				.data,
				.formula,
				.sample = !!.sample,
				.transcript = !!.transcript,
				.abundance = !!.abundance,
				.contrasts = .contrasts,
				method = method,
				significance_threshold = significance_threshold,
				minimum_counts = minimum_counts,
				minimum_proportion = minimum_proportion,
				fill_missing_values = fill_missing_values,
				scaling_method = scaling_method,
				omit_contrast_in_colnames = omit_contrast_in_colnames
			)
	}
	else stop("tidybulk says: the onyl methods supported at the moment are \"edgeR_quasi_likelihood\" (i.e., QLF), \"edgeR_likelihood_ratio\" (i.e., LRT)")

	if (action == "add"){

		.data %>%
			dplyr::left_join(.data_processed, by = quo_name(.transcript)) %>%

			# Arrange
			ifelse_pipe(.contrasts %>% is.null,
									~ .x %>% arrange(FDR))	%>%

			# Attach attributes
			reattach_internals(.data_processed)

	}
	else if (action == "get"){

		.data %>%

			# Selecting the right columns
			select(
				!!.transcript,
				get_x_y_annotation_columns(.data, !!.sample,!!.transcript, !!.abundance, NULL)$vertical_cols
			) %>%
			distinct() %>%

			dplyr::left_join(.data_processed, by = quo_name(.transcript)) %>%

			# Arrange
			ifelse_pipe(.contrasts %>% is.null,
									~ .x %>% arrange(FDR))	%>%

			# Attach attributes
			reattach_internals(.data_processed)

	}
	else if (action == "only") .data_processed
	else
		stop(
			"tidybulk says: action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
		)
}

#' test_differential_abundance
#' @inheritParams test_differential_abundance
#' 
#' @docType methods
#' @rdname test_differential_abundance-methods
#' 
#' @return A `tbl` with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
setMethod("test_differential_abundance",
					"spec_tbl_df",
					.test_differential_abundance)

#' test_differential_abundance
#' @inheritParams test_differential_abundance
#' 
#' @docType methods
#' @rdname test_differential_abundance-methods
#' 
#' @return A `tbl` with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
setMethod("test_differential_abundance",
					"tbl_df",
					.test_differential_abundance)

#' test_differential_abundance
#' @inheritParams test_differential_abundance
#' 
#' @docType methods
#' @rdname test_differential_abundance-methods
#' 
#' @return A `tbl` with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
setMethod("test_differential_abundance",
					"tidybulk",
					.test_differential_abundance)



.test_differential_abundance_se = function(.data,
																					 .formula,
																					 .sample = NULL,
																					 .transcript = NULL,
																					 .abundance = NULL,
																					 .contrasts = NULL,
																					 method = "edgeR_quasi_likelihood",
																					 significance_threshold = 0.05,
																					 minimum_counts = 10,
																					 minimum_proportion = 0.7,
																					 fill_missing_values = FALSE,
																					 scaling_method = "TMM",
																					 omit_contrast_in_colnames = FALSE,
																					 action = "add")
{
	# Make col names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	.data %>%

		# Convert to tidybulk
		tidybulk() %>%

		# Apply scale method
		test_differential_abundance(
			.formula,
			.sample = !!.sample,
			.transcript = !!.transcript,
			.abundance = !!.abundance,
			.contrasts = .contrasts,
			method = method,
			significance_threshold = significance_threshold,
			minimum_counts = minimum_counts,
			minimum_proportion = minimum_proportion,
			fill_missing_values = fill_missing_values,
			scaling_method = scaling_method,
			omit_contrast_in_colnames = omit_contrast_in_colnames,
			action = action
		) %>%

		# Convert to SummaizedExperiment
		tidybulk_to_SummarizedExperiment()

}

#' test_differential_abundance
#' @inheritParams test_differential_abundance
#' 
#' @docType methods
#' @rdname test_differential_abundance-methods
#' 
#' @return A `SummarizedExperiment` object
#'
setMethod(
	"test_differential_abundance",
	"SummarizedExperiment",
	.test_differential_abundance_se
)

#' test_differential_abundance
#' @inheritParams test_differential_abundance
#' 
#' @docType methods
#' @rdname test_differential_abundance-methods
#' 
#' @return A `SummarizedExperiment` object
#'
setMethod(
	"test_differential_abundance",
	"RangedSummarizedExperiment",
	.test_differential_abundance_se
)



#' Keep variable transcripts
#'
#' \lifecycle{maturing}
#'
#' @description keep_variable() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and returns a `tbl` with additional columns for the statistics from the hypothesis test.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name keep_variable
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
#' 	keep_variable(
#' 	tidybulk::counts_mini,
#' 	    sample,
#' 	    transcript,
#' 	    `count`,
#' 	    top = 500
#' 	)
#'
#'
#' @docType methods
#' @rdname keep_variable-methods
#' @export
#'
setGeneric("keep_variable", function(.data,
																			 .sample = NULL,
																			 .transcript = NULL,
																			 .abundance = NULL,
																			 top = 500,
																			 log_transform = TRUE)
	standardGeneric("keep_variable"))

# Set internal
.keep_variable = 		function(.data,
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
	validation(.data, !!.sample, !!.transcript, !!.abundance)

	keep_variable_transcripts(
		.data,
		.sample = !!.sample,
		.transcript = !!.transcript,
		.abundance = !!.abundance,
		top = top,
		log_transform = log_transform
	)
}

#' keep_variable
#' @inheritParams keep_variable
#' 
#' @docType methods
#' @rdname keep_variable-methods
#' 
#' @return A `tbl` with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
setMethod("keep_variable", "spec_tbl_df", .keep_variable)

#' keep_variable
#' @inheritParams keep_variable
#' 
#' @docType methods
#' @rdname keep_variable-methods
#' 
#' @return A `tbl` with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
setMethod("keep_variable", "tbl_df", .keep_variable)

#' keep_variable
#' @inheritParams keep_variable
#' 
#' @docType methods
#' @rdname keep_variable-methods
#' 
#' @return A `tbl` with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
setMethod("keep_variable", "tidybulk", .keep_variable)

.keep_variable_se = function(.data,
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

		# Convert to tidybulk
		tidybulk() %>%

		# Apply scale method
		keep_variable(
			.sample = !!.sample,
			.transcript = !!.transcript,
			.abundance = !!.abundance,
			top = top,
			log_transform = log_transform
		) %>%

		# Convert to SummaizedExperiment
		tidybulk_to_SummarizedExperiment()

}

#' keep_variable
#' @inheritParams keep_variable
#' 
#' @docType methods
#' @rdname keep_variable-methods
#' 
#' @return A `SummarizedExperiment` object
#'
setMethod("keep_variable",
					"SummarizedExperiment",
					.keep_variable_se)

#' keep_variable
#' @inheritParams keep_variable
#' 
#' @docType methods
#' @rdname keep_variable-methods
#' 
#' @return A `SummarizedExperiment` object
#'
setMethod("keep_variable",
					"RangedSummarizedExperiment",
					.keep_variable_se)



#' Keep abundant transcripts
#'
#' \lifecycle{maturing}
#'
#' @description keep_abundant() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and returns a `tbl` with additional columns for the statistics from the hypothesis test.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#' @importFrom dplyr filter
#'
#' @name keep_abundant
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param factor_of_interest The name of the column of the factor of interest. This is used for defining sample groups for the filtering process.
#' @param minimum_counts A real positive number. It is the threshold of count per million that is used to filter transcripts/genes out from the scaling procedure.
#' @param minimum_proportion A real positive number between 0 and 1. It is the threshold of proportion of samples for each transcripts/genes that have to be characterised by a cmp bigger than the threshold to be included for scaling procedure.
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
#' 	keep_abundant(
#' 	tidybulk::counts_mini,
#' 	    sample,
#' 	    transcript,
#' 	    `count`
#' 	)
#'
#'
#' @docType methods
#' @rdname keep_abundant-methods
#' @export
#'
setGeneric("keep_abundant", function(.data,
																			 .sample = NULL,
																			 .transcript = NULL,
																			 .abundance = NULL,
																			 factor_of_interest = NULL,
																			 minimum_counts = 10,
																			 minimum_proportion = 0.7)
	standardGeneric("keep_abundant"))

# Set internal
.keep_abundant = 		function(.data,
															.sample = NULL,
															.transcript = NULL,
															.abundance = NULL,
															factor_of_interest = NULL,
															minimum_counts = 10,
															minimum_proportion = 0.7)
{
	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
	.sample = col_names$.sample
	.transcript = col_names$.transcript
	.abundance = col_names$.abundance

	factor_of_interest = enquo(factor_of_interest)

	# Validate data frame
	validation(.data, !!.sample, !!.transcript, !!.abundance)


	.data %>%

		# Filter
		ifelse_pipe("lowly_abundant" %in% colnames((.)),

								# If column is present use this instead of doing more work
								~ {
									#message("tidybulk says: \"lowly_abundant\" column is present. Using this to filter data")
									.x %>% dplyr::filter(!lowly_abundant)
								},

								# If no column is present go
								~ {
									gene_to_exclude =
										add_scaled_counts_bulk.get_low_expressed(
											.data,
											.sample = !!.sample,
											.transcript = !!.transcript,
											.abundance = !!.abundance,
											factor_of_interest = !!factor_of_interest,
											minimum_counts = minimum_counts,
											minimum_proportion = minimum_proportion
										)

									.x %>% dplyr::filter(!!.transcript %in% gene_to_exclude %>% `!`)
								})	%>%

		# Attach attributes
		reattach_internals(.data)
}

#' keep_abundant
#' @inheritParams keep_abundant
#' 
#' @docType methods
#' @rdname keep_abundant-methods
#' 
#' @return A `tbl` with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
setMethod("keep_abundant", "spec_tbl_df", .keep_abundant)

#' keep_abundant
#' @inheritParams keep_abundant
#' 
#' @docType methods
#' @rdname keep_abundant-methods
#' 
#' @return A `tbl` with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
setMethod("keep_abundant", "tbl_df", .keep_abundant)

#' keep_abundant
#' @inheritParams keep_abundant
#' 
#' @docType methods
#' @rdname keep_abundant-methods
#' 
#' @return A `tbl` with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
setMethod("keep_abundant", "tidybulk", .keep_abundant)

.keep_abundant_se = function(.data,
															 .sample = NULL,
															 .transcript = NULL,
															 .abundance = NULL,
															 factor_of_interest = NULL,
															 minimum_counts = 10,
															 minimum_proportion = 0.7)
{
	# Make col names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	factor_of_interest = enquo(factor_of_interest)

	.data %>%

		# Convert to tidybulk
		tidybulk() %>%

		# Apply scale method
		keep_abundant(
			.sample = !!.sample,
			.transcript = !!.transcript,
			.abundance = !!.abundance,
			factor_of_interest = !!factor_of_interest,
			minimum_counts = minimum_counts,
			minimum_proportion = minimum_proportion
		) %>%

		# Convert to SummaizedExperiment
		tidybulk_to_SummarizedExperiment()

}

#' keep_abundant
#' @inheritParams keep_abundant
#' 
#' @docType methods
#' @rdname keep_abundant-methods
#' 
#' @return A `SummarizedExperiment` object
#'
setMethod("keep_abundant",
					"SummarizedExperiment",
					.keep_abundant_se)

#' keep_abundant
#' @inheritParams keep_abundant
#' 
#' @docType methods
#' @rdname keep_abundant-methods
#' 
#' @return A `SummarizedExperiment` object
#'
setMethod("keep_abundant",
					"RangedSummarizedExperiment",
					.keep_abundant_se)



#' analyse gene enrichment with EGSEA
#'
#' \lifecycle{maturing}
#'
#' @description test_gene_enrichment() takes as imput a `tbl` formatted as | <SAMPLE> | <ENSEMBL_ID> | <COUNT> | <...> | and returns a `tbl` with the additional transcript symbol column
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name test_gene_enrichment
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .formula A formula with no response variable, representing the desired linear model
#' @param .sample The name of the sample column
#' @param .entrez The ENTREZ ID of the transcripts/genes
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
#' \dontrun{
#'
#' df_entrez = symbol_to_entrez(tidybulk::counts_mini, .transcript = transcript, .sample = sample)
#' df_entrez = aggregate_duplicates(df_entrez, aggregation_function = sum, .sample = sample, .transcript = entrez, .abundance = count)
#'
#' library("EGSEA")
#'
#' 	test_gene_enrichment(
#'			df_entrez,
#'			~ condition,
#'			.sample = sample,
#'			.entrez = entrez,
#'			.abundance = count,
#'			species="human",
#'			cores = 1
#'		)
#'
#'}
#'
#' @docType methods
#' @rdname test_gene_enrichment-methods
#' @export
#'
#'
setGeneric("test_gene_enrichment", function(.data,
																							 .formula,
																							 .sample = NULL,
																							 .entrez,
																							 .abundance = NULL,
																							 .contrasts = NULL,
																							 species,
																							 cores = 10)
	standardGeneric("test_gene_enrichment"))

# Set internal
.test_gene_enrichment = 		function(.data,
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
	col_names = get_sample_counts(.data, .sample, .abundance)
	.sample = col_names$.sample
	.abundance = col_names$.abundance
	
	.entrez = enquo(.entrez)

	# Validate data frame
	validation(.data, !!.sample, !!.entrez, !!.abundance)

	test_gene_enrichment_bulk_EGSEA(
		.data,
		.formula,
		.sample = !!.sample,
		.entrez = !!.entrez,
		.abundance = !!.abundance,
		.contrasts = .contrasts,
		species = species,
		cores = cores
	)



}

#' test_gene_enrichment
#' @inheritParams test_gene_enrichment
#'
#' @docType methods
#' @rdname test_gene_enrichment-methods
#' 
#' @return A `tbl` object
setMethod("test_gene_enrichment",
					"spec_tbl_df",
					.test_gene_enrichment)

#' test_gene_enrichment
#' @inheritParams test_gene_enrichment
#' 
#' @docType methods
#' @rdname test_gene_enrichment-methods
#' 
#' @return A `tbl` object
setMethod("test_gene_enrichment",
					"tbl_df",
					.test_gene_enrichment)

#' test_gene_enrichment
#' @inheritParams test_gene_enrichment
#' 
#' @docType methods
#' @rdname test_gene_enrichment-methods
#' 
#' @return A `tbl` object
setMethod("test_gene_enrichment",
					"tidybulk",
					.test_gene_enrichment)

#' analyse gene over-representation with GSEA
#'
#' \lifecycle{maturing}
#'
#' @description test_gene_overrepresentation() takes as imput a `tbl` formatted as | <SAMPLE> | <ENSEMBL_ID> | <COUNT> | <...> | and returns a `tbl` with the GSEA statistics
#'
#' @importFrom rlang enquo
#' @importFrom rlang quo_is_missing
#' @importFrom magrittr "%>%"
#'
#' @name test_gene_overrepresentation
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .sample The name of the sample column
#' @param .entrez The ENTREZ ID of the transcripts/genes
#' @param .do_test A boolean column name symbol. It indicates the transcript to check
#' @param species A character. For example, human or mouse. MSigDB uses the latin species names (e.g., \"Mus musculus\", \"Homo sapiens\")
#'
#' @details This wrapper execute gene enrichment analyses of the dataset using a list of transcripts and GSEA. This wrapper uses clusterProfiler on the backend.
#'
#' @return A `tbl` object
#'
#'
#'
#'
#' @examples
#'
#' df_entrez = symbol_to_entrez(tidybulk::counts_mini, .transcript = transcript, .sample = sample)
#' df_entrez = aggregate_duplicates(df_entrez, aggregation_function = sum, .sample = sample, .transcript = entrez, .abundance = count)
#' df_entrez = mutate(df_entrez, do_test = transcript %in% c("TNFRSF4", "PLCH2", "PADI4", "PAX7"))
#'
#' 	test_gene_overrepresentation(
#'			df_entrez,
#'			.sample = sample,
#'			.entrez = entrez,
#'			.do_test = do_test,
#'			species="Homo sapiens"
#'		)
#'
#'
#' @docType methods
#' @rdname test_gene_overrepresentation-methods
#' @export
#'
#'
setGeneric("test_gene_overrepresentation", function(.data,
																										.sample = NULL,
																										.entrez,
																										.do_test,
																										species)
	standardGeneric("test_gene_overrepresentation"))

# Set internal
.test_gene_overrepresentation = 		function(.data,
																					 .sample = NULL,
																					 .entrez,
																					 .do_test,
																					 species)	{
	
	# Comply with CRAN NOTES
	. = NULL
	

	# Get column names
	.sample = enquo(.sample)
	.sample =  get_sample(.data, .sample)$.sample
	.do_test = enquo(.do_test)
	.entrez = enquo(.entrez)
	
	# Check if entrez is set
	if(quo_is_missing(.entrez))
		stop("tidybulk says: the .entrez parameter appears to no be set")
	
	# Check column type
	if (.data %>% distinct(!!.do_test) %>% sapply(class) %in% c("logical") %>% `!` %>% any)
		stop("tidybulk says: .do_test column must be logical (i.e., TRUE or FALSE)")
	
	# Check packages msigdbr
	# Check if package is installed, otherwise install
	if (find.package("msigdbr", quiet = TRUE) %>% length %>% equals(0)) {
		message("msigdbr not installed. Installing.")
		BiocManager::install("msigdbr", ask = FALSE)
	}
	
	# Check is correct species name
	if(species %in% msigdbr::msigdbr_show_species() %>% `!`)
		stop(sprintf("tidybulk says: wrong species name. MSigDB uses the latin species names (e.g., %s)", paste(msigdbr::msigdbr_show_species(), collapse=", ")))
	
	#m_df <- msigdbr(species = species)
	
	.data %>% 
		#filter(!!.entrez %in% unique(m_df$entrez_gene)) %>%
		filter(!!.do_test) %>% 
		distinct(!!.entrez) %>% 
		pull(!!.entrez) %>%
		entrez_rank_to_gsea(species)
	
	
}

#' test_gene_overrepresentation
#' @inheritParams test_gene_overrepresentation
#' 
#' @docType methods
#' @rdname test_gene_overrepresentation-methods
#' 
#' @return A `tbl` object
setMethod("test_gene_overrepresentation",
					"spec_tbl_df",
					.test_gene_overrepresentation)

#' test_gene_overrepresentation
#' @inheritParams test_gene_overrepresentation
#' 
#' @docType methods
#' @rdname test_gene_overrepresentation-methods
#' 
#' @return A `tbl` object
setMethod("test_gene_overrepresentation",
					"tbl_df",
					.test_gene_overrepresentation)

#' test_gene_overrepresentation
#' @inheritParams test_gene_overrepresentation
#' 
#' @docType methods
#' @rdname test_gene_overrepresentation-methods
#' 
#' @return A `tbl` object
setMethod("test_gene_overrepresentation",
					"tidybulk",
					.test_gene_overrepresentation)


#' Extract sample-wise information
#'
#' \lifecycle{maturing}
#'
#' @description pivot_sample() takes as imput a `tbl` formatted as | <SAMPLE> | <ENSEMBL_ID> | <COUNT> | <...> | and returns a `tbl` with only sample-related columns
#'
#' @importFrom magrittr "%>%"
#'
#' @name pivot_sample
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .sample The name of the sample column
#'
#'
#' @details This functon extracts only sample-related information for downstream analysis (e.g., visualisation). It is disruptive in the sense that it cannot be passed anymore to tidybulk function.
#'
#' @return A `tbl` object
#'
#'
#'
#'
#' @examples
#'
#'
#' 	pivot_sample(
#'			tidybulk::counts_mini,
#'			.sample = sample
#'		)
#'
#'
#' @docType methods
#' @rdname pivot_sample-methods
#' @export
#'
#'
setGeneric("pivot_sample", function(.data,
																						.sample = NULL)
	standardGeneric("pivot_sample"))

# Set internal
.pivot_sample = 		function(.data,
																	 .sample = NULL)	{
	# Make col names
	.sample = enquo(.sample)
	col_names = get_sample(.data, .sample)
	.sample = col_names$.sample

	.data %>%

		# Selecting the right columns
		select(
			!!.sample,
			get_specific_annotation_columns(.data, !!.sample)
		) %>%
		distinct() %>%

		drop_class(c("tidybulk", "tt")) %>%
		drop_internals()


}

#' pivot_sample
#' @inheritParams pivot_sample
#' 
#' @docType methods
#' @rdname pivot_sample-methods
#' 
#' @return A `tbl` object
setMethod("pivot_sample",
					"spec_tbl_df",
					.pivot_sample)

#' pivot_sample
#' @inheritParams pivot_sample
#' 
#' @docType methods
#' @rdname pivot_sample-methods
#' 
#' @return A `tbl` object
setMethod("pivot_sample",
					"tbl_df",
					.pivot_sample)

#' pivot_sample
#' @inheritParams pivot_sample
#' 
#' @docType methods
#' @rdname pivot_sample-methods
#' 
#' @return A `tbl` object
setMethod("pivot_sample",
					"tidybulk",
					.pivot_sample)

#' Extract transcript-wise information
#'
#' \lifecycle{maturing}
#'
#' @description pivot_transcript() takes as imput a `tbl` formatted as | <SAMPLE> | <ENSEMBL_ID> | <COUNT> | <...> | and returns a `tbl` with only sample-related columns
#'
#' @importFrom magrittr "%>%"
#'
#' @name pivot_transcript
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .transcript The name of the transcript column
#'
#'
#' @details This functon extracts only transcript-related information for downstream analysis (e.g., visualisation). It is disruptive in the sense that it cannot be passed anymore to tidybulk function.
#'
#' @return A `tbl` object
#'
#'
#'
#'
#' @examples
#'
#'
#' 	pivot_transcript(
#'			tidybulk::counts_mini,
#'			.transcript = transcript
#'		)
#'
#'
#' @docType methods
#' @rdname pivot_transcript-methods
#' @export
#'
#'
setGeneric("pivot_transcript", function(.data,
																		.transcript = NULL)
	standardGeneric("pivot_transcript"))

# Set internal
.pivot_transcript = 		function(.data,
													 .transcript = NULL)	{
	# Make col names
	.transcript = enquo(.transcript)
	col_names = get_transcript(.data, .transcript)
	.transcript = col_names$.transcript

	.data %>%

		# Selecting the right columns
		select(
			!!.transcript,
			get_specific_annotation_columns(.data, !!.transcript)
		) %>%
		distinct() %>%

		drop_class(c("tidybulk", "tt")) %>%
		drop_internals()


}

#' pivot_transcript
#' @inheritParams pivot_transcript
#' 
#' @docType methods
#' @rdname pivot_transcript-methods
#' 
#' @return A `tbl` object
setMethod("pivot_transcript",
					"spec_tbl_df",
					.pivot_transcript)

#' pivot_transcript
#' @inheritParams pivot_transcript
#' 
#' @docType methods
#' @rdname pivot_transcript-methods
#' 
#' @return A `tbl` object
setMethod("pivot_transcript",
					"tbl_df",
					.pivot_transcript)

#' pivot_transcript
#' @inheritParams pivot_transcript
#' 
#' @docType methods
#' @rdname pivot_transcript-methods
#' 
#' @return A `tbl` object
setMethod("pivot_transcript",
					"tidybulk",
					.pivot_transcript)



#' Impute transcript abundance if missing from sample-transcript pairs
#'
#' \lifecycle{maturing}
#'
#' @description impute_abundance() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and returns a `tbl` with an edditional adjusted abundance column. This method uses scaled counts if present.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name impute_abundance
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .formula A formula with no response variable, representing the desired linear model where the first covariate is the factor of interest and the second covariate is the unwanted variation (of the kind ~ factor_of_intrest + batch)
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#'
#' @details This function imputes the abundance of missing sample-transcript pair using the median of the sample group defined by the formula
#'
#' @return A `tbl` non-sparse abundance
#'
#'
#'
#'
#' @examples
#'
#'
#' res =
#' 	impute_abundance(
#' 		tidybulk::counts_mini,
#' 	~ condition,
#' 	.sample = sample,
#' 	.transcript = transcript,
#' 	.abundance = count
#' )
#'
#'
#' @docType methods
#' @rdname impute_abundance-methods
#'
#' @export
#'
#'
setGeneric("impute_abundance", function(.data,
																				.formula,
																				.sample = NULL,
																				.transcript = NULL,
																				.abundance = NULL)
	standardGeneric("impute_abundance"))

# Set internal
.impute_abundance = 	function(.data,
															.formula,
															.sample = NULL,
															.transcript = NULL,
															.abundance = NULL)
{
	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
	.sample = col_names$.sample
	.transcript = col_names$.transcript
	.abundance = col_names$.abundance

	# Get scaled abundance if present, otherwise get abundance
	.abundance_scaled = NULL
	if(
		.data %>% get_tt_columns() %>% is.null %>% `!` &&
		".abundance_scaled" %in% (.data %>% get_tt_columns() %>% names) &&
		quo_name(.data %>% get_tt_columns() %$% .abundance_scaled) %in% (.data %>% colnames) &&
		quo_name(.data %>% get_tt_columns() %$% .abundance_scaled) != quo_name(.abundance)
	)
		.abundance_scaled = get_tt_columns(.data)$.abundance_scaled

	# Validate data frame
	validation(.data, !!.sample, !!.transcript, !!.abundance)

	.data_processed =
		fill_NA_using_formula(
			.data,
			.formula,
			.sample = !!.sample,
			.transcript = !!.transcript,
			.abundance = !!.abundance,
			.abundance_scaled = !!.abundance_scaled) %>%

		# Reattach internals
		reattach_internals(.data)

}

#' impute_abundance
#' @inheritParams impute_abundance
#' 
#' @docType methods
#' @rdname impute_abundance-methods
#' 
#' @return A `tbl` with imputed abundnce
setMethod("impute_abundance", "spec_tbl_df", .impute_abundance)

#' impute_abundance
#' @inheritParams impute_abundance
#' 
#' @docType methods
#' @rdname impute_abundance-methods
#' 
#' @return A `tbl` with imputed abundnce
setMethod("impute_abundance", "tbl_df", .impute_abundance)

#' impute_abundance
#' @inheritParams impute_abundance
#' 
#' @docType methods
#' @rdname impute_abundance-methods
#' 
#' @return A `tbl` with imputed abundnce
setMethod("impute_abundance", "tidybulk", .impute_abundance)

.impute_abundance_se = function(.data,
																.formula,
																.sample = NULL,
																.transcript = NULL,
																.abundance = NULL) {
	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	.data %>%

		# Convert to tidybulk
		tidybulk() %>%

		# Apply scale method
		impute_abundance(
			.formula = .formula,
			.sample = !!.sample,
			.transcript = !!.transcript,
			.abundance = !!.abundance
		) %>%

		# Convert to SummaizedExperiment
		tidybulk_to_SummarizedExperiment()

}

#' impute_abundance
#' @inheritParams impute_abundance
#' 
#' @docType methods
#' @rdname impute_abundance-methods
#' 
#' @return A `SummarizedExperiment` object
#'
setMethod("impute_abundance",
					"SummarizedExperiment",
					.impute_abundance_se)

#' impute_abundance
#' @inheritParams impute_abundance
#' 
#' @docType methods
#' @rdname impute_abundance-methods
#' 
#' @return A `SummarizedExperiment` object
#'
setMethod("impute_abundance",
					"RangedSummarizedExperiment",
					.impute_abundance_se)



