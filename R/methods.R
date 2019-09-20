#' Creates a `tt` object from a `tbl``
#'
#' @description create_tt() creates a `tt` object from a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name create_ttBulk
#' @rdname create_ttBulk
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
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
#' \donttest{
#'
#'
#'
#' my_tt =
#'     ttBulk::counts %>%
#'     create_ttBulk(sample, transcript, counts)
#'
#' class(my_tt)
#'
#' }
#'
create_ttBulk <- function(.data,
													.sample,
													.transcript,
													.abundance) {
	UseMethod("create_ttBulk", .data)
}
#' @export
create_ttBulk.default <- function(.data,
																	.sample,
																	.transcript,
																	.abundance)
{
	print("This function cannot be applied to this object")
}
#' @export
create_ttBulk.tbl_df <- function(.data,
																 .sample,
																 .transcript,
																 .abundance)
{
	# Make col names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	create_tt_from_tibble_bulk(.data, !!.sample, !!.transcript,  !!.abundance)

}

#' Normalise the counts of transcripts/genes
#'
#' @description normalise_counts() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and normalises the data for the library size (e.g., with TMM algorithm, Robinson and Oshlack doi.org/10.1186/gb-2010-11-3-r25).
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name normalise_counts
#' @rdname normalise_counts
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
#' \donttest{
#'
#'
#'
#' counts %>%
#'     normalise_counts(sample, transcript, `count`)
#'
#'
#'}
#'
#' @export
normalise_counts <- function(.data,
														 .sample = NULL,
														 .transcript = NULL,
														 .abundance = NULL,
														 cpm_threshold = 0.5,
														 prop = 3 / 4,
														 method = "TMM",
														 reference_selection_function = median,
														 action = "add") {
	UseMethod("normalise_counts", .data)
}
#' @export
normalise_counts.default <-  function(.data,
																			.sample = NULL,
																			.transcript = NULL,
																			.abundance = NULL,
																			cpm_threshold = 0.5,
																			prop = 3 / 4,
																			method = "TMM",
																			reference_selection_function = median,
																			action = "add")
{
	print("This function cannot be applied to this object")
}
#' @export
normalise_counts.tbl_df = normalise_counts.ttBulk <-
	function(.data,
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

		if (action == "add")
			add_normalised_counts_bulk(.data,
																 !!.sample,
																 !!.transcript,
																 !!.abundance,
																 cpm_threshold = cpm_threshold,
																 prop = prop,
																 method = method,
																 reference_selection_function = reference_selection_function)
		else if (action == "get")
			get_normalised_counts_bulk(.data,
																 !!.sample,
																 !!.transcript,
																 !!.abundance,
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

#' Get clusters of elements (e.g., samples or transcripts)
#'
#' @description annotate_clusters() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and identify clusters in the data.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name annotate_clusters
#' @rdname annotate_clusters
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .element The name of the element column (normally samples).
#' @param .feature The name of the feature column (normally transcripts/genes)
#' @param .value The name of the column including the numerical value the clustering is based on (normally transcript abundance)
#'
#' @param number_of_clusters A integer indicating how many clusters you are seeking for the method k-means.
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
#' \donttest{
#'
#'
#'
#' counts %>%
#'     annotate_clusters(sample, transcript, count,	number_of_clusters = 2)
#'
#'     }
#'
#' @export
#'
annotate_clusters <- function(.data,
															.element = NULL,
															.feature = NULL,
															.value,
															number_of_clusters,
															method = "kmeans",
															of_samples = T,
															log_transform = T,
															action = "add", ...) {
	UseMethod("annotate_clusters", .data)
}

#' @export
annotate_clusters.default <-  function(.data,
																			 .element = NULL,
																			 .feature = NULL,
																			 .value,
																			 number_of_clusters,
																			 method = "kmeans",
																			 of_samples = T,
																			 log_transform = T,
																			 action = "add", ...)
{
	print("This function cannot be applied to this object")
}

#' @export
annotate_clusters.tbl_df = annotate_clusters.ttBulk <-  function(.data,
																																 .element = NULL,
																																 .feature = NULL,
																																 .value,
																																 number_of_clusters,
																																 method = "kmeans",
																																 of_samples = T,
																																 log_transform = T,
																																 action = "add", ...)
{
	# Make col names
	.value = enquo(.value)
	.element = enquo(.element)
	.feature = enquo(.feature)

	if(method == "kmeans"){
		if (action == "add")
			add_clusters_kmeans_bulk(
				.data,
				.value = !!.value,
				number_of_clusters = number_of_clusters,
				.element = !!.element,
				.feature = !!.feature,
				of_samples = of_samples,
				log_transform = log_transform,
				...
			)
		else if (action == "get")
			get_clusters_kmeans_bulk(
				.data,
				.value = !!.value,
				number_of_clusters = number_of_clusters,
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
		stop("the only supported method is \"kmeans\" ")

}


#' Dimension reduction of the transcript abundance data
#'
#' @description reduce_dimensions() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and calculates the reduced dimensional space of the transcript abundance.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name reduce_dimensions
#' @rdname reduce_dimensions
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .element The name of the element column (normally samples).
#' @param .feature The name of the feature column (normally transcripts/genes)
#' @param .value The name of the column including the numerical value the clustering is based on (normally transcript abundance)
#'
#' @param method A character string. The cluster algorithm to use, ay the moment k-means is the only algorithm included.
#' @param of_samples A boolean. In case the input is a ttBulk object, it indicates Whether the element column will be sample or transcript column
#' @param components A list of integer vectors corresponding to principal components of interest (e.g., list(1:2, 3:4, 5:6))

#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#' @param scale A boolean for method="PCA", this will be passed to the `prcomp` function. It is not included in the ... argument because although the default for `prcomp` if FALSE, it is advisable to set it as TRUE.
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#' @param ... Further parameters passed to the function prcomp if you choose method="PCA"
#'
#' @details This function reduces the dimensions of the transcript abundances.
#' It can use multi-dimensional scaling (MDS) of principal component analysis (PCA).
#'
#' @return A tbl object with additional columns for the reduced dimensions
#'
#'
#' @examples
#' \donttest{
#'
#'
#' library(GGally)
#'
#' counts.MDS =
#'     counts %>%
#'     reduce_dimensions(sample, transcript, count, method="MDS", components = 1:3)
#'
#' counts.MDS %>%
#'     select(contains("Dim"), sample, `Cell type`) %>%
#'     distinct() %>%
#'     GGally::ggpairs(columns = 1:3, ggplot2::aes(colour=`Cell type`))
#'
#' counts.PCA =
#'     counts.norm %>%
#'     reduce_dimensions(sample, transcript, count, method="PCA", components = 1:3)
#'
#'counts.PCA %>%
#'    select(contains("PC"), sample, `Cell type`) %>%
#'    distinct() %>%
#'    GGally::ggpairs(columns = 1:3, ggplot2::aes(colour=`Cell type`))
#'
#'}
#'
#' @export
#'
#'
reduce_dimensions <- function(.data,
															.value,
															method,
															components = 1:2,
															.element = NULL,
															.feature = NULL,
															of_samples = T,
															log_transform = T,
															scale = T,
															action = "add",
															...) {
	UseMethod("reduce_dimensions", .data)
}

#' @export
reduce_dimensions.default <-  function(.data,
																			 .value,
																			 method,
																			 components = 1:2,
																			 .element = NULL,
																			 .feature = NULL,
																			 of_samples = T,
																			 log_transform = T,
																			 scale = T,
																			 action = "add",
																			 ...)
{
	print("This function cannot be applied to this object")
}
#' @export
reduce_dimensions.tbl_df = reduce_dimensions.ttBulk <-
	function(.data,
					 .value,
					 method,
					 components = 1:2,
					 .element = NULL,
					 .feature = NULL,
					 of_samples = T,
					 log_transform = T,
					 scale = T,
					 action = "add",
					 ...)
	{
		# Make col names
		.value = enquo(.value)
		.element = enquo(.element)
		.feature = enquo(.feature)

		if (method == "MDS") {
			if (action == "add")
				add_reduced_dimensions_MDS_bulk(
					.data,
					.value = !!.value,
					components = components,
					.element = !!.element,
					.feature = !!.feature,
					of_samples = of_samples,
					log_transform = log_transform,
					...
				)
			else if (action == "get")
				get_reduced_dimensions_MDS_bulk(
					.data,
					.value = !!.value,
					components = components,
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
		else if (method == "PCA") {
			if (action == "add")
				add_reduced_dimensions_PCA_bulk(
					.data,
					.value = !!.value,
					components = components,
					.element = !!.element,
					.feature = !!.feature,
					of_samples = of_samples,
					log_transform = log_transform,
					scale = scale,
					...
				)
			else if (action == "get")
				get_reduced_dimensions_PCA_bulk(
					.data,
					.value = !!.value,
					components = components,
					.element = !!.element,
					.feature = !!.feature,
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



#' Rotate two dimensions (e.g., principal components) of an arbitrary angle
#'
#' @description rotate_dimensions() takes as imput a `tbl` formatted as | <DIMENSION 1> | <DIMENSION 2> | <...> | and calculates the rotated dimensional space of the transcript abundance.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name rotate_dimensions
#' @rdname rotate_dimensions
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
#' \donttest{
#'
#'
#'
#' counts.MDS.rotated =
#'     counts.MDS %>%
#'     rotate_dimensions(`Dim 1`, `Dim 2`, rotation_degrees = 45, .element = sample)
#'
#' counts.MDS.rotated %>%
#'     distinct(sample, `Dim 1`,`Dim 2`, `Cell type`) %>%
#'     ggplot(aes(x=`Dim 1`, y=`Dim 2`, color=`Cell type` )) +
#'     geom_point()
#'}
#'
#' @export
#'
rotate_dimensions <- function(.data,
															dimension_1_column,
															dimension_2_column,
															rotation_degrees,
															.element = NULL,
															of_samples = T,
															dimension_1_column_rotated = NULL,
															dimension_2_column_rotated = NULL,
															action = "add") {
	UseMethod("rotate_dimensions", .data)
}

#' @export
rotate_dimensions.default <-  function(.data,
																			 dimension_1_column,
																			 dimension_2_column,
																			 rotation_degrees,
																			 .element = NULL,
																			 of_samples = T,
																			 dimension_1_column_rotated = NULL,
																			 dimension_2_column_rotated = NULL,
																			 action = "add")
{
	print("This function cannot be applied to this object")
}

#' @export
rotate_dimensions.tbl_df = rotate_dimensions.ttBulk <-
	function(.data,
					 dimension_1_column,
					 dimension_2_column,
					 rotation_degrees,
					 .element = NULL,
					 of_samples = T,
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


#' Drop redundant elements (e.g., samples) for which feature (e.g., transcript/gene) aboundances are correlated
#'
#' @description drop_redundant() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | for correlation method or | <DIMENSION 1> | <DIMENSION 2> | <...> | for reduced_dimensions method, and returns a `tbl` with dropped elements (e.g., samples).
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name drop_redundant
#' @rdname drop_redundant
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .element The name of the element column (normally samples).
#' @param .feature The name of the feature column (normally transcripts/genes)
#' @param .value The name of the column including the numerical value the clustering is based on (normally transcript abundance)
#'
#' @param method A character string. The cluster algorithm to use, ay the moment k-means is the only algorithm included.
#' @param of_samples A boolean. In case the input is a ttBulk object, it indicates Whether the element column will be sample or transcript column
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#' @param correlation_threshold A real number between 0 and 1. For correlation based calculation.
#' @param Dim_a_column A character string. For reduced_dimension based calculation. The column of one principal component
#' @param Dim_b_column A character string. For reduced_dimension based calculation. The column of another principal component
#'
#' @details This function removes redundant elements from the original data set (e.g., samples or transcripts). For example, if we want to define cell-type specific signatures with low sample redundancy. This function returns a tibble with dropped recundant elements (e.g., samples). Two redundancy estimation approaches are supported: (i) removal of highly correlated clusters of elements (keeping a representative) with method="correlation"; (ii) removal of most proximal element pairs in a reduced dimensional space.
#'
#' @return A tbl object with with dropped recundant elements (e.g., samples).
#'
#' @examples
#' \donttest{
#'
#'
#'
#' counts %>%
#'     drop_redundant(
#' 	   .element = sample,
#' 	   .feature = transcript,
#' 	   	.value =  count,
#' 	   	method = "correlation"
#' 	   	)
#'
#' counts %>%
#'     drop_redundant(
#' 	   .element = sample,
#' 	   .feature = transcript,
#' 	   	.value = count,
#' 	   	method = "reduced_dimensions"
#' 	   	)
#'}
#'
#' @export
#'
#'
drop_redundant <- function(.data,
													 method,
													 .element = NULL,
													 of_samples = T,

													 .value,
													 .feature = NULL,
													 correlation_threshold = 0.9,
													 log_transform = F,

													 Dim_a_column,
													 Dim_b_column) {
	UseMethod("drop_redundant", .data)
}
#' @export
drop_redundant.default <-  function(.data,
																		method,
																		.element = NULL,
																		of_samples = T,

																		.value,
																		.feature = NULL,
																		correlation_threshold = 0.9,
																		log_transform = F,

																		Dim_a_column,
																		Dim_b_column)
{
	print("This function cannot be applied to this object")
}
#' @export
drop_redundant.tbl_df = drop_redundant.ttBulk <-  function(.data,
																													 method,
																													 .element = NULL,
																													 of_samples = T,

																													 .value,
																													 .feature = NULL,
																													 correlation_threshold = 0.9,
																													 log_transform = F,

																													 Dim_a_column = NULL,
																													 Dim_b_column = NULL)
{
	# Make col names
	.value = enquo(.value)
	.element = enquo(.element)
	.feature = enquo(.feature)

	Dim_a_column = enquo(Dim_a_column)
	Dim_b_column = enquo(Dim_b_column)

	if (method == "correlation")
		drop_redundant_elements_through_correlation(
			.data,
			.value = !!.value,
			.element = !!.element,
			.feature = !!.feature,
			correlation_threshold = correlation_threshold,
			of_samples = of_samples,
			log_transform = log_transform
		)
	else if (method == "reduced_dimensions")
		drop_redundant_elements_though_reduced_dimensions(
			.data,
			Dim_a_column = !!Dim_a_column,
			Dim_b_column = !!Dim_b_column,
			.element = !!.element,
			of_samples = of_samples
		)
	else
		stop(
			"method must be either \"correlation\" for dropping correlated elements or \"reduced_dimension\" to drop the closest pair according to two dimensions (e.g., PCA)"
		)

}



#' Adjust transcript abundance for unwanted variation
#'
#' @description adjust_abundance() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and returns a `tbl` with an edditional adjusted abundance column..
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name adjust_abundance
#' @rdname adjust_abundance
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
#' \donttest{
#'
#'
#'
#' adjust_abundance(
#'     ~ factor_of_interest + batch,
#'     sample,
#'     transcript,
#'     `count`
#' )
#'
#'}
#'
#' @export
#'
#'
adjust_abundance <- function(.data,
													.formula,
													.sample = NULL,
													.transcript = NULL,
													.abundance = NULL,
													log_transform = T,
													action = "add",
													...) {
	UseMethod("adjust_abundance", .data)
}
#' @export
adjust_abundance.default <-  function(.data,
																	 .formula,
																	 .sample = NULL,
																	 .transcript = NULL,
																	 .abundance = NULL,
																	 log_transform = T,
																	 action = "add",
																	 ...)
{
	print("This function cannot be applied to this object")
}
#' @export
adjust_abundance.tbl_df = adjust_abundance.ttBulk <-  function(.data,
																												 .formula,
																												 .sample = NULL,
																												 .transcript = NULL,
																												 .abundance = NULL,
																												 log_transform = T,
																												 action = "add",
																												 ...)
{
	# Make col names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

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


#' Aggregates multiple counts from the same samples (e.g., from isoforms), concatenates other character columns, and averages other numeric columns
#'
#' @description aggregate_duplicates() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and returns a `tbl` with aggregated transcripts that were duplicated.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name aggregate_duplicates
#' @rdname aggregate_duplicates
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
#' \donttest{
#'
#'
#'
#' counts %>%
#'     aggregate_duplicates(
#'     sample,
#'     transcript,
#'     `count`,
#'     aggregation_function = sum
#'     )
#'
#'}
#'
#' @export
#'
#'
aggregate_duplicates <- function(.data,
																 aggregation_function = sum,
																 .sample = NULL,
																 .transcript = NULL,
																 .abundance = NULL,
																 keep_integer = T) {
	UseMethod("aggregate_duplicates", .data)
}

#' @export
aggregate_duplicates.default <-  function(.data,
																					aggregation_function = sum,
																					.sample = NULL,
																					.transcript = NULL,
																					.abundance = NULL,
																					keep_integer = T)
{
	print("This function cannot be applied to this object")
}

#' @export
aggregate_duplicates.tbl_df = aggregate_duplicates.ttBulk <-
	function(.data,
					 aggregation_function = sum,
					 .sample = NULL,
					 .transcript = NULL,
					 .abundance = NULL,
					 keep_integer = T)  {
		# Make col names
		.sample = enquo(.sample)
		.transcript = enquo(.transcript)
		.abundance = enquo(.abundance)


		aggregate_duplicated_transcripts_bulk(
			.data,
			aggregation_function = aggregation_function,
			.sample = !!.sample,
			.transcript = !!.transcript,
			.abundance = !!.abundance,
			keep_integer = T
		)
	}


#' Get cell type proportions from samples
#'
#' @description annotate_cell_type() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and returns a `tbl` with the estimated cell type abundance for each sample
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name annotate_cell_type
#' @rdname annotate_cell_type
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#'
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
#' \donttest{
#'
#'
#'
#' 	counts %>%
#' 	    annotate_cell_type(sample, transcript, `count`)
#'
#'}
#'
#' @export
#'
annotate_cell_type <- function(.data,
															 .sample = NULL,
															 .transcript = NULL,
															 .abundance = NULL,
															 action = "add",
															 ...) {
	UseMethod("annotate_cell_type", .data)
}
#' @export
annotate_cell_type.default <-  function(.data,
																				.sample = NULL,
																				.transcript = NULL,
																				.abundance = NULL,
																				action = "add",
																				...)
{
	print("This function cannot be applied to this object")
}
#' @export
annotate_cell_type.tbl_df = annotate_cell_type.ttBulk <-
	function(.data,
					 .sample = NULL,
					 .transcript = NULL,
					 .abundance = NULL,
					 action = "add", ...)  {
		# Make col names
		.sample = enquo(.sample)
		.transcript = enquo(.transcript)
		.abundance = enquo(.abundance)

		if (action == "add")
			add_cell_type_proportions(
				.data,
				.sample = !!.sample,
				.transcript = !!.transcript,
				.abundance = !!.abundance,
				...
			)
		else if (action == "get")
			get_cell_type_proportions(
				.data,
				.sample = !!.sample,
				.transcript = !!.transcript,
				.abundance = !!.abundance,
				...
			)
		else
			stop(
				"action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
			)
	}



#' Add transcript symbol column from ensembl id
#'
#' @description annotate_symbol() takes as imput a `tbl` formatted as | <SAMPLE> | <ENSEMBL_ID> | <COUNT> | <...> | and returns a `tbl` with the the additional transcript symbol column
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name annotate_symbol
#' @rdname annotate_symbol
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
#' \donttest{
#'
#'
#'
#' 	counts_ensembl %>% annotate_symbol(ens)
#'
#'}
#'
#' @export
#'
#'
annotate_symbol <- function(.data,
														.ensembl,
														action = "add") {
	UseMethod("annotate_symbol", .data)
}
#' @export
annotate_symbol.default <-  function(.data,
																		 .ensembl,
																		 action = "add")
{
	print("This function cannot be applied to this object")
}
#' @export
annotate_symbol.tbl_df = annotate_symbol.ttBulk <-
	function(.data,
					 .ensembl,
					 action = "add")
	{
		# Make col names
		.ensembl = enquo(.ensembl)


		if (action == "add")
			add_symbol_from_ensembl(.data,!!.ensembl)

		else if (action == "get")
			get_symbol_from_ensembl(.data,!!.ensembl)

		else
			stop(
				"action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
			)

	}



#' Add differential transcription information to a tbl using edgeR.
#'
#' @description test_differential_transcription() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and returns a `tbl` with additional columns for the statistics from the hypothesis test.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name test_differential_transcription
#' @rdname test_differential_transcription
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .formula A formula with no response variable, representing the desired linear model
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#'
#' @param significance_threshold A real between 0 and 1 (usually 0.05).
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
#'\donttest{
#'
#'
#'
#' 	test_differential_transcription(
#' 	    ~ condition,
#' 	    sample,
#' 	    transcript,
#' 	    `count`
#' 	)
#'
#'}
#'
#' @export
#'
test_differential_transcription <- function(.data,
																								.formula,
																								.sample = NULL,
																								.transcript = NULL,
																								.abundance = NULL,
																								significance_threshold = 0.05,
																								action = "add") {
	UseMethod("test_differential_transcription", .data)
}
#' @export
test_differential_transcription.default <-  function(.data,
																												 .formula,
																												 .sample = NULL,
																												 .transcript = NULL,
																												 .abundance = NULL,
																												 significance_threshold = 0.05,
																												 action = "add")
{
	print("This function cannot be applied to this object")
}
#' @export
test_differential_transcription.tbl_df = test_differential_transcription.ttBulk <-
	function(.data,
					 .formula,
					 .sample = NULL,
					 .transcript = NULL,
					 .abundance = NULL,
					 significance_threshold = 0.05,
					 action = "add")
	{
		# Make col names
		.sample = enquo(.sample)
		.transcript = enquo(.transcript)
		.abundance = enquo(.abundance)

		if (action == "add")
			add_differential_transcript_abundance_bulk(
				.data,
				.formula,
				.sample = !!.sample,
				.transcript = !!.transcript,
				.abundance = !!.abundance,
				significance_threshold = significance_threshold
			)
		else if (action == "get")
			get_differential_transcript_abundance_bulk(
				.data,
				.formula,
				.sample = !!.sample,
				.transcript = !!.transcript,
				.abundance = !!.abundance,
				significance_threshold = significance_threshold
			)
		else
			stop(
				"action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
			)
	}
