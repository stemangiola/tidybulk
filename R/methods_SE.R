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


.scale_abundance_se = function(.data,
															 .sample = NULL,
															 .transcript = NULL,
															 .abundance = NULL,
															 method = "TMM",
															 reference_selection_function = median,
															 action = "add") {
	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	
	.data %>%
		
		# Convert to tidybulk
		tidybulk() %>%
		
		# Apply scale method
		scale_abundance(
			!!.sample,
			!!.transcript,
			!!.abundance,
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
	# Get column namest
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



.test_differential_abundance_se = function(.data,
																					 .formula,
																					 .sample = NULL,
																					 .transcript = NULL,
																					 .abundance = NULL,
																					 .contrasts = NULL,
																					 method = "edgeR_quasi_likelihood",
																					 scaling_method = "TMM",
																					 omit_contrast_in_colnames = FALSE,
																					 prefix = "",
																					 
																					 action = "add",
																					 
																					 # DEPRECATED
																					 significance_threshold = 0.05,
																					 fill_missing_values = FALSE)
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

.identify_abundant_se = function(.data,
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
		identify_abundant(
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

#' identify_abundant
#' @inheritParams identify_abundant
#'
#' @docType methods
#' @rdname identify_abundant-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("identify_abundant",
					"SummarizedExperiment",
					.identify_abundant_se)

#' identify_abundant
#' @inheritParams identify_abundant
#'
#' @docType methods
#' @rdname identify_abundant-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("identify_abundant",
					"RangedSummarizedExperiment",
					.identify_abundant_se)




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




.fill_missing_abundance_se = function(.data,
																			.sample = NULL,
																			.transcript= NULL,
																			.abundance= NULL,
																			fill_with) {
	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	
	.data %>%
		
		# Convert to tidybulk
		tidybulk() %>%
		
		# Apply scale method
		fill_missing_abundance(
			.sample = !!.sample,
			.transcript = !!.transcript,
			.abundance = !!.abundance,
			fill_with = fill_with
		) %>%
		
		# Convert to SummaizedExperiment
		tidybulk_to_SummarizedExperiment()
	
}

#' fill_missing_abundance
#' @inheritParams fill_missing_abundance
#'
#' @docType methods
#' @rdname fill_missing_abundance-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("fill_missing_abundance",
					"SummarizedExperiment",
					.fill_missing_abundance_se)

#' fill_missing_abundance
#' @inheritParams fill_missing_abundance
#'
#' @docType methods
#' @rdname fill_missing_abundance-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("fill_missing_abundance",
					"RangedSummarizedExperiment",
					.fill_missing_abundance_se)



.impute_missing_abundance_se = function(.data,
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
		impute_missing_abundance(
			.formula = .formula,
			.sample = !!.sample,
			.transcript = !!.transcript,
			.abundance = !!.abundance
		) %>%
		
		# Convert to SummaizedExperiment
		tidybulk_to_SummarizedExperiment()
	
}

#' impute_missing_abundance
#' @inheritParams impute_missing_abundance
#'
#' @docType methods
#' @rdname impute_missing_abundance-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("impute_missing_abundance",
					"SummarizedExperiment",
					.impute_missing_abundance_se)

#' impute_missing_abundance
#' @inheritParams impute_missing_abundance
#'
#' @docType methods
#' @rdname impute_missing_abundance-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("impute_missing_abundance",
					"RangedSummarizedExperiment",
					.impute_missing_abundance_se)



.test_differential_cellularity_se = function(.data,
																						 .formula,
																						 .sample = NULL,
																						 .transcript = NULL,
																						 .abundance = NULL,
																						 method = "cibersort",
																						 reference = X_cibersort,
																						 significance_threshold = 0.05,
																						 ...)
{
	# Make col names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	
	.data %>%
		
		# Convert to tidybulk
		tidybulk() %>%
		
		# Apply scale method
		test_differential_cellularity_(
			.data,
			.formula = .formula,
			.sample = .sample,
			.transcript = .transcript,
			.abundance = .abundance,
			method = method,
			reference = reference,
			significance_threshold = significance_threshold,
			...
		)
	
}

#' test_differential_cellularity
#' @inheritParams test_differential_cellularity
#'
#' @docType methods
#' @rdname test_differential_cellularity-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod(
	"test_differential_cellularity",
	"SummarizedExperiment",
	.test_differential_cellularity_se
)

#' test_differential_cellularity
#' @inheritParams test_differential_cellularity
#'
#' @docType methods
#' @rdname test_differential_cellularity-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod(
	"test_differential_cellularity",
	"RangedSummarizedExperiment",
	.test_differential_cellularity_se
)




