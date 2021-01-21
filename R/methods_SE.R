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

#' @importFrom magrittr multiply_by
#' @importFrom magrittr divide_by
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment colData
#' 
.scale_abundance_se = function(.data,
															 method = "TMM",
															 reference_sample = NULL) {
	
	
	# Check if package is installed, otherwise install
	if (find.package("edgeR", quiet = TRUE) %>% length %>% equals(0)) {
		message("Installing edgeR needed for analyses")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("edgeR", ask = FALSE)
	}
	
	# Check that reference sample exists
	if(!is.null(reference_sample) && !reference_sample %in% (.data %>% colnames))
		stop("tidybulk says: your reference sample is not among the samples in your data frame")
	
	

	.data_filtered = filter_if_abundant_were_identified(.data) 
	
	my_assay = assays(.data_filtered) %>% as.list() %>% .[1]
	my_counts_filtered = my_assay[[1]]
	library_size_filtered = my_counts_filtered %>% colSums()

	# Set column name for value scaled
	value_scaled = my_assay %>% names() %>% paste0(scaled_string)  
	
	# Get reference
	reference <-
		reference_sample %>%
		when(
			!is.null(.) ~ (.),
			
			# If not specified take most abundance sample
			library_size_filtered %>% 
				sort() %>% 
				tail(1) %>% 
				names()
		)
	
	# Calculate TMM
	nf <-
		edgeR::calcNormFactors(
			my_counts_filtered,
			refColumn = reference,
			method = method
		)
	
	# Calculate multiplier
	multiplier = 
		1 %>% 
		divide_by(library_size_filtered * nf) %>%
		
		# NOT HELPING - Put everything to the reference sample scale
		multiply_by(library_size_filtered[reference]) 
	
		# At the moment no because would be different from TIBBLE behaviour
		# %>%
		# 
		# # Make reference == 1
		# divide_by(.[reference])
		
	# Add to sample info
	colData(.data)$TMM = nf
	colData(.data)$multiplier = multiplier
	
	my_counts_scaled = 
		list(
			assays(.data) %>% 
				as.list() %>%
				.[[1]] %>%
				multiply_by(
					rep(multiplier, rep(nrow(.),length(multiplier)))
				)
			) %>%
		setNames(value_scaled)
	
	# Add the assay
	assays(.data) =  assays(.data) %>% c(my_counts_scaled)
	
	.data %>%
		
		# Add methods
		memorise_methods_used(c("edger", "tmm")) %>%
		
		# Attach column internals
		add_tt_columns(.abundance_scaled = !!(function(x, v)	enquo(v))(x,!!as.symbol(value_scaled))) 

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
																method ,
																of_samples = TRUE,
																log_transform = TRUE,
																...) {
	
	my_assay = 
		.data %>% 
		# Filter abundant if performed
		filter_if_abundant_were_identified() %>%
		assays() %>% 
		as.list() %>% 
		.[[get_assay_scaled_if_exists_SE(.data)]] 
	
	my_cluster_function  = 
		method %>%
			when(
				(.) == "kmeans" ~ get_clusters_kmeans_bulk_SE,
				(.) == "SNN" ~ get_clusters_SNN_bulk_SE,
				~ stop("tidybulk says: the only supported methods are \"kmeans\" or \"SNN\" ")
			)
	
	my_clusters = 
		my_cluster_function(
			my_assay,
			of_samples = of_samples,
			log_transform = log_transform,
			...
		) %>%
		as.character() %>%
		as.factor()
	
	my_cluster_column = paste("cluster", method, sep="_")
		
	.data %>%
		
		# Add clusters to metadata
		when(
			of_samples ~ {.x = (.); colData(.x)[,my_cluster_column] = my_clusters; .x},
			~ {.x = (.); rowData(.x)[,my_cluster_column] = my_clusters; .x}
		) %>%
		
		# Add bibliography
		when(
			method == "kmeans" ~ memorise_methods_used(., "stats"),
			method == "SNN" ~ memorise_methods_used(., "seurat"),
			~ stop("tidybulk says: the only supported methods are \"kmeans\" or \"SNN\" ")
		)
	
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
																 method,
																 .dims = 2,
																 top = 500,
																 of_samples = TRUE,
																 log_transform = TRUE,
																 scale = TRUE,
																 ...) {
	
	my_assay = 
		.data %>% 
		
		# Filter abundant if performed
		filter_if_abundant_were_identified() %>%
		
		assays() %>% 
		as.list() %>% 
		.[[get_assay_scaled_if_exists_SE(.data)]] %>%
		
		# Filter most variable genes
		keep_variable_transcripts_SE(top = top, log_transform = log_transform) %>%
		
		# Check if log transform is needed
		when(log_transform ~ log1p(.), ~ (.) ) 
	
	my_reduction_function  = 
		method %>%
		when(
			(.) == "MDS" ~ get_reduced_dimensions_MDS_bulk_SE,
			(.) == "PCA" ~ get_reduced_dimensions_PCA_bulk_SE,
			(.) == "tSNE" ~ get_reduced_dimensions_TSNE_bulk_SE,
			~ stop("tidybulk says: method must be either \"MDS\" or \"PCA\" or \"tSNE\"")
		)
	
	# Both dataframe and raw result object are returned
	reduced_dimensions = 
		my_reduction_function(
			my_assay,
			.dims = .dims,
			top = top,
			of_samples = of_samples,
			log_transform = log_transform,
			...
		)
	
	.data %>%
		
		# Add dimensions to metadata
		when(
			of_samples ~ {.x = (.); colData(.x) = colData(.x) %>% cbind(reduced_dimensions$result); .x},
			~ {.x = (.); rowData(.x) = rowData(.x) %>% cbind(reduced_dimensions$result); .x}
		) %>%
	
		# Add bibliography
		when(
			method == "MDS" ~ memorise_methods_used(., "limma"),
			method == "PCA" ~ memorise_methods_used(., "stats"),
			method == "tSNE" ~ memorise_methods_used(., "rtsne"),
			~ stop("tidybulk says: the only supported methods are \"kmeans\" or \"SNN\" ")
		) %>%
		
		# Attach edgeR for keep variable filtering
		memorise_methods_used(c("edger")) %>%
		
		# Add raw object
		attach_to_internals(reduced_dimensions$raw_result, method) %>%
		
		# Communicate the attribute added
		{
			message(sprintf("tidybulk says: to access the raw results do `attr(..., \"internals\")$%s`", method))
			(.)
		}

	
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
																 of_samples = TRUE,
																 dimension_1_column_rotated = NULL,
																 dimension_2_column_rotated = NULL) {
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
	
	# Sanity check of the angle selected
	if (rotation_degrees %>% between(-360, 360) %>% not())
		stop("tidybulk says: rotation_degrees must be between -360 and 360")
	
	# Return
	my_rotated_dimensions = 
		.data %>%
		
		# Select correct annotation
		when(
			of_samples ~ colData(.),
			~ rowData(.)
		) %>%
		
		# Select dimensions
		.[,c(quo_name(dimension_1_column), quo_name(dimension_2_column))] %>%
		as.matrix() %>%
		t() %>%
		rotation(rotation_degrees) %>%
		t() %>%
		as.data.frame() %>%
		setNames(c(
			quo_name(dimension_1_column_rotated),
			quo_name(dimension_2_column_rotated)
		))

	
	.data %>%
			
		# Add dimensions to metadata
		when(
			of_samples ~ {.x = (.); colData(.x) = colData(.x) %>% cbind(my_rotated_dimensions); .x},
			~ {.x = (.); rowData(.x) = rowData(.x) %>% cbind(my_rotated_dimensions); .x}
		) 
	
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
																 method,
																 of_samples = TRUE,
																 correlation_threshold = 0.9,
																 top = Inf,
																 log_transform = FALSE,
																 
																 Dim_a_column = NULL,
																 Dim_b_column = NULL) {

	
	Dim_a_column = enquo(Dim_a_column)
	Dim_b_column = enquo(Dim_b_column)
	
	redundant_elements = 
		method %>%
		when(
			. == "correlation" ~ {
				
				# Get counts
				my_assay = 
					.data %>% 
					
					# Filter abundant if performed
					filter_if_abundant_were_identified() %>%
					
					assays() %>% 
					as.list() %>% 
					.[[get_assay_scaled_if_exists_SE(.data)]] %>%
					
					# Filter most variable genes
					keep_variable_transcripts_SE(top = top, log_transform = log_transform) %>%
					
					# Check if log transform is needed
					when(log_transform ~ log1p(.), ~ (.) ) 
				
				# Get correlated elements
				remove_redundancy_elements_through_correlation_SE(
					my_assay,
					correlation_threshold = correlation_threshold,
					of_samples = of_samples
				)
			}	,
			. == "reduced_dimensions" ~ {
				
				# Get dimensions
				my_dims = 
					of_samples %>%
					when(
						of_samples ~ colData(.data)[,c(quo_name(Dim_a_column), quo_name(Dim_b_column))],
						~ rowData(.data)[,c(quo_name(Dim_a_column), quo_name(Dim_b_column))]
					)
				
				# Get correlated elements
				remove_redundancy_elements_though_reduced_dimensions_SE(
					my_dims
				)
			} ,
			~ stop(
				"tidybulk says: method must be either \"correlation\" for dropping correlated elements or \"reduced_dimension\" to drop the closest pair according to two dimensions (e.g., PCA)"
			)
		) 
	
		.data %>%
		
			# Condition on of_samples
			when(
				of_samples ~ (.)[,!colnames(.) %in% redundant_elements], 	
				~ (.)[-!rownames(.) %in% redundant_elements,]
			) %>%
				
			# Add bibliography
			when(
				method == "correlation" ~ memorise_methods_used(., "widyr"),
				method == "reduced_dimensions" ~ (.),
				~ stop("tidybulk says: method must be either \"correlation\" for dropping correlated elements or \"reduced_dimension\" to drop the closest pair according to two dimensions (e.g., PCA)")
			)
	
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
																log_transform = TRUE,
																...) {

	# Check if package is installed, otherwise install
	if (find.package("sva", quiet = TRUE) %>% length %>% equals(0)) {
		message("Installing sva - Combat needed for adjustment for unwanted variation")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("sva", ask = FALSE)
	}
	
	
	# Check that .formula includes at least two covariates
	if (parse_formula(.formula) %>% length %>% st(2))
		stop(
			"The .formula must contain two covariates, the first being the factor of interest, the second being the factor of unwanted variation"
		)
	
	# Check that .formula includes no more than two covariates at the moment
	if (parse_formula(.formula) %>% length %>% gt(3))
		warning("tidybulk says: Only the second covariate in the .formula is adjusted for, at the moment")
	
	# Create design matrix
	design =
		model.matrix(object = as.formula("~" %>% paste0(parse_formula(.formula)[1])),
								 # get first argument of the .formula
								 data = colData(.data)) %>%
		set_colnames(c("(Intercept)", parse_formula(.formula)[1]))
	
	
	
	my_batch = colData(.data)[, parse_formula(.formula)[2]]
	
	
	my_assay =
		.data %>%
		
		assays() %>%
		as.list() %>%
		.[[get_assay_scaled_if_exists_SE(.data)]] %>%
		
		# Check if log transform is needed
		when(log_transform ~ log1p(.), ~ (.))
	
	
	# Set column name for value scaled
	value_adjusted = get_assay_scaled_if_exists_SE(.data) %>% paste0(adjusted_string)
	
	# Calculate adjusted assay
	my_assay_adjusted =
		
		my_assay %>%
		
		# Add little noise to avoid all 0s for a covariate that would error combat code (not statistics that would be fine)
		`+` (rnorm(length(.), 0, 0.000001)) %>%
		
		# Run combat
		sva::ComBat(batch = my_batch,
								mod = design,
								prior.plots = FALSE,
								...)  %>%
		
		# Check if log transform needs to be reverted
		when(log_transform ~ expm1(.), ~ (.))
	
	
	# Add the assay
	my_assay_scaled = 
		list(my_assay_adjusted) %>% setNames(value_adjusted)
	
	assays(.data) =  assays(.data) %>% c(my_assay_scaled)
	
	# Return
	.data %>%
		
		# Add methods
		memorise_methods_used("sva") %>%
		
		# Attach column internals
		add_tt_columns(.abundance_adjusted = !!(function(x, v)
			enquo(v))(x, !!as.symbol(value_adjusted)))
	
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
																			reference = X_cibersort,
																			method = "cibersort",
																			prefix = "",
																			...) {

	my_assay =
		.data %>%
		
		assays() %>%
		as.list() %>%
		.[[get_assay_scaled_if_exists_SE(.data)]] 
	
	# Get the dots arguments
	dots_args = rlang::dots_list(...)
	
	my_proportions = 
		my_assay %>%
		
		# Run Cibersort or llsr through custom function, depending on method choice
		when(
			
			# Execute do.call because I have to deal with ...
			method %>% tolower %>% equals("cibersort") 	~ {
				
				# Check if package is installed, otherwise install
				if (find.package("class", quiet = TRUE) %>% length %>% equals(0)) {
					message("Installing class needed for Cibersort")
					install.packages("class", repos = "https://cloud.r-project.org", dependencies = c("Depends", "Imports"))
				}
				
				# Check if package is installed, otherwise install
				if (find.package("e1071", quiet = TRUE) %>% length %>% equals(0)) {
					message("Installing e1071 needed for Cibersort")
					install.packages("e1071", repos = "https://cloud.r-project.org", dependencies = c("Depends", "Imports"))
				}
				
				# Check if package is installed, otherwise install
				if (find.package("preprocessCore", quiet = TRUE) %>% length %>% equals(0)) {
					message("Installing preprocessCore needed for Cibersort")
					if (!requireNamespace("BiocManager", quietly = TRUE))
						install.packages("BiocManager", repos = "https://cloud.r-project.org")
					BiocManager::install("preprocessCore", ask = FALSE)
					
				}
				
				# Choose reference
				reference = reference %>% when(is.null(.) ~ X_cibersort, ~ .)
				
				# Validate reference
				validate_signature_SE(.data, reference, !!.transcript)
				
				do.call(my_CIBERSORT, list(Y = ., X = reference) %>% c(dots_args)) %$%
					proportions %>%
					as_tibble(rownames = quo_name(.sample)) %>%
					select(-`P-value`,-Correlation,-RMSE) 
			},
			
			# Don't need to execute do.call
			method %>% tolower %>% equals("llsr") ~ {
				
				# Choose reference
				reference = reference %>% when(is.null(.) ~ X_cibersort, ~ .)
				
				# Validate reference
				validate_signature_SE(.data, reference, !!.transcript)
				
				(.) %>%
					run_llsr(reference) %>%
					as_tibble(rownames = quo_name(.sample)) 
			},
			
			# Don't need to execute do.call
			method %>% tolower %>% equals("epic") ~ {
				
				# Choose reference
				reference = reference %>% when(is.null(.) ~ "BRef", ~ .)
				
				(.) %>%
					run_epic(reference) %>%
					as_tibble(rownames = quo_name(.sample)) 
			},
			
			# Other (hidden for the moment) methods using third party wrapper https://icbi-lab.github.io/immunedeconv	
			method %>% tolower %in% c("mcp_counter", "quantiseq", "xcell") ~ {
				
				# Check if package is installed, otherwise install
				if (find.package("immunedeconv", quiet = TRUE) %>% length %>% equals(0)) {
					message("Installing immunedeconv")
					devtools::install_github("icbi-lab/immunedeconv", upgrade = FALSE)
				}
				
				if(method == "xcell" & !"immunedeconv" %in% (.packages()))
					stop("tidybulk says: for xcell deconvolution you should have the package immunedeconv attached. Please execute library(immunedeconv)")
				
				(.) %>%
					immunedeconv::deconvolute(method %>% tolower, tumor = FALSE) %>%
					gather(!!.sample, .proportion, -cell_type) %>%
					spread(cell_type,  .proportion)
			},
			
			~ stop(
				"tidybulk says: please choose between cibersort, llsr and epic methods"
			)
		)	 %>%
		
		# Parse results and return
		setNames(c(
			"sample",
			(.) %>% select(-1) %>% colnames() %>% sprintf("%s%s", prefix, .)
			
		)) 
	
	# Att proportions
	colData(.data) = colData(.data) %>% cbind(
		my_proportions %>% 
			as_matrix(rownames = "sample") %>%
			.[rownames(colData(.data)),]
		)
	
	.data %>%
		
		# Attach attributes
		memorise_methods_used(tolower(method))
	
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
																					 .contrasts = NULL,
																					 method = "edgeR_quasi_likelihood",
																					 scaling_method = "TMM",
																					 omit_contrast_in_colnames = FALSE,
																					 prefix = "")
{

	# Clearly state what counts are used
	message("=====================================
tidybulk says: All testing methods use raw counts, irrespective of if scale_abundance 
or adjust_abundance have been calculated. Therefore, it is essential to add covariates 
such as batch effects (if applicable) in the formula.
=====================================")
	
	
	my_differential_abundance = 
		.data %>%
		
		# Filter abundant if performed
		filter_if_abundant_were_identified() %>%
		
		# Choose method
		when(
			
			# edgeR
			tolower(method) %in% c("edger_quasi_likelihood", "edger_likelihood_ratio", "edger_robust_likelihood_ratio") ~ 
				get_differential_transcript_abundance_bulk_SE(
					.,
					.formula,
					.contrasts = .contrasts,
					colData(.data),
					method = method,
					scaling_method = scaling_method,
					omit_contrast_in_colnames = omit_contrast_in_colnames,
					prefix = prefix
				),
			
			# Voom
			grepl("voom", method) ~ get_differential_transcript_abundance_bulk_voom_SE(
				.,
				.formula,
				.contrasts = .contrasts,
				colData(.data),
				method = method,
				scaling_method = scaling_method,
				omit_contrast_in_colnames = omit_contrast_in_colnames,
				prefix = prefix
			),
			
			# DESeq2
			tolower(method)=="deseq2" ~ get_differential_transcript_abundance_deseq2_SE(
				.,
				.formula,
				.contrasts = .contrasts,
				method = method,
				scaling_method = scaling_method,
				omit_contrast_in_colnames = omit_contrast_in_colnames,
				prefix = prefix
			),
			
			# Else error
			TRUE ~  stop("tidybulk says: the only methods supported at the moment are \"edgeR_quasi_likelihood\" (i.e., QLF), \"edgeR_likelihood_ratio\" (i.e., LRT), \"limma_voom\", \"limma_voom_sample_weights\", \"DESeq2\"")
		)

	# Add results
	rowData(.data) = rowData(.data) %>% cbind(
		my_differential_abundance$result %>%
			as_matrix(rownames = "transcript") %>%
			.[match(rownames(rowData(.data)), rownames(.)),]
	)
	
	
	.data %>%
		
		
		# Add bibliography
		when(
			tolower(method) ==  "edger_likelihood_ratio" ~ (.) %>% memorise_methods_used(c("edger", "edgeR_likelihood_ratio")),
			tolower(method) ==  "edger_quasi_likelihood" ~ (.) %>% memorise_methods_used(c("edger", "edgeR_quasi_likelihood")),
			tolower(method) ==  "edger_robust_likelihood_ratio" ~ (.) %>% memorise_methods_used(c("edger", "edger_robust_likelihood_ratio")),
			tolower(method) == "limma_voom" ~ (.) %>% memorise_methods_used("voom"),
			tolower(method) == "limma_voom_sample_weights" ~ (.) %>% memorise_methods_used("voom_sample_weights"),
			tolower(method) == "deseq2" ~ (.) %>% memorise_methods_used("DESeq2"),
			~ stop("tidybulk says: method must be either \"correlation\" for dropping correlated elements or \"reduced_dimension\" to drop the closest pair according to two dimensions (e.g., PCA)")
		) %>%
		
		attach_to_internals(my_differential_abundance$result_raw, method) %>%
		
		# Communicate the attribute added
		{
			message(
				sprintf("tidybulk says: to access the raw results (fitted GLM) do `attr(..., \"internals\")$%s`", method)
			)
			(.)
		}
	

	
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
														 top = 500,
														 log_transform = TRUE)
{

	
	variable_transcripts = 
		.data %>% 
		
		# Filter abundant if performed
		filter_if_abundant_were_identified() %>%
		
		assays() %>% 
		as.list() %>% 
		.[[get_assay_scaled_if_exists_SE(.data)]] %>%
		
		# Filter most variable genes
		keep_variable_transcripts_SE(top = top, log_transform = log_transform) %>%
		
		# Take gene names
		rowname()
	
	.data[variable_transcripts]
	
	
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
														 factor_of_interest = NULL,
														 minimum_counts = 10,
														 minimum_proportion = 0.7)
{
	
	
	# # Make col names
	# .sample = enquo(.sample)
	# .transcript = enquo(.transcript)
	# .abundance = enquo(.abundance)
	factor_of_interest = enquo(factor_of_interest)
	
browser()
	
	
	# if(!is.null(factor_of_interest))
	# 	factor_of_interest = enquo(factor_of_interest)
	# 
	# # Check factor_of_interest
	# if(
	# 	!is.null(factor_of_interest) && 
	# 	quo_name(factor_of_interest) %in% colnames(colData(.data)) %>% not()
	# )
	# 	stop(sprintf("tidybulk says: the column %s is not present in colData", quo_name(factor_of_interest)))
	# 
	# 
	# if (minimum_counts < 0)
	# 	stop("The parameter minimum_counts must be > 0")
	# if (minimum_proportion < 0 |	minimum_proportion > 1)
	# 	stop("The parameter minimum_proportion must be between 0 and 1")
	# 
	# # If column is present use this instead of doing more work
	# if(".abundant" %in% colnames(colData(.data))){
	# 	message("tidybulk says: the column .abundant already exists in colData. Nothing was done")
	# 	
	# 	# Return
	# 	return(.data)
	# }
	# 
	# 			
	# Check if factor_of_interest is continuous and exists
	
	#colData(.data) %>% as_tibble() %>% select(!!factor_of_interest)
	
	# string_factor_of_interest =
	# 
	# 	factor_of_interest %>%
	# 	when(
	# 		!is.null(factor_of_interest) &&
	# 			quo_is_symbol(factor_of_interest) &&
	# 			colData(.data)[, quo_name(factor_of_interest)] %>%
	# 			class %in%
	# 			c("numeric", "integer", "double") ~ {
	# 				message("tidybulk says: The factor of interest is continuous (e.g., integer,numeric, double). The data will be filtered without grouping.")
	# 				NULL
	# 			},
	# 		!is.null(factor_of_interest) &&
	# 			quo_is_symbol(factor_of_interest) ~
	# 			colData(.data)[, quo_name(factor_of_interest)],
	# 		~ NULL
	# 	)
	# 
	# # Get gene to exclude
	# gene_to_exclude =
	# 	.data %>%
	# 	
	# 	# Extract assay
	# 	assays() %>%
	# 	as.list() %>%
	# 	.[[1]] %>%
	# 	
	# 	# Call edgeR
	# 	edgeR::filterByExpr(
	# 		min.count = minimum_counts,
	# 		group = string_factor_of_interest,
	# 		min.prop = minimum_proportion
	# 	) %>%
	# 	not() %>%
	# 	which %>%
	# 	names 
	# 
	# rowData(.data)$.abundant = (rownames(rowData(.data)) %in% gene_to_exclude) %>% not()
	# 
	# # Return
	# .data
	
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
														 factor_of_interest = NULL,
														 minimum_counts = 10,
														 minimum_proportion = 0.7)
{
	
	
	factor_of_interest = enquo(factor_of_interest)
	
	.data = 
		.data %>%
		
		# Apply scale method
		identify_abundant(
			factor_of_interest = !!factor_of_interest,
			minimum_counts = minimum_counts,
			minimum_proportion = minimum_proportion
		) 
	
	.data[rowData(.data)$.abundant,]
	
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
	
	# DEPRECATED maybe
	#deprecate_warn("1.3.0", "tidybulk::fill_missing_abundance()")
	
	
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
																				.formula) {


	
	# Get scaled abundance if present, otherwise get abundance
	.abundance_scaled = NULL
	if(
		.data %>% get_tt_columns() %>% is.null %>% not() &&
		".abundance_scaled" %in% (.data %>% get_tt_columns() %>% names) &&
		quo_name(.data %>% get_tt_columns() %$% .abundance_scaled) %in% (.data %>% colnames) &&
		quo_name(.data %>% get_tt_columns() %$% .abundance_scaled) != quo_name(.abundance)
	)
		.abundance_scaled = get_tt_columns(.data)$.abundance_scaled
	
	# Validate data frame
	if(do_validate())  validation(.data, !!.sample, !!.transcript, !!.abundance)
	
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
																						 method = "cibersort",
																						 reference = X_cibersort,
																						 ...)
{
	if (find.package("broom", quiet = TRUE) %>% length %>% equals(0)) {
		message("Installing broom needed for analyses")
		install.packages("broom", repos = "https://cloud.r-project.org")
	}
	
	deconvoluted = 
		.data %>%
		
		# Deconvolution
		deconvolve_cellularity(
			method=method,
			prefix = sprintf("%s:", method),
			reference = reference,
			...
		) %>%
		colData()
	
	min_detected_proportion = 
		deconvoluted %>%
		select(starts_with(method)) %>%
		gather(cell_type, prop) %>%
		filter(prop > 0) %>%
		pull(prop) %>%
		min()
	
	
	# Check if test is univaiable or multivariable
	.formula %>%
		when(
			
			# Univariable
			format(.) %>%
				str_split("~") %>%
				.[[1]] %>%
				map_lgl( ~ grepl("\\. | \\.", .x)) %>%
				which	== 1 ~ {
					# Parse formula
					.my_formula =
						.formula %>%
						when(
							# If I have the dot, needed definitely for censored
							format(.) %>% grepl("\\.", .) %>% any ~ format(.) %>% str_replace("([-\\+\\*~ ]?)(\\.)", "\\1.proportion_0_corrected"),
							
							# If normal formula
							~ sprintf(".proportion_0_corrected%s", format(.))
						) %>%
						
						as.formula
					
					# Test
					univariable_differential_tissue_composition(deconvoluted,
																											method,
																											.my_formula,
																											min_detected_proportion) %>%
						
						# Attach attributes
						reattach_internals(.data) %>%
						
						# Add methods used
						when(
							grepl("Surv", .my_formula) ~ (.) %>% memorise_methods_used(c("survival", "boot")),
							~ (.) %>% memorise_methods_used("betareg")
						)
				},
			
			# Multivariable
			~ {
				# Parse formula
				covariates =
					deconvoluted %>%
					select(starts_with(method)) %>%
					colnames() %>%
					gsub(sprintf("%s:", method), "", .) %>%
					str_replace_all("[ \\(\\)]", "___")
				
				# Parse formula
				.my_formula =
					.formula %>%
					when(
						# If I have the dot, needed definitely for censored
						format(.) %>% grepl("\\.", .) %>% any ~
							format(.formula) %>%
							str_replace("([-\\+\\*~ ])(\\.)",
													sprintf(
														"\\1%s", paste(covariates, collapse = " + ")
													)),
						
						# If normal formula
						~ sprintf(".proportion_0_corrected%s", format(.))
					) %>%
					
					as.formula
				
				# Test
				multivariable_differential_tissue_composition(deconvoluted,
																											method,
																											.my_formula,
																											min_detected_proportion) %>%
					
					# Attach attributes
					reattach_internals(.data) %>%
					
					# Add methods used
					when(grepl("Surv", .my_formula) ~ (.) %>% memorise_methods_used(c("survival", "boot")),
							 ~ (.))
				
			}) %>%
		
		# Eliminate prefix
		mutate(.cell_type = str_remove(.cell_type, sprintf("%s:", method)))
	
	
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




