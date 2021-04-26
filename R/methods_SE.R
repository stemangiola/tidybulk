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
#' @importFrom utils tail 
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
			scale=scale,
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
																 .element = NULL,
																 
																 of_samples = TRUE,
																 dimension_1_column_rotated = NULL,
																 dimension_2_column_rotated = NULL,
																 action = "add") {

	# Parse other colnames
	dimension_1_column = enquo(dimension_1_column)
	dimension_2_column = enquo(dimension_2_column)
	dimension_1_column_rotated = enquo(dimension_1_column_rotated)
	dimension_2_column_rotated = enquo(dimension_2_column_rotated)
	
	# Set default col names for rotated dimensions if not set
	if (quo_is_null(dimension_1_column_rotated))
		dimension_1_column_rotated = as.symbol(sprintf(
			"%s_rotated_%s",
			quo_name(dimension_1_column),
			rotation_degrees
		))
	if (quo_is_null(dimension_2_column_rotated))
		dimension_2_column_rotated = as.symbol(sprintf(
			"%s_rotated_%s",
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

	
	Dim_a_column = enquo(Dim_a_column)
	Dim_b_column = enquo(Dim_b_column)
	
	# Check if .data has more than one element
	if(
		(nrow(.data) <= 1 & of_samples == FALSE) |
		(ncol(.data) <= 1 & of_samples == TRUE) 
	)
		stop("tidybulk says: You must have more than one element (trancripts if of_samples == FALSE) to perform remove_redundancy")
	
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
								 data = colData(.data)) 

	# Maybe not needed and causing trouble if more columns that in the formula
	# %>%
	#set_colnames(c("(Intercept)", parse_formula(.formula)[1]))
	
	
	
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
				
				if(method %in% c("mcp_counter", "quantiseq", "xcell") & !"immunedeconv" %in% (.packages()))
					stop("tidybulk says: for xcell, mcp_counter, or quantiseq deconvolution you should have the package immunedeconv attached. Please execute library(immunedeconv)")
				
				(.) %>%
					deconvolute(method %>% tolower, tumor = FALSE) %>%
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
																					 test_above_log2_fold_change = NULL,
																					 scaling_method = "TMM",
																					 omit_contrast_in_colnames = FALSE,
																					 prefix = "",
																					 ...)
{

	# Clearly state what counts are used
	message("=====================================
tidybulk says: All testing methods use raw counts, irrespective of if scale_abundance 
or adjust_abundance have been calculated. Therefore, it is essential to add covariates 
such as batch effects (if applicable) in the formula.
=====================================")
	
	# Test test_above_log2_fold_change
	if(!is.null(test_above_log2_fold_change) && test_above_log2_fold_change < 0)
		stop("tidybulk says: test_above_log2_fold_change should be a positive real or NULL")
	
	
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
					test_above_log2_fold_change = test_above_log2_fold_change,
					scaling_method = scaling_method,
					omit_contrast_in_colnames = omit_contrast_in_colnames,
					prefix = prefix,
					...
				),
			
			# Voom
			grepl("voom", method) ~ get_differential_transcript_abundance_bulk_voom_SE(
				.,
				.formula,
				.contrasts = .contrasts,
				colData(.data),
				method = method,
				test_above_log2_fold_change = test_above_log2_fold_change,
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
	    
	    when(
			!is.null(test_above_log2_fold_change) ~ (.) %>% memorise_methods_used("treat"),
			~ (.)
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
		rownames()
	
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
																 .sample = NULL,
																 .transcript = NULL,
																 .abundance = NULL,
																 factor_of_interest = NULL,
																 minimum_counts = 10,
																 minimum_proportion = 0.7)
{
	
	
	factor_of_interest = enquo(factor_of_interest)
	
	
	# Check factor_of_interest
	if(
		!is.null(factor_of_interest) &&
		quo_is_symbol(factor_of_interest) &&
		(quo_name(factor_of_interest) %in% colnames(colData(.data)) %>% not())
	)
		stop(sprintf("tidybulk says: the column %s is not present in colData", quo_name(factor_of_interest)))
	
	if (minimum_counts < 0)
		stop("The parameter minimum_counts must be > 0")
	
	if (minimum_proportion < 0 |	minimum_proportion > 1)
		stop("The parameter minimum_proportion must be between 0 and 1")
	
	# If column is present use this instead of doing more work
	if(".abundant" %in% colnames(colData(.data))){
		message("tidybulk says: the column .abundant already exists in colData. Nothing was done")
		
		# Return
		return(.data)
	}
	
	
	string_factor_of_interest =
		
		factor_of_interest %>%
		when(
			quo_is_symbol(factor_of_interest) &&
				(
					colData(.data)[, quo_name(factor_of_interest)] %>%
						class %in% c("numeric", "integer", "double")) ~ 
				{
					message("tidybulk says: The factor of interest is continuous (e.g., integer,numeric, double). The data will be filtered without grouping.")
					NULL
				},
			quo_is_symbol(factor_of_interest) ~
				colData(.data)[, quo_name(factor_of_interest)],
			~ NULL
		)
	
	# Check if package is installed, otherwise install
	if (find.package("edgeR", quiet = TRUE) %>% length %>% equals(0)) {
		message("Installing edgeR needed for analyses")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("edgeR", ask = FALSE)
	}
	
	# Get gene to exclude
	gene_to_exclude =
		.data %>%
		
		# Extract assay
		assays() %>%
		as.list() %>%
		.[[1]] %>%
		
		# Call edgeR
		edgeR::filterByExpr(
			min.count = minimum_counts,
			group = string_factor_of_interest,
			min.prop = minimum_proportion
		) %>%
		not() %>%
		which %>%
		names
	
	rowData(.data)$.abundant = (rownames(rowData(.data)) %in% gene_to_exclude) %>% not()
	
	# Return
	.data
	
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




.test_gene_enrichment_SE = 		function(.data,
																			.formula,
																			.sample = NULL,
																			.entrez,
																			.abundance = NULL,
																			.contrasts = NULL,
																			method = c("camera" ,    "roast" ,     "safe",       "gage"  ,     "padog" ,     "globaltest",  "ora" ),
																			species,
																			cores = 10)	{

	.entrez = enquo(.entrez)

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
		stop("tidybulk says: You need at least two replicated for each condition for EGSEA to work")
	
	
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
		.contrasts %>%
		when(
			length(.) > 0 ~ limma::makeContrasts(contrasts = ., levels = design),
			~ NULL
			)
	
	# Check if package is installed, otherwise install
	if (find.package("EGSEA", quiet = TRUE) %>% length %>% equals(0)) {
		stop("
				 EGSEA not installed. Please install it. EGSEA require manual installation for not overwelming the user in case it is not needed.
				 BiocManager::install(\"EGSEA\", ask = FALSE)
				 ")
	}
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
	
	idx =  buildIdx(entrezIDs = rownames(dge), species = species)
	
	# Due to a bug in kegg, this data set is run without report
	# http://supportupgrade.bioconductor.org/p/122172/#122218
	message("tidybulk says: due to a bug in the call to KEGG database (http://supportupgrade.bioconductor.org/p/122172/#122218), the analysis for this database is run without report production.")
	
	res_kegg =
		dge %>%
		
		# Calculate weights
		limma::voom(design, plot = FALSE) %>%
		
		# Execute EGSEA
		egsea(
			contrasts = my_contrasts,
			gs.annots = idx %>% .["kegg"],
			baseGSEAs = method,
			sort.by = "med.rank",
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
		arrange(med.rank) %>%
		select(data_base, pathway, everything())
	
	res =
		dge %>%
		
		# Calculate weights
		limma::voom(design, plot = FALSE) %>%
		
		# Execute EGSEA
		egsea(
			contrasts = my_contrasts,
			gs.annots = idx[which(names(idx)!="kegg")],
			baseGSEAs = method,
			sort.by = "med.rank",
			num.threads = cores,
		)
	
	gsea_web_page = "https://www.gsea-msigdb.org/gsea/msigdb/cards/%s.html"
	
	res_formatted_all =
		res@results %>%
		map2_dfr(
			(.) %>% names,
			~ .x[[1]][[1]] %>%
				as_tibble(rownames = "pathway") %>%
				mutate(data_base = .y)
		) %>%
		arrange(med.rank) %>%
		
		# Add webpage
		mutate(web_page = sprintf(gsea_web_page, pathway)) %>%
		select(data_base, pathway, web_page, med.rank, everything()) 
	
	
	bind_rows(res_formatted_all, res_formatted_kegg)
	
	
}

#' test_gene_enrichment
#' @inheritParams test_gene_enrichment
#'
#' @docType methods
#' @rdname test_gene_enrichment-methods
#'
#' @return A `tbl` object
setMethod("test_gene_enrichment",
					"SummarizedExperiment",
					.test_gene_enrichment_SE)

#' test_gene_enrichment
#' @inheritParams test_gene_enrichment
#'
#' @docType methods
#' @rdname test_gene_enrichment-methods
#'
#' @return A `tbl` object
setMethod("test_gene_enrichment",
					"RangedSummarizedExperiment",
					.test_gene_enrichment_SE)


# Set internal
.test_gene_overrepresentation_SE = 		function(.data,
																					 .entrez,
																					 .do_test,
																					 species,
																					 .sample = NULL,
																					 gene_collections = NULL,
																					 gene_set = NULL  # DEPRECATED
																					 )	{

	# Comply with CRAN NOTES
	. = NULL
	
	# DEPRECATION OF reference function
	if (is_present(gene_set) & !is.null(gene_set)) {
		
		# Signal the deprecation to the user
		deprecate_warn("1.3.1", "tidybulk::.test_gene_overrepresentation(gene_set = )", details = "The argument gene_set is now deprecated please use gene_collections.")
		gene_collections = gene_set
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
	if (.data %>% rowData() %>% as_tibble() %>% distinct(!!.do_test) %>% sapply(class) %in% c("logical") %>% not() %>% any)
		stop("tidybulk says: .do_test column must be logical (i.e., TRUE or FALSE)")
	
	# Check packages msigdbr
	# Check if package is installed, otherwise install
	if (find.package("msigdbr", quiet = TRUE) %>% length %>% equals(0)) {
		message("msigdbr not installed. Installing.")
		BiocManager::install("msigdbr", ask = FALSE)
	}
	
	# Check is correct species name
	if(species %in% msigdbr::msigdbr_species()$species_name %>% not())
		stop(sprintf("tidybulk says: wrong species name. MSigDB uses the latin species names (e.g., %s)", paste(msigdbr::msigdbr_species()$species_name, collapse=", ")))
	
	.data %>%
		pivot_transcript() %>%
		filter(!!.do_test) %>%
		distinct(!!.entrez) %>%
		pull(!!.entrez) %>%
		entrez_over_to_gsea(species, gene_collections = gene_collections)
	
	
}

#' test_gene_overrepresentation
#' @inheritParams test_gene_overrepresentation
#'
#' @docType methods
#' @rdname test_gene_overrepresentation-methods
#'
#' @return A `SummarizedExperiment` object
setMethod("test_gene_overrepresentation",
					"SummarizedExperiment",
					.test_gene_overrepresentation_SE)

#' test_gene_overrepresentation
#' @inheritParams test_gene_overrepresentation
#'
#' @docType methods
#' @rdname test_gene_overrepresentation-methods
#'
#' @return A `RangedSummarizedExperiment` object
setMethod("test_gene_overrepresentation",
					"RangedSummarizedExperiment",
					.test_gene_overrepresentation_SE)


# Set internal
.test_gene_rank_SE = 		function(.data,
																.entrez,
																.arrange_desc,
																species,
																.sample = NULL,
																gene_collections = NULL,
																gene_set = NULL  # DEPRECATED
																)	{
	
	# Comply with CRAN NOTES
	. = NULL
	
	# DEPRECATION OF reference function
	if (is_present(gene_set) & !is.null(gene_set)) {
		
		# Signal the deprecation to the user
		deprecate_warn("1.3.1", "tidybulk::test_gene_rank(gene_set = )", details = "The argument gene_set is now deprecated please use gene_collections.")
		gene_collections = gene_set
		
	}
	
	# Get column names
	.arrange_desc = enquo(.arrange_desc)
	.entrez = enquo(.entrez)
	# 
	# expr <- rlang::quo_get_expr(.do_test)
	# env <- quo_get_env(x)
	# 
	
	# Check if entrez is set
	if(quo_is_missing(.entrez))
		stop("tidybulk says: the .entrez parameter appears to no be set")
	
	# Check packages msigdbr
	# Check if package is installed, otherwise install
	if (find.package("msigdbr", quiet = TRUE) %>% length %>% equals(0)) {
		message("msigdbr not installed. Installing.")
		BiocManager::install("msigdbr", ask = FALSE)
	}
	
	# Check is correct species name
	if(species %in% msigdbr::msigdbr_species()$species_name %>% not())
		stop(sprintf("tidybulk says: wrong species name. MSigDB uses the latin species names (e.g., %s)", paste(msigdbr::msigdbr_species()$species_name, collapse=", ")))
	
	.data %>%
		pivot_transcript() %>%
		arrange(desc(!!.arrange_desc)) %>%
		select(!!.entrez, !!.arrange_desc) %>%
		deframe() %>%
		entrez_rank_to_gsea(species, gene_collections = gene_collections)
	
	
}

#' test_gene_rank
#' @inheritParams test_gene_rank
#'
#' @docType methods
#' @rdname test_gene_rank-methods
#'
#' @return A `SummarizedExperiment` object
setMethod("test_gene_rank",
					"SummarizedExperiment",
					.test_gene_rank_SE)

#' test_gene_rank
#' @inheritParams test_gene_rank
#'
#' @docType methods
#' @rdname test_gene_rank-methods
#'
#' @return A `RangedSummarizedExperiment` object
setMethod("test_gene_rank",
					"RangedSummarizedExperiment",
					.test_gene_rank_SE)




# Set internal
.pivot_sample = 		function(.data,
													 .sample = NULL)	{
	
	colData(.data) %>%
		
		# If reserved column names are present add .x
		setNames(
			colnames(.) %>% 
				str_replace("^sample$", "sample.x")
		) %>%
		
		# Convert to tibble
		tibble::as_tibble(rownames="sample")
	

	
	
}

#' pivot_sample
#' @inheritParams pivot_sample
#'
#' @docType methods
#' @rdname pivot_sample-methods
#'
#' @return A `tbl` object
setMethod("pivot_sample",
					"SummarizedExperiment",
					.pivot_sample)

#' pivot_sample
#' @inheritParams pivot_sample
#'
#' @docType methods
#' @rdname pivot_sample-methods
#'
#' @return A `tbl` object
setMethod("pivot_sample",
					"RangedSummarizedExperiment",
					.pivot_sample)



# Set internal
.pivot_transcript = 		function(.data,
															 .transcript = NULL)	{

	range_info <-
		get_special_datasets(.data) %>%
		reduce(left_join, by="transcript")
	
	gene_info <-
		rowData(.data) %>%
		
		# If reserved column names are present add .x
		setNames(
			colnames(.) %>% 
				str_replace("^transcript$", "transcript.x")
		) %>%
		
		# Convert to tibble
		tibble::as_tibble(rownames="transcript")
	
	gene_info %>%
		when(
			nrow(range_info) > 0 ~ (.) %>% left_join(range_info, by="transcript"), 
			~ (.)
		) 
}

#' pivot_transcript
#' @inheritParams pivot_transcript
#'
#' @docType methods
#' @rdname pivot_transcript-methods
#'
#' @return A `tbl` object
setMethod("pivot_transcript",
					"SummarizedExperiment",
					.pivot_transcript)

#' pivot_transcript
#' @inheritParams pivot_transcript
#'
#' @docType methods
#' @rdname pivot_transcript-methods
#'
#' @return A `tbl` object
setMethod("pivot_transcript",
					"RangedSummarizedExperiment",
					.pivot_transcript)


.impute_missing_abundance_se = function(.data,
																				.formula) {


	
	
	
	col_formula =
		colData(.data) %>%
		as_tibble() %>%
		select(parse_formula(.formula)) %>%
		distinct() %>%
		select_if(function(x) is.character(x) | is.logical(x) | is.factor(x)) %>%
		colnames
	
	# Create NAs for missing sample/transcript pair
	assays(.data) = 
		assays(.data) %>%
		as.list() %>%
		map(~{
			.my_data = 
				.x %>%
				as.matrix() %>%
				as_tibble(rownames = "transcript") %>%
				gather(sample, abundance, -transcript) %>%
				
				# Attach annotation
				left_join(
					colData(.data) %>%
						as_tibble(rownames="sample") %>%
						select(sample, col_formula),
					by="sample"
				) 
			
			# Data used for filtering
			NA_data = 
				.my_data %>% 
				filter(abundance %>% is.na) %>% 
				select( transcript, col_formula) %>%
				distinct()
			
			# If no missing just return the same matrix
			if(nrow(NA_data) == 0) return(.x)
			
			.data_OK = 
				.my_data %>%
				anti_join(NA_data, by = c("transcript", col_formula))
			
			
			.data_FIXED = 
				.my_data %>%
				inner_join(NA_data, by = c("transcript", col_formula)) %>%
				
				# Group by covariate
				nest(cov_data = -c(col_formula, transcript)) %>%
				mutate(cov_data = map(cov_data, ~
																.x %>%
																mutate(abundance = 
																			 	case_when(
																			 		is.na(abundance) ~ median(abundance, na.rm=TRUE),	
																			 		TRUE ~ abundance
																			 	)
																) %>%
																
																# Through warning if group of size 1
																when(
																	nrow(.) %>% st(2) ~ warning("tidybulk says: According to your design matrix, u have sample groups of size < 2, so you your dataset could still be sparse."), 
																	~ (.)
																)
				)) %>%
				unnest(cov_data) 
			
			.data_OK %>%
				bind_rows(.data_FIXED) %>%
				select(-col_formula) %>%
				spread(sample, abundance) %>%
				as_matrix(rownames = transcript)
			
		})
	
	.data
	
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
		as_tibble(rownames = "sample") %>%
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
					univariable_differential_tissue_composition_SE(deconvoluted,
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
					as_tibble(rownames = "sample") %>%
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
				multivariable_differential_tissue_composition_SE(deconvoluted,
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

# Set internal
.test_stratification_cellularity_SE = 		function(.data,
																							.formula,
																							.sample = NULL,
																							.transcript = NULL,
																							.abundance = NULL,
																							method = "cibersort",
																							reference = X_cibersort,
																							...)
{
	
	# Validate formula
	if(.formula %>% format() %>% grepl(" \\.|\\. ", .) %>% not)
		stop("tidybulk says: in the formula a dot must be present in either these forms \". ~\" or \"~ .\" with a white-space after or before respectively")

	deconvoluted = 
		.data %>%
		
		# Deconvolution
		deconvolve_cellularity(
			method=method,
			prefix = sprintf("%s:", method),
			reference = reference,
			...
		) 
	
	# Check if test is univaiable or multivariable
	.formula %>%
		{
			# Parse formula
			.my_formula =
				format(.formula) %>% 
				str_replace("([~ ])(\\.)", "\\1.high_cellularity") %>%
				as.formula
			
			# Test
			univariable_differential_tissue_stratification_SE(deconvoluted,
																										 method,
																										 .my_formula) %>%
				
				# Attach attributes
				reattach_internals(.data) %>%
				
				# Add methods used
				memorise_methods_used(c("survival", "boot", "survminer"))
		} %>%
		
		# Eliminate prefix
		mutate(.cell_type = str_remove(.cell_type, sprintf("%s:", method)))
	
}

#' test_stratification_cellularity
#' @inheritParams test_stratification_cellularity
#'
#' @docType methods
#' @rdname test_stratification_cellularity-methods
#'
#' @return A `tbl` with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
setMethod("test_stratification_cellularity",
					"SummarizedExperiment",
					.test_stratification_cellularity_SE)

#' test_stratification_cellularity
#' @inheritParams test_stratification_cellularity
#'
#' @docType methods
#' @rdname test_stratification_cellularity-methods
#'
#' @return A `tbl` with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
setMethod("test_stratification_cellularity",
					"RangedSummarizedExperiment",
					.test_stratification_cellularity_SE)




#' get_bibliography
#' @inheritParams get_bibliography
#'
#' @docType methods
#' @rdname get_bibliography-methods
#'
#' @return A `tbl` with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
setMethod("get_bibliography",
					"SummarizedExperiment",
					.get_bibliography)

#' get_bibliography
#' @inheritParams get_bibliography
#'
#' @docType methods
#' @rdname get_bibliography-methods
#'
#' @return A `tbl` with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
setMethod("get_bibliography",
					"RangedSummarizedExperiment",
					.get_bibliography)

#' describe_transcript
#' 
#' @importFrom SummarizedExperiment rowData
#' @importFrom tibble enframe
#' 
#' @inheritParams describe_transcript
#'
#' @docType methods
#' @rdname describe_transcript-methods
#'
#' @return A `SummarizedExperiment` object
#'
#'
.describe_transcript_SE = function(.data,
															 .transcript = NULL) {
	
	# Check if package is installed, otherwise install
	if (find.package("org.Hs.eg.db", quiet = TRUE) %>% length %>% equals(0)) {
		message("Installing org.Hs.eg.db needed for differential transcript abundance analyses")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("org.Hs.eg.db", ask = FALSE)
	}
	
	# Check if package is installed, otherwise install
	if (find.package("org.Mm.eg.db", quiet = TRUE) %>% length %>% equals(0)) {
		message("Installing org.Mm.eg.db needed for differential transcript abundance analyses")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("org.Mm.eg.db", ask = FALSE)
	}
	
	# Check if package is installed, otherwise install
	if (find.package("AnnotationDbi", quiet = TRUE) %>% length %>% equals(0)) {
		message("Installing AnnotationDbi needed for differential transcript abundance analyses")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("AnnotationDbi", ask = FALSE)
	}
	
	.transcript = enquo(.transcript)

	# Transcript rownames by default
	my_transcripts = 
		.transcript %>%
		when(
			quo_is_null(.) ~ rownames(.data),
			~ rowData(.data)[,quo_name(.transcript)]
		)
		
	description_df = 
		# Human
		tryCatch(suppressMessages(AnnotationDbi::mapIds(
			org.Hs.eg.db::org.Hs.eg.db,
			keys = my_transcripts,  #ensembl_symbol_mapping$transcript %>% unique,
			column = "GENENAME",
			keytype = "SYMBOL",
			multiVals = "first"
		))  %>%
			.[!is.na(.)], error = function(x){}) %>%
		
		# Mouse
		c(
			tryCatch(suppressMessages(AnnotationDbi::mapIds(
				org.Mm.eg.db::org.Mm.eg.db,
				keys = my_transcripts,  #ensembl_symbol_mapping$transcript %>% unique,
				column = "GENENAME",
				keytype = "SYMBOL",
				multiVals = "first"
			)) %>% .[!is.na(.)], error = function(x){})
			
		) %>%
		
		# Parse
		unlist() %>%
		#unique() %>%
		enframe(name = "transcript", value = "description") %>%
		
		# Select just one per transcript
		distinct() %>%
		group_by(transcript) %>%
		slice(1) %>%
		ungroup()
	
	rowData(.data) = rowData(.data) %>% cbind(
		tibble(transcript = rownames(!!.data)) %>%
			left_join(description_df, by = "transcript") %>%
			select(description)
	)

	.data
}

#' describe_transcript
#' @inheritParams describe_transcript
#'
#' @docType methods
#' @rdname describe_transcript-methods
#'
#' @return A `tbl` object including additional columns for transcript symbol
setMethod("describe_transcript", "SummarizedExperiment", .describe_transcript_SE)

#' describe_transcript
#' @inheritParams describe_transcript
#'
#' @docType methods
#' @rdname describe_transcript-methods
#'
#' @return A `tbl` object including additional columns for transcript symbol
setMethod("describe_transcript", "RangedSummarizedExperiment", .describe_transcript_SE)
