.tidybulk_se = function(.data,
												.sample,
												.transcript,
												.abundance,
												.abundance_scaled = NULL) {

  # Fix NOTEs
  . = NULL

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

	.as_tibble_optimised(.data) %>%

	# mutate_if(is.character, as.factor) %>%
	tidybulk(
		!!as.symbol(sample__$name),
		!!as.symbol(feature__$name),
		!!as.symbol(SummarizedExperiment::assays(.data)[1] %>%  names	),
		!!norm_col # scaled counts if any
	)



}

#' tidybulk
#'
#' @importFrom tibble as_tibble
#' @importFrom purrr reduce
#'
#' @export
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
#'
#' @export
#'
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
#' @importFrom stats na.omit
#'
.scale_abundance_se = function(.data,
                               .sample = NULL,
                               .transcript = NULL,
                               .abundance = NULL,
															 method = "TMM",
															 reference_sample = NULL,
															 .subset_for_scaling = NULL,
															 action = NULL,

															 # DEPRECATED
															 reference_selection_function = NULL) {


  # Fix NOTEs
  . = NULL

	# Check if package is installed, otherwise install
	if (find.package("edgeR", quiet = TRUE) %>% length %>% equals(0)) {
		message("Installing edgeR needed for analyses")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("edgeR", ask = FALSE)
	}

  # DEPRECATION OF reference function
  if (is_present(reference_selection_function) & !is.null(reference_selection_function)) {

    # Signal the deprecation to the user
    deprecate_warn("1.1.8", "tidybulk::scale_abundance(reference_selection_function = )", 
                   details = "The argument reference_selection_function is now deprecated please use reference_sample. By default the reference selection function is max()")

  }

	# Check that reference sample exists
	if(!is.null(reference_sample) && !reference_sample %in% (.data %>% colnames))
		stop("tidybulk says: your reference sample is not among the samples in your data frame")

  .subset_for_scaling = enquo(.subset_for_scaling)

	.data_filtered =
	  filter_if_abundant_were_identified(.data) %>%

	  # Filter based on user condition
	  when(
	    !quo_is_null(.subset_for_scaling) ~ filter_genes_on_condition(., !!.subset_for_scaling),
	    ~ (.)
	  ) %>%

	  # Check I have genes left
	  when(nrow(.) == 0 ~ stop("tidybulk says: there are 0 genes that passes the filters (.abundant and/or .subset_for_scaling). Please check your filtering or your data."), ~ (.))

	my_assay = assays(.data_filtered) %>% as.list() %>% .[1]

	# Drop genes with NAs, as edgeR::calcNormFactors does not accept them
	my_counts_filtered = my_assay[[1]] %>% na.omit()
	library_size_filtered = my_counts_filtered %>% colSums(na.rm  = TRUE)

	# If not enough genes, warning
	if(nrow(my_counts_filtered)<100) warning(warning_for_scaling_with_few_genes)

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

	# Communicate the reference if chosen by default
	if(is.null(reference_sample)) {
	  message(sprintf("tidybulk says: the sample with largest library size %s was chosen as reference for scaling", reference))
	}
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
		add_tt_columns(.abundance_scaled = !!(function(x, v) enquo(v))(x,!!as.symbol(value_scaled)))

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



#' @importFrom magrittr multiply_by
#' @importFrom magrittr divide_by
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment colData
#' @importFrom utils tail
#' @importFrom stats na.omit
#'
.quantile_normalise_abundance_se = function(.data,
                               .sample = NULL,
                               .transcript = NULL,
                               .abundance = NULL,
                               method = "limma_normalize_quantiles",
                               action = NULL) {


  # Fix NOTEs
  . = NULL

  .abundance = enquo(.abundance)

  # Set column name for value scaled

  # If no assay is specified take first
  my_assay = ifelse(
    quo_is_symbol(.abundance),
    quo_name(.abundance),
    .data |>
      assayNames() |>
      extract2(1)
  )

  # Set column name for value scaled
  value_scaled = my_assay %>% paste0(scaled_string)

  # Check if the matrix is empty and avoid error
  if(.data |> assay(my_assay) |> dim() |> min() == 0)
    .data_norm =
      .data |>
      assay(my_assay) |>
      list() |>
      setNames(value_scaled)

  else if(tolower(method) == "limma_normalize_quantiles"){

    # Check if package is installed, otherwise install
    if (find.package("limma", quiet = TRUE) %>% length %>% equals(0)) {
      message("tidybulk says: Installing limma needed for analyses")
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install("limma", ask = FALSE)
    }

    .data_norm <-
      .data %>%
      assay(my_assay) |>
      limma::normalizeQuantiles() |>
      list() |>
      setNames(value_scaled)

  }
  else if(tolower(method) == "preprocesscore_normalize_quantiles_use_target"){

    # Check if package is installed, otherwise install
    if (find.package("preprocessCore", quiet = TRUE) %>% length %>% equals(0)) {
      message("tidybulk says: Installing preprocessCore needed for analyses")
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install("preprocessCore", ask = FALSE)
    }

    .data_norm =
      .data |>
      assay(my_assay) |>
      as.matrix()

    .data_norm =
      .data_norm |>
      preprocessCore::normalize.quantiles.use.target(
        target = preprocessCore::normalize.quantiles.determine.target(.data_norm)
      )

    colnames(.data_norm) = .data |> assay(my_assay) |> colnames()
    rownames(.data_norm) = .data |> assay(my_assay) |> rownames()

    .data_norm =
      .data_norm |>
      list() |>
      setNames(value_scaled)

  } else stop("tidybulk says: the methods must be limma_normalize_quantiles or preprocesscore")

  # Add the assay
  assays(.data) =  assays(.data) %>% c(.data_norm)

  .data %>%

    # Add methods
    memorise_methods_used(c("quantile")) %>%

    # Attach column internals
    add_tt_columns(.abundance_scaled = !!(function(x, v)	enquo(v))(x,!!as.symbol(value_scaled)))

}

#' quantile_normalise_abundance
#' @inheritParams quantile_normalise_abundance
#'
#' @docType methods
#' @rdname quantile_normalise_abundance-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("quantile_normalise_abundance",
          "SummarizedExperiment",
          .quantile_normalise_abundance_se)

#' quantile_normalise_abundance
#' @inheritParams quantile_normalise_abundance
#'
#' @docType methods
#' @rdname quantile_normalise_abundance-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("quantile_normalise_abundance",
          "RangedSummarizedExperiment",
          .quantile_normalise_abundance_se)



.cluster_elements_se = function(.data,
																method ,
																of_samples = TRUE,
																transform = log1p,
																...) {

  # Fix NOTEs
  . = NULL

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
				(.) == "SNN" ~  stop("tidybulk says: Matrix package (v1.3-3) causes an ",
				"error with Seurat::FindNeighbors used in this method. We are trying to ",
				"solve this issue. At the moment this option in unaviable."), #get_clusters_SNN_bulk_SE,
				~ stop("tidybulk says: the only supported methods are \"kmeans\" or \"SNN\" ")
			)

	my_clusters =
		my_cluster_function(
			my_assay,
			of_samples = of_samples,
			transform = transform,
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
#' @importFrom rlang inform
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
                                 .abundance = NULL,

																 method,
																 .dims = 2,
																 top = 500,
																 of_samples = TRUE,
																 transform = log1p,
																 scale = TRUE,
																 ...) {

  # Fix NOTEs
  . = NULL

  .abundance = enquo(.abundance)

  if(.abundance |> quo_is_symbolic()) my_assay = quo_name(.abundance)
  else my_assay = get_assay_scaled_if_exists_SE(.data)

  # adjust top for the max number of features I have
  if(top > nrow(.data)){
    warning(sprintf(
      "tidybulk says: the \"top\" argument %s is higher than the number of features %s",
      top,
      nrow(.data)
    ))

    top = min(top, nrow(.data))
  }

	my_assay =
		.data %>%

		# Filter abundant if performed
		filter_if_abundant_were_identified() %>%

		assay(my_assay) %>%

		# Filter most variable genes
		keep_variable_transcripts_SE(top = top, transform = transform) %>%

		# Check if log transform is needed
		transform()

	my_reduction_function  =
		method %>%
		when(
			tolower(.) == tolower("MDS") ~ get_reduced_dimensions_MDS_bulk_SE,
			tolower(.) == tolower("PCA") ~ get_reduced_dimensions_PCA_bulk_SE,
			tolower(.) == tolower("tSNE") ~ get_reduced_dimensions_TSNE_bulk_SE,
			tolower(.) == tolower("UMAP") ~ get_reduced_dimensions_UMAP_bulk_SE,
			~ stop("tidybulk says: method must be either \"MDS\" or \"PCA\" or \"tSNE\", or \"UMAP\" ")
		)

	# Both dataframe and raw result object are returned
	reduced_dimensions =
		my_reduction_function(
			my_assay,
			.dims = .dims,
			top = top,
			of_samples = of_samples,
			transform = transform,
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
		  tolower(method) == tolower("MDS") ~ memorise_methods_used(., "limma"),
			tolower(method) == tolower("PCA") ~ memorise_methods_used(., "stats"),
			tolower(method) == tolower("tSNE") ~ memorise_methods_used(., "rtsne"),
			tolower(method) == tolower("UMAP") ~ memorise_methods_used(., "uwot"),
			~ stop("tidybulk says: method must be either \"MDS\" or \"PCA\" or \"tSNE\", or \"UMAP\" ")
		) %>%

		# Attach edgeR for keep variable filtering
		memorise_methods_used(c("edger")) %>%

		# Add raw object
		attach_to_internals(reduced_dimensions$raw_result, method) %>%

		# Communicate the attribute added
		{

		  rlang::inform(sprintf("tidybulk says: to access the raw results do `attr(..., \"internals\")$%s`", method), 
		                .frequency_id = sprintf("Access %s results", method),  .frequency = "always")

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

  # Fix NOTEs
  . = NULL

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


.remove_redundancy_se <- function(.data,
																 .element = NULL,
																 .feature = NULL,
																 .abundance = NULL,
																 method,
																 of_samples = TRUE,
																 correlation_threshold = 0.9,
																 top = Inf,
																 transform = identity,

																 Dim_a_column = NULL,
																 Dim_b_column = NULL,

																 # DEPRECATED
																 log_transform = NULL) {


  # Fix NOTEs
  . = NULL

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
					keep_variable_transcripts_SE(top = top, transform = transform) %>%

					# Check if log transform is needed
					transform()

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
				"tidybulk says: method must be either \"correlation\" for dropping ",
				"correlated elements or \"reduced_dimension\" to drop the closest pair according to two dimensions (e.g., PCA)"
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
				~ stop("tidybulk says: method must be either \"correlation\" for ",
				       "dropping correlated elements or \"reduced_dimension\" to drop ",
				       "the closest pair according to two dimensions (e.g., PCA)")
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
#' @importFrom rlang quo
#'
#' @docType methods
#' @rdname remove_redundancy-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("remove_redundancy",
					"RangedSummarizedExperiment",
					.remove_redundancy_se)


.adjust_abundance_se <- function(.data,
                                # DEPRECATED
                                .formula = NULL,
                                .factor_unwanted = NULL,
                                .factor_of_interest = NULL,
                                .abundance = NULL,
                                method = "combat_seq",
																...,

																# DEPRECATED
																transform = NULL,
																inverse_transform = NULL
																) {
  # Fix NOTEs
  . = NULL

  .abundance = enquo(.abundance)

	# Check if package is installed, otherwise install
	if (find.package("sva", quiet = TRUE) %>% length %>% equals(0)) {
		message("Installing sva - Combat needed for adjustment for unwanted variation")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("sva", ask = FALSE)
	}


  # DEPRECATION OF log_transform
  if (
    (is_present(transform) & !is.null(transform)) |
    is_present(inverse_transform) & !is.null(inverse_transform)
  ) {

    # Signal the deprecation to the user
    deprecate_warn("1.11.6", "tidybulk::test_differential_abundance(transform = )", 
      details = "The argument transform and inverse_transform is now deprecated, 
      please use method argument instead specifying \"combat\" or \"combat_seq\".")

  }

  # Set column name for value scaled
  value_adjusted = get_assay_scaled_if_exists_SE(.data) %>% paste0(adjusted_string)


  # DEPRECATION OF .formula
  if (is_present(.formula) & !is.null(.formula)) {

    # Signal the deprecation to the user
    deprecate_warn("1.11.6", "tidybulk::test_differential_abundance(.formula = )", 
      details = "The argument .formula is now deprecated, 
      please use factor_unwanted and factor_of_interest. 
      Using the formula, the first factor is of interest and the second is unwanted")

    # Check that .formula includes at least two covariates
    if (parse_formula(.formula) %>% length %>% st(2))
      stop(
        "The .formula must contain two covariates, 
        the first being the factor of interest, 
        the second being the factor of unwanted variation"
      )

    # Check that .formula includes no more than two covariates at the moment
    if (parse_formula(.formula) %>% length %>% gt(3))
      warning("tidybulk says: Only the second covariate in the .formula is adjusted for")


    .factor_of_interest = quo(!!as.symbol(parse_formula(.formula)[1]))
    .factor_unwanted = quo(!!as.symbol(parse_formula(.formula)[2]))

  } else {

    .factor_of_interest = enquo(.factor_of_interest)
    .factor_unwanted = enquo(.factor_unwanted)
  }

	# Create design matrix
  design =
    model.matrix(
      object = as.formula(sprintf("~ %s", colData(.data) |> as_tibble() |> 
                           select(!!.factor_of_interest) |> 
                             colnames() |> str_c(collapse = '+'))),
      # get first argument of the .formula
      data = colData(.data)
    )

	my_batch = colData(.data) |> as_tibble() |> select(!!.factor_unwanted)



	# If no assay is specified take first
	my_assay = ifelse(
	  quo_is_symbol(.abundance),
	  quo_name(.abundance),
	  get_assay_scaled_if_exists_SE(.data)
	)

	if(tolower(method) == "combat"){

	  my_assay_adjusted =
	    .data |>
	    assay(my_assay) |> # Check if log transform is needed
	   log1p() %>%
	    # Add little noise to avoid all 0s for a covariate that would error combat 
	    # code (not statistics that would be fine)
	    `+` (rnorm(length(.), 0, 0.000001))


	  for(i in colnames(my_batch)){
	    my_assay_adjusted =
	      my_assay_adjusted %>%
	      # Run combat
	      sva::ComBat(
	        batch = my_batch[,i] |> pull(1),
          mod = design,
          prior.plots = FALSE,
          ...
	       )
	  }

	  # Tranfrom back
	  my_assay_adjusted =
	    my_assay_adjusted |>
	      expm1() |>
	      (\(x) apply(x, 2, pmax, 0))()
	}
	else if (tolower(method) == "combat_seq"){

	  my_assay_adjusted =
	    .data %>%

	    assay(my_assay)

	  for(i in ncol(my_batch)){

	    my_assay_adjusted =
	      my_assay_adjusted |>
	      sva::ComBat_seq(batch = my_batch[,i] |> pull(1),
	                    covar_mod = design,
	                    full_mod=TRUE,
	                    ...)
	  }

	}
	else if(tolower(method) == "limma_remove_batch_effect") {

	  unwanted_covariate_matrix =
	    model.matrix(
	      object = as.formula(sprintf("~ 0 + %s", colData(.data) |> as_tibble() |> 
	                       select(!!.factor_unwanted) |> colnames() |> 
	                       str_c(collapse = '+'))),
	      # get first argument of the .formula
	      data = colData(.data)
	    )

	  my_assay_adjusted =
	    .data |>
	    assay(my_assay) |>
	    edgeR::cpm(log = T) |>
	    limma::removeBatchEffect(
	      design = design,
	      covariates = unwanted_covariate_matrix,
	      ...
	    ) |>
	    expm1() |>
	    apply(2, pmax, 0)

	} else {
	  stop("tidybulk says: the argument \"method\" must be \"combat_seq\", 
	       \"combat\", or \"limma_remove_batch_effect\"")
	}


	# Add the assay
	my_assay_scaled = list(my_assay_adjusted) %>% setNames(value_adjusted)

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


#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom GenomicRanges makeGRangesListFromDataFrame
#' @importFrom dplyr setdiff
.aggregate_duplicates_se <- function(.data,
																		.sample = NULL,
																		.transcript = NULL,
																		.abundance = NULL,
																		aggregation_function = sum,
																		keep_integer = TRUE) {
  # Fix NOTEs
  . = NULL

  # Make col names
  .transcript = enquo(.transcript)


  if(quo_is_null(.transcript)) stop("tidybulk says: using SummarizedExperiment 
     with aggregate_duplicates, you need to specify .transcript parameter. 
     It should be a feature-wise column (e.g. gene symbol) that you want to 
     collapse he features with (e.g. ensembl). It cannot be the representation 
     of rownames(SummarizedExperiment), as those are unique by definition, 
     and not part of rowData per-se.")

  if(!quo_name(.transcript) %in% colnames( .data %>% rowData()))
    stop("tidybulk says: the .transcript argument must be a feature-wise column 
          names. The feature-wise information can be found with rowData()")
  if(!is.null(.sample) | !is.null(.abundance))
    warning("tidybulk says: for SummarizedExperiment objects only the argument 
            .transcript (feature ID to collapse) is considered")

  collapse_function = function(x){ x %>% unique() %>% paste(collapse = "___")	}


  # Non standard column classes
  non_standard_columns =
    .data %>%
    rowData() %>%
    as_tibble() %>%
    select_if(select_non_standard_column_class) %>%
    colnames()

  # GRanges
  columns_to_collapse =
    .data %>%
    rowData() %>%
    colnames() %>%
    outersect(non_standard_columns) %>%
    setdiff(quo_name(.transcript)) %>%
    c(feature__$name)
    # when(
    #   !is.null(rownames(.data)) ~ c(., feature__$name),
    #   ~ (.)
    # )

  # Row data
  new_row_data =
    .data %>%
    rowData() %>%
    as_tibble(rownames = feature__$name) %>%
    group_by(!!as.symbol(quo_name(.transcript))) %>%
    summarise(
      across(columns_to_collapse, ~ .x %>% collapse_function()),
      across(non_standard_columns, ~ .x[1]),
      merged_transcripts = n()
    ) %>%

    arrange(!!as.symbol(feature__$name)) %>%
    as.data.frame() %>%
    {
      .x = (.)
      rownames(.x) = .x[,feature__$name]
      .x = .x %>% select(-feature__$name)
      .x
    }

  # If no duplicate exit
  if(!nrow(new_row_data)<nrow(rowData(.data))){
    message(sprintf("tidybulk says: your object does not have duplicates along 
        the %s column. The input dataset is returned.", quo_name(.transcript)))
    return(.data)
  }

  # Counts
  new_count_data =
    .data %>%
    assays() %>%
    as.list() %>%
    map(
      ~ {
        is_data_frame = .x %>% is("data.frame")
        if(is_data_frame) .x = .x %>% as.matrix()

        # Gove duplicated rownames
        rownames(.x) = rowData(.data)[,quo_name(.transcript)]

        # Combine
        if(rownames(.x) |> is.na() |> which() |> length() |> gt(0))
          stop(sprintf("tidybulk says: you have some %s that are NAs", 
                       quo_name(.transcript)))

        .x =  combineByRow(.x, aggregation_function)
        .x = .x[match(new_row_data[,quo_name(.transcript)], rownames(.x)),,drop=FALSE]
        rownames(.x) = rownames(new_row_data)

        if(is_data_frame) .x = .x %>% as.data.frame()
        .x
      }
    )

  if(!is.null(rowRanges(.data))){

    new_range_data = rowRanges(.data) %>% as_tibble()

    # If GRangesList & and .transcript is not there add .transcript
    if(is(rowRanges(.data), "CompressedGRangesList") & !quo_name(.transcript) %in% colnames(new_range_data)){

      new_range_data =
        new_range_data %>% left_join(
        rowData(.data)[,quo_name(.transcript),drop=FALSE] %>%
          as_tibble(rownames = feature__$name) ,
        by=c("group_name" = feature__$name)
      ) %>%
        select(-group_name, -group)
    }

    # Through warning if there are logicals of factor in the data frame
    # because they cannot be merged if they are not unique
    if (length(non_standard_columns)>0 & new_range_data %>% 
        pull(!!.transcript) %>% duplicated() %>% which() %>% length() %>% gt(0) ) {
      warning(paste(capture.output({
        cat(crayon::blue("tidybulk says: If duplicates exist from the following 
          columns, only the first instance was taken (lossy behaviour), 
          as aggregating those classes with concatenation is not possible.\n"))
        print(rowData(.data)[1,non_standard_columns,drop=FALSE])
      }), collapse = "\n"))
    }

    new_range_data = new_range_data %>%

      # I have to use this trick because rowRanges() and rowData() share @elementMetadata
      select(-one_of(colnames(new_row_data) %>% outersect(quo_name(.transcript)))) %>%
      suppressWarnings()


    #if(is(rr, "CompressedGRangesList") | nrow(new_row_data)<nrow(rowData(.data))) {
    new_range_data = makeGRangesListFromDataFrame(
        new_range_data,
        split.field = quo_name(.transcript),
        keep.extra.columns = TRUE
      )

      # Give back rownames
      new_range_data = new_range_data %>%  .[match(new_row_data[,quo_name(.transcript)], names(.))]
      #names(new_range_data) = rownames(new_row_data)
    #}
    # else if(is(rr, "GRanges")) new_range_data = makeGRangesFromDataFrame(new_range_data, keep.extra.columns = TRUE)
    # else stop("tidybulk says: riowRanges should be either GRanges or CompressedGRangesList. Or am I missing something?")

  }

  # Build the object
  .data_collapsed =
    SummarizedExperiment(
      assays = new_count_data,
      colData = colData(.data)
    )

  if(!is.null(rowRanges(.data))) rowRanges(.data_collapsed) = new_range_data

  rowData(.data_collapsed) = new_row_data

  .data_collapsed

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



#' @importFrom rlang quo_is_symbolic
.deconvolve_cellularity_se <- function(.data,
																			reference = X_cibersort,
																			method = "cibersort",
																			prefix = "",
																			...) {
  # Fix NOTEs
  . = NULL

  .transcript = enquo(.transcript)
  .sample = s_(.data)$symbol

	my_assay =
		.data %>%

		assays() %>%
		as.list() %>%
		.[[get_assay_scaled_if_exists_SE(.data)]] %>%

	  # Change row names
	  when(quo_is_symbolic(.transcript) ~ {
  	    .x = (.)
  	    rownames(.x) = .data %>% pivot_transcript() %>% pull(!!.transcript)
  	    .x
  	  },
  	  ~ (.)
	  )

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
					install.packages("class", repos = "https://cloud.r-project.org", 
					                 dependencies = c("Depends", "Imports"))
				}

				# Check if package is installed, otherwise install
				if (find.package("e1071", quiet = TRUE) %>% length %>% equals(0)) {
					message("Installing e1071 needed for Cibersort")
					install.packages("e1071", repos = "https://cloud.r-project.org", 
					                 dependencies = c("Depends", "Imports"))
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
				validate_signature_SE(., reference)

				do.call(my_CIBERSORT, list(Y = ., X = reference, QN=FALSE) %>% 
				  c(dots_args)) %$%
					proportions %>%
					as_tibble(rownames = quo_name(.sample)) %>%
					select(-`P-value`,-Correlation,-RMSE)
			},

			# Don't need to execute do.call
			method %>% tolower %>% equals("llsr") ~ {

				# Choose reference
				reference = reference %>% when(is.null(.) ~ X_cibersort, ~ .)

				# Validate reference
				validate_signature_SE(., reference)

				(.) %>%
					run_llsr(reference, ...) %>%
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
					stop("tidybulk says: for xcell, mcp_counter, 
					     or quantiseq deconvolution you should have the package 
					     immunedeconv attached. Please execute library(immunedeconv)")

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
			quo_name(.sample),
			(.) %>% select(-1) %>% colnames() %>% sprintf("%s%s", prefix, .)

		))

	# Att proportions
	colData(.data) = colData(.data) %>% cbind(
		my_proportions %>%
			as_matrix(rownames = .sample) %>%
		  .[match(rownames(colData(.data)), rownames(.)),]
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
#' @importFrom rlang inform
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


#' @importFrom rlang inform
.test_differential_abundance_se = function(.data,
																					 .formula,
																					 .abundance = NULL,
																					 contrasts = NULL,
																					 method = "edgeR_quasi_likelihood",
																					 test_above_log2_fold_change = NULL,
																					 scaling_method = "TMM",
																					 omit_contrast_in_colnames = FALSE,
																					 prefix = "",
																					 ...) {

  .abundance = enquo(.abundance)

  # Fix NOTEs
  . = NULL

  # DEPRECATION OF .constrasts
  if (is_present(.contrasts) & !is.null(.contrasts)) {

    # Signal the deprecation to the user
    deprecate_warn("1.7.4", "tidybulk::test_differential_abundance(.contrasts = )", 
     details = "The argument .contrasts is now deprecated please use contrasts (without the dot).")

    contrasts = .contrasts
  }

  # Clearly state what counts are used
  rlang::inform("=====================================
tidybulk says: All testing methods use raw counts, irrespective of if scale_abundance
or adjust_abundance have been calculated. Therefore, it is essential to add covariates
such as batch effects (if applicable) in the formula.
=====================================", 
    .frequency_id = "All testing methods use raw counts",  .frequency = "once")


	# Test test_above_log2_fold_change
	if(!is.null(test_above_log2_fold_change) && test_above_log2_fold_change < 0)
		stop("tidybulk says: test_above_log2_fold_change should be a positive real or NULL")

  # Filter abundant if performed
  .data = filter_if_abundant_were_identified(.data)

  if(tolower(method) %in% c("edger_quasi_likelihood", "edger_likelihood_ratio", 
                            "edger_robust_likelihood_ratio"))
  	my_differential_abundance =
      get_differential_transcript_abundance_bulk_SE(
        .data,
        .formula,
        .abundance = !!.abundance,
        .contrasts = contrasts,
        sample_annotation = colData(.data),
        method = method,
        test_above_log2_fold_change = test_above_log2_fold_change,
        scaling_method = scaling_method,
        omit_contrast_in_colnames = omit_contrast_in_colnames,
        prefix = prefix,
        ...
      )

  else if (grepl("voom", method))
  my_differential_abundance =
  get_differential_transcript_abundance_bulk_voom_SE(
    .data,
    .formula,
    .abundance = !!.abundance,
    .contrasts = contrasts,
    sample_annotation = colData(.data),
    method = method,
    test_above_log2_fold_change = test_above_log2_fold_change,
    scaling_method = scaling_method,
    omit_contrast_in_colnames = omit_contrast_in_colnames,
    prefix = prefix,
    ...
  )

  else if(tolower(method)=="deseq2")
  my_differential_abundance =
    get_differential_transcript_abundance_deseq2_SE(
      .data,
      .formula,
      .abundance = !!.abundance,
      .contrasts = contrasts,
      method = method,
      test_above_log2_fold_change = test_above_log2_fold_change,
      scaling_method = scaling_method,
      omit_contrast_in_colnames = omit_contrast_in_colnames,
      prefix = prefix,
      ...
    )


  else if(	tolower(method) %in% c("glmmseq_lme4", "glmmseq_glmmtmb"))
  my_differential_abundance =
    get_differential_transcript_abundance_glmmSeq_SE(
      .data,
    .formula,
    .abundance = !!.abundance,
    .contrasts = contrasts,
    sample_annotation = colData(.data),
    method = method,
    test_above_log2_fold_change = test_above_log2_fold_change,
    scaling_method = scaling_method,
    omit_contrast_in_colnames = omit_contrast_in_colnames,
    prefix = prefix,
    ...
  )
  else
    stop("tidybulk says: the only methods supported at the moment are 
         \"edgeR_quasi_likelihood\" (i.e., QLF), \"edgeR_likelihood_ratio\" 
         (i.e., LRT), \"limma_voom\", \"limma_voom_sample_weights\", \"DESeq2\", 
         \"glmmseq_lme4\", \"glmmseq_glmmTMB\"")


  statistics =
    my_differential_abundance$result %>%
    as_matrix(rownames = "transcript") %>%
    .[match(rownames(rowData(.data)), rownames(.)),,drop=FALSE]

  # If action is get just return the statistics
  if(action == "get") return(statistics)

	# Add results
	rowData(.data) = rowData(.data) %>% cbind(statistics)


	.data %>%


		# Add bibliography
		when(
			tolower(method) ==  "edger_likelihood_ratio" ~ (.) %>% 
			  memorise_methods_used(c("edger", "edgeR_likelihood_ratio")),
			tolower(method) ==  "edger_quasi_likelihood" ~ (.) %>% 
			  memorise_methods_used(c("edger", "edgeR_quasi_likelihood")),
			tolower(method) ==  "edger_robust_likelihood_ratio" ~ (.) %>% 
			  memorise_methods_used(c("edger", "edger_robust_likelihood_ratio")),
			tolower(method) == "limma_voom" ~ (.) %>% memorise_methods_used("voom"),
			tolower(method) == "limma_voom_sample_weights" ~ (.) %>% 
			  memorise_methods_used("voom_sample_weights"),
			tolower(method) == "deseq2" ~ (.) %>% memorise_methods_used("deseq2"),
			tolower(method) %in% c("glmmseq_lme4", "glmmseq_glmmtmb") ~ (.) %>% 
			  memorise_methods_used("glmmseq"),
			~ stop("tidybulk says: method not supported")
		) %>%

	    when(
			!is.null(test_above_log2_fold_change) ~ (.) %>% 
			  memorise_methods_used("treat"),
			~ (.)
		) %>%

		attach_to_internals(my_differential_abundance$result_raw, method) %>%

		# Communicate the attribute added
		{
		  rlang::inform(sprintf("tidybulk says: to access the raw results (fitted GLM) 
		      do `attr(..., \"internals\")$%s`", method), 
		      .frequency_id = sprintf("Access DE results %s", method),  
		      .frequency = "always")

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
														 transform = log1p)
{

  # Fix NOTEs
  . = NULL


	variable_transcripts =
		.data %>%

		# Filter abundant if performed
		filter_if_abundant_were_identified() %>%

		assays() %>%
		as.list() %>%
		.[[get_assay_scaled_if_exists_SE(.data)]] %>%

		# Filter most variable genes
		keep_variable_transcripts_SE(top = top, transform = transform) %>%

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
#' @importFrom purrr map_chr
#' @importFrom tidyr unite
#' @importFrom Matrix colSums
#'
#' @docType methods
#' @rdname keep_variable-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("keep_variable",
					"RangedSummarizedExperiment",
					.keep_variable_se)

.identify_abundant_se <- function(.data,
																 .sample = NULL,
																 .transcript = NULL,
																 .abundance = NULL,
																 factor_of_interest = NULL,
																 minimum_counts = 10,
																 minimum_proportion = 0.7) {
  # Fix NOTEs
  . = NULL


	factor_of_interest = enquo(factor_of_interest)
	.abundance = enquo(.abundance)

	# If character fail
	if(
	  !is.null(factor_of_interest) &&
	  !factor_of_interest |> quo_is_null() &&
	  !factor_of_interest |> quo_is_symbolic()
	) stop("tidybulk says: factor_of_interest must be symbolic 
	  (i.e. column name/s not surrounded by single or double quotes) 
	  and not a character.")


	# Check factor_of_interest
	if(
		!is.null(factor_of_interest) &&
		quo_is_symbolic(factor_of_interest) &&
		(quo_names(factor_of_interest) %in% colnames(colData(.data)) |> all() %>% not())
	)
		stop(sprintf("tidybulk says: the column %s is not present in colData", 
		             quo_names(factor_of_interest)))

	if (minimum_counts < 0)
		stop("The parameter minimum_counts must be > 0")

	if (minimum_proportion < 0 |	minimum_proportion > 1)
		stop("The parameter minimum_proportion must be between 0 and 1")

	# If column is present use this instead of doing more work
	if(".abundant" %in% colnames(colData(.data))){
		message("tidybulk says: the column .abundant already exists in colData. 
		        Nothing was done")

		# Return
		return(.data)
	}

	if(
	  !is.null(factor_of_interest) &&
	  ( enquo(factor_of_interest) |> quo_is_symbolic() | is.character(factor_of_interest) )
	){

	  # DEPRECATION OF symbolic factor_of_interest
	  # factor_of_interest must be a character now because we identified
	  # a edge case for which if the column name is the same as an existing function,
	  # such as time the column name would not be registered as such but would be
	  # registered as that function

	  # # Signal the deprecation to the user
	  # warning(
	  #   "The `factor_of_interest` argument of `test_differential_abundance() is changed as of tidybulk 1.11.5",
	  #   details = "The argument factor_of_interest must now be a character array. This because we identified a edge case for which if the column name is the same as an existing function, such as time the column name would not be registered as such but would be registered as that function"
	  # )

	  factor_of_interest = factor_of_interest |> enquo() |> quo_names()


	  # If is numeric ERROR
	  if(
	    colData(.data)[, factor_of_interest, drop=FALSE] |>
	          as_tibble() |>
	          map(~class(.x)) |>
	          unlist() %in% c("numeric", "integer", "double") |>
	          any()
	  )
	    stop("tidybulk says: The factor(s) of interest must not include continuous 
	         variables (e.g., integer,numeric, double).")

	  string_factor_of_interest =
	    colData(.data)[, factor_of_interest, drop=FALSE] |>
	    as_tibble() |>
	    unite("factor_of_interest") |>
	    select(factor_of_interest) |>
	    pull(1)


	} else {
	  string_factor_of_interest = NULL
	}



	# Check if package is installed, otherwise install
	if (find.package("edgeR", quiet = TRUE) %>% length %>% equals(0)) {
		message("Installing edgeR needed for analyses")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("edgeR", ask = FALSE)
	}

	# If no assay is specified take first
	my_assay = ifelse(
	  quo_is_symbol(.abundance),
	  quo_name(.abundance),
	  .data |>
	    assayNames() |>
	    extract2(1)
	)

	# Get gene to exclude
	gene_to_exclude =
	  .data |>
	  assay(my_assay) %>%

		# Call edgeR
		edgeR::filterByExpr(
			min.count = minimum_counts,
			group = string_factor_of_interest,
			min.prop = minimum_proportion,
			lib.size = Matrix::colSums(., na.rm=TRUE)
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

  # Fix NOTEs
  . = NULL


    factor_of_interest = factor_of_interest |> enquo()
    .abundance = enquo(.abundance)

	.data =
		.data %>%

		# Apply scale method
		identify_abundant(
			factor_of_interest = !!factor_of_interest,
			minimum_counts = minimum_counts,
			minimum_proportion = minimum_proportion,
			.abundance = !!.abundance
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



#' @importFrom lifecycle deprecate_warn
#' @importFrom stringr str_replace
#' @importFrom dplyr everything
#'
#'
#'
.test_gene_enrichment_SE <- function(
    .data,
		.formula,
		.sample = NULL,
		.entrez,
		.abundance = NULL,
		contrasts = NULL,
		methods = c("camera", "roast", "safe", "gage", "padog", "globaltest","ora"),
		gene_sets = c("h", "c1", "c2", "c3", "c4", "c5", "c6", "c7", 
		              "kegg_disease", "kegg_metabolism", "kegg_signaling"),
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
		deprecate_warn("1.3.2", "tidybulk::test_gene_enrichment(method = )", 
		      details = "The argument method is now deprecated please use methods")
		methods = method
	}

  # DEPRECATION OF .constrasts
  if (is_present(.contrasts) & !is.null(.contrasts)) {

    # Signal the deprecation to the user
    deprecate_warn("1.7.4", "tidybulk::test_differential_abundance(.contrasts = )", 
     details = "The argument .contrasts is now deprecated please use contrasts (without the dot).")

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
		stop("tidybulk says: There are duplicated .entrez IDs. 
		     Please use aggregate_duplicates(.transcript = entrez).")

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
	if (find.package("EGSEA", quiet = TRUE) %>% length %>% equals(0)) {
		stop("
		 EGSEA not installed. Please install it. EGSEA require manual installation 
		   for not overwelming the user in case it is not needed.
		 BiocManager::install(\"EGSEA\", ask = FALSE)
				 ")
	}
	if (!"EGSEA" %in% (.packages())) {
		stop("EGSEA package not loaded. Please run library(\"EGSEA\"). 
		     With this setup, EGSEA require manual loading, for technical reasons.")
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


    	idx =  buildIdx(entrezIDs = rownames(dge), species = species, 
    	                msigdb.gsets = msigdb.gsets,
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
    	message("tidybulk says: due to a bug in the call to KEGG database 
    	        (http://supportupgrade.bioconductor.org/p/122172/#122218), 
    	        the analysis for this database is run without report production.")

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

#' test_gene_enrichment
#' @inheritParams test_gene_enrichment
#'
#' @docType methods
#' @rdname test_gene_enrichment-methods
#'
#' @return A consistent object (to the input)
setMethod("test_gene_enrichment",
					"SummarizedExperiment",
					.test_gene_enrichment_SE)

#' test_gene_enrichment
#' @inheritParams test_gene_enrichment
#'
#' @docType methods
#' @rdname test_gene_enrichment-methods
#'
#' @return A consistent object (to the input)
setMethod("test_gene_enrichment",
					"RangedSummarizedExperiment",
					.test_gene_enrichment_SE)


# Set internal
.test_gene_overrepresentation_SE = 		function(.data,
																					 .entrez,
																					 .do_test,
																					 species,
																					 .sample = NULL,
																					 gene_sets = NULL,
																					 gene_set = NULL  # DEPRECATED
																					 )	{

	# Comply with CRAN NOTES
	. = NULL

	# DEPRECATION OF reference function
	if (is_present(gene_set) & !is.null(gene_set)) {

		# Signal the deprecation to the user
		deprecate_warn("1.3.1", "tidybulk::.test_gene_overrepresentation(gene_set = )", 
		 details = "The argument gene_set is now deprecated please use gene_sets.")
		gene_sets = gene_set
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
	if (.data %>% rowData() %>% as_tibble(rownames = f_(.data)$name) %>% 
	  mutate(my_do_test = !!.do_test) %>% pull(my_do_test) |> is("logical") %>% not())
		stop("tidybulk says: .do_test column must be logical (i.e., TRUE or FALSE)")

	# Check packages msigdbr
	# Check if package is installed, otherwise install
	if (find.package("msigdbr", quiet = TRUE) %>% length %>% equals(0)) {
		message("msigdbr not installed. Installing.")
		BiocManager::install("msigdbr", ask = FALSE)
	}

	# Check is correct species name
	if(species %in% msigdbr::msigdbr_species()$species_name %>% not())
		stop(sprintf("tidybulk says: wrong species name. MSigDB uses the latin species 
		    names (e.g., %s)", paste(msigdbr::msigdbr_species()$species_name, 
		                             collapse=", ")))

	# # Check if missing entrez
	# if(.data %>% filter(!!.entrez %>% is.na) %>% nrow() %>% gt(0) ){
	# 	warning("tidybulk says: there are .entrez that are NA. Those will be removed")
	# 	.data = .data %>%	filter(!!.entrez %>% is.na %>% not())
	# }

	.data %>%
		pivot_transcript(!!.entrez) %>%
		filter(!!.do_test) %>%
		distinct(!!.entrez) %>%
		pull(!!.entrez) %>%
		entrez_over_to_gsea(species, gene_collections = gene_sets) %>%

	  # Add methods used
	  memorise_methods_used(c("clusterProfiler", "msigdbr", "msigdb"), 
	                        object_containing_methods = .data)


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
																gene_sets = NULL,
																gene_set = NULL  # DEPRECATED
																)	{

	# Comply with CRAN NOTES
	. = NULL

	# DEPRECATION OF reference function
	if (is_present(gene_set) & !is.null(gene_set)) {

		# Signal the deprecation to the user
		deprecate_warn("1.3.1", "tidybulk::test_gene_rank(gene_set = )", 
		  details = "The argument gene_set is now deprecated please use gene_sets.")
		gene_sets = gene_set

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
		stop(sprintf("tidybulk says: wrong species name. MSigDB uses the 
		  latin species names (e.g., %s)", 
		  paste(msigdbr::msigdbr_species()$species_name, collapse=", ")))

	.data %>%
		pivot_transcript() %>%
		arrange(desc(!!.arrange_desc)) %>%
		select(!!.entrez, !!.arrange_desc) %>%
		deframe() %>%
		entrez_rank_to_gsea(species, gene_collections = gene_sets)%>%

	  # Add methods used. It is here and not in functions because I need the original .data
	  memorise_methods_used(c("clusterProfiler", "enrichplot"), 
	                        object_containing_methods = .data) %>%
	  when(
	    gene_sets %>% is("character") ~ (.) %>% memorise_methods_used("msigdbr"),
	    ~ (.)
	  )


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
#' @importFrom stringr str_replace
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
		tibble::as_tibble(rownames=sample__$name)




}

#' pivot_sample
#' @inheritParams pivot_sample
#'
#' @docType methods
#' @rdname pivot_sample-methods
#'
#' @return A consistent object (to the input)
setMethod("pivot_sample",
					"SummarizedExperiment",
					.pivot_sample)

#' pivot_sample
#' @inheritParams pivot_sample
#'
#' @docType methods
#' @rdname pivot_sample-methods
#'
#' @importFrom stringr str_replace
#'
#'
#' @return A consistent object (to the input)
setMethod("pivot_sample",
					"RangedSummarizedExperiment",
					.pivot_sample)



# Set internal
.pivot_transcript = 		function(.data,
															 .transcript = NULL)	{

  # Fix NOTEs
  . = NULL

	range_info <-
		get_special_datasets(.data) %>%
		reduce(left_join, by=feature__$name)

	gene_info <-
		rowData(.data) %>%

		# If reserved column names are present add .x
		setNames(
			colnames(.) %>%
				str_replace("^feature$", "feature.x")
		) %>%

		# Convert to tibble
		tibble::as_tibble(rownames=feature__$name)

	gene_info %>%
		when(
			nrow(range_info) > 0 ~ (.) %>% left_join(range_info, by=feature__$name),
			~ (.)
		)
}

#' pivot_transcript
#' @inheritParams pivot_transcript
#'
#' @docType methods
#' @rdname pivot_transcript-methods
#'
#' @return A consistent object (to the input)
setMethod("pivot_transcript",
					"SummarizedExperiment",
					.pivot_transcript)

#' pivot_transcript
#' @inheritParams pivot_transcript
#'
#' @docType methods
#' @rdname pivot_transcript-methods
#'
#' @return A consistent object (to the input)
setMethod("pivot_transcript",
					"RangedSummarizedExperiment",
					.pivot_transcript)


.impute_missing_abundance_se <- function(.data,
																				.formula,
																				.sample = NULL,
																				.transcript = NULL,
																				.abundance  = NULL,
																				suffix = "",
																				force_scaling = FALSE) {
  # Fix NOTEs
  . = NULL

  .abundance = enquo(.abundance)

  .assay_to_impute =
    .abundance %>%
    when(
      quo_is_symbolic(.) ~ assays(.data)[quo_names(.abundance)],
      ~ assays(.data)
    )


  # Split data by formula and impute
  imputed_dataframe =
    map2(

      # Capture assay names as we need to know if scaled is in the name
      as.list(.assay_to_impute), names(.assay_to_impute),
      ~ {

        # Pseudo-scale if not scaled
        if(!grepl("_scaled", .y) & force_scaling) {
            library_size = colSums(.x, na.rm = TRUE)
           .x = .x / library_size
        }
        else message(sprintf("tidybulk says: %s appears not to be scaled for 
          sequencing depth (missing _scaled suffix; if you think this column is 
          idependent of sequencing depth ignore this message), therefore the 
          imputation can produce non meaningful results if sequencing depth for 
          samples are highly variable. If you use force_scaling = TRUE library 
          size will be used for eliminatig some sequencig depth effect before 
          imputation", .y))

        # Log
        need_log = max(.x, na.rm=TRUE) > 50
        if(need_log) .x = log1p(.x)

        # Imputation
        .x = fill_NA_matrix_with_factor_colwise(
          .x,
          # I split according to the formula
          colData(.data)[,parse_formula(.formula)]
        )

        # Exp back
        if(need_log) .x = exp(.x)-1

        # Scale back if pseudoscaled
        if(!grepl("_scaled", .y) & force_scaling) .x = .x * library_size

        # Return
        .x
      }
    ) %>%

    # Add imputed to the name
    setNames(sprintf("%s%s", names(.assay_to_impute), suffix))

  .assays_name_to_port = names(assays(.data)) %>% setdiff(names(.assay_to_impute))

  assays(.data) =
    as.list(assays(.data))[.assays_name_to_port] %>%
    c(imputed_dataframe ) %>%

    # Add .imputed column
    c(list(.imputed =  which_NA_matrix(.assay_to_impute[[1]] ))) %>%

    # Make names unique
    setNames(names(.) %>% make.unique())

  .data %>%

    # Reattach internals
    reattach_internals(.data)

}



#' impute_missing_abundance
#' @inheritParams impute_missing_abundance
#'
#' @docType methods
#' @rdname impute_missing_abundance-methods
#'
#' @importFrom stringr str_replace
#'
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
#' @importFrom stringr str_replace
#'
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("impute_missing_abundance",
					"RangedSummarizedExperiment",
					.impute_missing_abundance_se)



.test_differential_cellularity_se <- function(.data,
																						 .formula,
																						 method = "cibersort",
																						 reference = X_cibersort,
																						 ...) {

  # Fix NOTEs
  . = NULL

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
							format(.) %>% grepl("\\.", .) %>% any ~ format(.) %>% 
							  str_replace("([-\\+\\*~ ]?)(\\.)", "\\1.proportion_0_corrected"),

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
							grepl("Surv", .my_formula) ~ (.) %>% 
							  memorise_methods_used(c("survival", "boot")),
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
					when(grepl("Surv", .my_formula) ~ (.) %>% 
					       memorise_methods_used(c("survival", "boot")),
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
#' @importFrom stringr str_replace
.test_stratification_cellularity_SE <- function(.data,
																							.formula,
																							.sample = NULL,
																							.transcript = NULL,
																							.abundance = NULL,
																							method = "cibersort",
																							reference = X_cibersort,
																							...) {

  # Fix NOTEs
  . = NULL

	# Validate formula
	if(.formula %>% format() %>% grepl(" \\.|\\. ", .) %>% not)
		stop("tidybulk says: in the formula a dot must be present in either these 
		forms \". ~\" or \"~ .\" with a white-space after or before respectively")

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
#' @return A consistent object (to the input) with additional columns for the 
#' statistics from the hypothesis test (e.g.,  log fold change, p-value and 
#' false discovery rate).
setMethod("test_stratification_cellularity",
					"SummarizedExperiment",
					.test_stratification_cellularity_SE)

#' test_stratification_cellularity
#' @inheritParams test_stratification_cellularity
#'
#' @docType methods
#' @rdname test_stratification_cellularity-methods
#'
#' @return A consistent object (to the input) with additional columns for the 
#' statistics from the hypothesis test (e.g.,  log fold change, p-value and 
#' false discovery rate).
setMethod("test_stratification_cellularity",
					"RangedSummarizedExperiment",
					.test_stratification_cellularity_SE)


#' get_bibliography
#' @inheritParams get_bibliography
#'
#' @docType methods
#' @rdname get_bibliography-methods
#'
#' @return A consistent object (to the input) with additional columns for the 
#' statistics from the hypothesis test (e.g.,  log fold change, p-value and 
#' false discovery rate).
setMethod("get_bibliography",
					"SummarizedExperiment",
					.get_bibliography)

#' get_bibliography
#' @inheritParams get_bibliography
#'
#' @docType methods
#' @rdname get_bibliography-methods
#'
#' @return A consistent object (to the input) with additional columns for the 
#' statistics from the hypothesis test (e.g.,  log fold change, p-value and 
#' false discovery rate).
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
.describe_transcript_SE <- function(.data,
															 .transcript = NULL) {

  # Fix NOTEs
  . = NULL

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
#' @return A consistent object (to the input) including additional 
#' columns for transcript symbol
setMethod("describe_transcript", "SummarizedExperiment", 
          .describe_transcript_SE)

#' describe_transcript
#' @inheritParams describe_transcript
#'
#' @docType methods
#' @rdname describe_transcript-methods
#'
#' @return A consistent object (to the input) including additional columns 
#' for transcript symbol
setMethod("describe_transcript", "RangedSummarizedExperiment", 
          .describe_transcript_SE)
