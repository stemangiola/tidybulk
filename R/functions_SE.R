#' Get K-mean clusters to a tibble
#'
#' @keywords internal
#' @noRd
#'
#' 
#' 
#' @import tibble
#' @importFrom stats kmeans
#' @importFrom rlang :=
#'
#' @param .data A tibble
#' @param .abundance A column symbol with the value the clustering is based on (e.g., `count`)
#' @param .feature A column symbol. The column that is represents entities to cluster (i.e., normally samples)
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally genes)
#' @param of_samples A boolean
#' @param transform A function that will tranform the counts, by default it is log1p for RNA sequencing data, but for avoinding tranformation you can use identity
#' @param ... Further parameters passed to the function kmeans
#'
#' @return A tibble with additional columns
#'
#'
get_clusters_kmeans_bulk_SE <-
	function(.data,
					 of_samples = TRUE,
					 transform = log1p,
					 ...) {

		# Check if centers is in dots
		dots_args = rlang::dots_list(...)
		if ("centers" %in% names(dots_args) %>% not())
			stop("tidybulk says: for kmeans you need to provide the \"centers\" integer argument")

		.data %>%

			# Check if log transform is needed
			transform() %>%

			# Decide if of samples or transcripts
			when(
				of_samples ~ t(.),
				~ (.)
			) %>%

			# Wrap the do.call because of the centrers check
			{
				do.call(kmeans, list(x = (.), iter.max = 1000) %>% c(dots_args))
			}	 %$%
			cluster

	}

#' Get SNN shared nearest neighbour clusters to a tibble
#'
#' @keywords internal
#' @noRd
#'
#' 
#' 
#' @import tibble
#' @importFrom rlang :=
#' @importFrom utils install.packages
#'
#' @param .data A tibble
#' @param .abundance A column symbol with the value the clustering is based on (e.g., `count`)
#' @param .feature A column symbol. The column that is represents entities to cluster (i.e., normally samples)
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally genes)
#' @param of_samples A boolean
#' @param transform A function that will tranform the counts, by default it is log1p for RNA sequencing data, but for avoinding tranformation you can use identity
#' @param ... Further parameters passed to the function kmeans
#'
#' @return A tibble with additional columns
#'
get_clusters_SNN_bulk_SE <-
	function(.data,
					 of_samples = TRUE,
					 transform = log1p,
					 ...) {


		# Check if package is installed, otherwise install
		if (find.package("cluster", quiet = TRUE) %>% length %>% equals(0)) {
			message("Installing cluster")
			install.packages("cluster", repos = "https://cloud.r-project.org")
		}
		if (find.package("Seurat", quiet = TRUE) %>% length %>% equals(0)) {
			message("Installing Seurat")
			install.packages("Seurat", repos = "https://cloud.r-project.org")
		}
		if (find.package("KernSmooth", quiet = TRUE) %>% length %>% equals(0)) {
			message("Installing KernSmooth")
			install.packages("KernSmooth", repos = "https://cloud.r-project.org")
		}

		ndims = min(c(nrow(.data), ncol(.data), 30))-1

		.data %>%
			Seurat::CreateSeuratObject() %>%
			Seurat::ScaleData(display.progress = TRUE,
												num.cores = 4,
												do.par = TRUE) %>%
			Seurat::FindVariableFeatures(selection.method = "vst") %>%
			Seurat::RunPCA(npcs = ndims) %>%
			Seurat::FindNeighbors(dims = 1:ndims) %>%
			Seurat::FindClusters(method = "igraph", ...) %>%
			.[["seurat_clusters"]] %$%
			seurat_clusters

	}

#' Get dimensionality information to a tibble using MDS
#'
#' @keywords internal
#' @noRd
#'
#' 
#' 
#' @import tibble
#' @importFrom purrr map_dfr
#' @importFrom rlang :=
#' @importFrom stats setNames
#'
#' @param .data A tibble
#' @param .abundance A column symbol with the value the clustering is based on (e.g., `count`)
#' @param .dims A integer vector corresponding to principal components of interest (e.g., 1:6)
#' @param .feature A column symbol. The column that is represents entities to cluster (i.e., normally genes)
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally samples)
#' @param top An integer. How many top genes to select
#' @param of_samples A boolean
#' @param transform A function that will tranform the counts, by default it is log1p for RNA sequencing data, but for avoinding tranformation you can use identity
#' @param scale A boolean
#'
#' @return A tibble with additional columns
#'
#'
get_reduced_dimensions_MDS_bulk_SE <-
	function(.data,
					 .dims = 2,
					 top = 500,
					 of_samples = TRUE,
					 transform = log1p,
					 scale = NULL # This is only a dummy argument for making it compatibble with PCA
					) {
		# Comply with CRAN NOTES
		. = NULL

		# Get components from dims
		components = 1:.dims


		# Convert components to components list
		if((length(components) %% 2) != 0 ) components = components %>% append(components[1])
		components_list = split(components, ceiling(seq_along(components)/2))

		# Loop over components list and calculate MDS. (I have to make this process more elegant)
		mds_object =
			components_list %>%
			map(
				~ .data %>%
					limma::plotMDS(dim.plot = .x, plot = FALSE, top = top)
			)

		# Return
		list(
			raw_result = mds_object,
			result =
				map2_dfr(
					mds_object, components_list,
					~ {

						# Change of function from Bioconductor 3_13 of plotMDS
						my_rownames = .x %>% when(
							"distance.matrix.squared" %in% names(.x) ~ .x$distance.matrix.squared,
							~ .x$distance.matrix
						) %>%
							rownames()

						tibble(my_rownames, .x$x, .x$y) %>%
							rename(
								sample := my_rownames,
								!!as.symbol(.y[1]) := `.x$x`,
								!!as.symbol(.y[2]) := `.x$y`
							) %>%
							gather(Component, `Component value`,-sample)

					}


				)  %>%
				distinct() %>%
				spread(Component, `Component value`) %>%
				setNames(c((.) %>% select(1) %>% colnames(),
									 paste0("Dim", (.) %>% select(-1) %>% colnames())
				)) %>%
				select(-sample)
		)


	}





#' Get principal component information to a tibble using PCA
#'
#' @keywords internal
#' @noRd
#'
#' 
#' 
#' @import tibble
#' @importFrom rlang :=
#' @importFrom stats prcomp
#' @importFrom utils capture.output
#' @importFrom magrittr divide_by
#'
#' @param .data A tibble
#' @param .abundance A column symbol with the value the clustering is based on (e.g., `count`)
#' @param .dims A integer vector corresponding to principal components of interest (e.g., 1:6)
#' @param .feature A column symbol. The column that is represents entities to cluster (i.e., normally genes)
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally samples)
#' @param top An integer. How many top genes to select
#' @param of_samples A boolean
#' @param transform A function that will tranform the counts, by default it is log1p for RNA sequencing data, but for avoinding tranformation you can use identity
#' @param scale A boolean
#' @param ... Further parameters passed to the function prcomp
#'
#' @return A tibble with additional columns
#'
#'
get_reduced_dimensions_PCA_bulk_SE <-
	function(.data,
					 .dims = 2,
					 top = 500,
					 of_samples = TRUE,
					 transform = log1p,
					 scale = FALSE,
					 ...) {
		# Comply with CRAN NOTES
		. = NULL

		# Get components from dims
		components = 1:.dims

		prcomp_obj =
			.data %>%

			# check that there are non-NA genes for enough samples
			when(# First condition
				(.) %>% nrow == 0 ~ stop(
					"tidybulk says: In calculating PCA there is no gene that have non NA values is all samples"
				),

				# Second condition
				(.) %>% nrow < 100 ~ {
					warning(
						"
						tidybulk says: In PCA correlation there is < 100 genes that have non NA values is all samples.
The correlation calculation would not be reliable,
we suggest to partition the dataset for sample clusters.
						"
					)
					(.)
				},
				~ (.)) %>%

			t() %>%

			# Calculate principal components
			prcomp(scale = scale, ...)

		# Return
		list(
			raw_result = prcomp_obj,
			result =
				prcomp_obj %>%

				# Anonymous function - Prints fraction of variance
				# input: PCA object
				# output: PCA object
				{
					message("Fraction of variance explained by the selected principal components")

					(.) %$% sdev %>% pow(2) %>% # Eigen value
						divide_by(sum(.)) %>%
						`[` (components) %>%
						enframe() %>%
						select(-name) %>%
						rename(`Fraction of variance` = value) %>%
						mutate(PC = components) %>%
						capture.output() %>% paste0(collapse = "\n") %>% message()

					(.)

				} %$%

				# Parse the PCA results to a tibble
				x %>%
				as_tibble(rownames = "sample") %>%
				select(sprintf("PC%s", components))
		)


	}

#' Get principal component information to a tibble using tSNE
#'
#' @keywords internal
#' @noRd
#'
#' 
#' 
#' @import tibble
#' @importFrom rlang :=
#' @importFrom stats setNames
#' @importFrom utils install.packages
#'
#' @param .data A tibble
#' @param .abundance A column symbol with the value the clustering is based on (e.g., `count`)
#' @param .dims A integer vector corresponding to principal components of interest (e.g., 1:6)
#' @param .feature A column symbol. The column that is represents entities to cluster (i.e., normally genes)
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally samples)
#' @param top An integer. How many top genes to select
#' @param of_samples A boolean
#' @param transform A function that will tranform the counts, by default it is log1p for RNA sequencing data, but for avoinding tranformation you can use identity
#' @param scale A boolean
#' @param ... Further parameters passed to the function Rtsne
#'
#' @return A tibble with additional columns
#'
get_reduced_dimensions_TSNE_bulk_SE <-
	function(.data,
					 .dims = 2,
					 top = 500,
					 of_samples = TRUE,
					 transform = log1p,
					 scale = NULL, # This is only a dummy argument for making it compatibble with PCA
					 ...) {
		# Comply with CRAN NOTES
		. = NULL

		# To avoid dplyr complications


		# Evaluate ...
		arguments <- list(...)
		if (!"check_duplicates" %in% names(arguments))
			arguments = arguments %>% c(check_duplicates = FALSE)
		if (!"verbose" %in% names(arguments))
			arguments = arguments %>% c(verbose = TRUE)
		if (!"dims" %in% names(arguments))
			arguments = arguments %>% c(dims = .dims)


		# Check if package is installed, otherwise install
		if (find.package("Rtsne", quiet = TRUE) %>% length %>% equals(0)) {
			message("Installing Rtsne")
			install.packages("Rtsne", repos = "https://cloud.r-project.org")
		}

		# Set perprexity to not be too high
		if (!"perplexity" %in% names(arguments))
		  arguments = arguments %>% c(perplexity = ((
		    .data %>% ncol() %>% sum(-1)
		  ) / 3 / 2) %>% floor() %>% min(30))

		# If not enough samples stop
		if (arguments$perplexity <= 2)
			stop("tidybulk says: You don't have enough samples to run tSNE")

		# Calculate the most variable genes, from plotMDS Limma
		tsne_obj = do.call(Rtsne::Rtsne, c(list(t(.data)), arguments))



		list(
			raw_result = tsne_obj,
			result = tsne_obj %$%
				Y %>%
				as_tibble(.name_repair = "minimal") %>%
				setNames(c("tSNE1", "tSNE2")) %>%

				# add element name
				dplyr::mutate(sample = !!.data %>% colnames) %>%
				select(-sample)
		)

	}

#' Get UMAP
#'
#' @keywords internal
#'
#' 
#' 
#' @import tibble
#' @importFrom rlang :=
#' @importFrom stats setNames
#' @importFrom utils install.packages
#'
#' @param .data A tibble
#' @param .abundance A column symbol with the value the clustering is based on (e.g., `count`)
#' @param .dims A integer vector corresponding to principal components of interest (e.g., 1:6)
#' @param .feature A column symbol. The column that is represents entities to cluster (i.e., normally genes)
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally samples)
#' @param top An integer. How many top genes to select
#' @param of_samples A boolean
#' @param transform A function that will tranform the counts, by default it is log1p for RNA sequencing data, but for avoinding tranformation you can use identity
#' @param calculate_for_pca_dimensions An integer of length one. The number of PCA dimensions to based the UMAP calculatio on. If NULL all variable features are considered
#' @param ... Further parameters passed to the function uwot
#'
#' @return A tibble with additional columns
#'
get_reduced_dimensions_UMAP_bulk_SE <-
  function(.data,
           .dims = 2,
           top = 500,
           of_samples = TRUE,
           transform = log1p,
           scale = NULL, # This is only a dummy argument for making it compatibble with PCA
           calculate_for_pca_dimensions = 20,
           ...) {
    # Comply with CRAN NOTES
    . = NULL

    # To avoid dplyr complications

    # Evaluate ...
    arguments <- list(...)
    # if (!"check_duplicates" %in% names(arguments))
    #   arguments = arguments %>% c(check_duplicates = FALSE)
    if (!"dims" %in% names(arguments))
      arguments = arguments %>% c(n_components = .dims)
    if (!"init" %in% names(arguments))
      arguments = arguments %>% c(init = "spca")


    # Check if package is installed, otherwise install
    if (find.package("uwot", quiet = TRUE) %>% length %>% equals(0)) {
      message("tidybulk says: Installing uwot")
      install.packages("uwot", repos = "https://cloud.r-project.org")
    }


    # Calculate based on PCA
    if(!is.null(calculate_for_pca_dimensions))
      df_UMAP =
      .data %>%

      t() %>%

      # Calculate principal components
      prcomp(scale = scale) %$%

      # Parse the PCA results to a tibble
      x %>%
      .[,1:calculate_for_pca_dimensions]

    # Calculate based on all features
    else
      df_UMAP = .data

    umap_obj = do.call(uwot::tumap, c(list(df_UMAP), arguments))

    list(
      raw_result = umap_obj,
      result = umap_obj  %>%
        as_tibble(.name_repair = "minimal") %>%
        setNames(c("UMAP1", "UMAP2")) %>%

        # add element name
        dplyr::mutate(sample = !!.data %>% colnames) %>%
        select(-sample)
    )

  }

counts_scaled_exist_SE = function(.data){

	("tt_columns" %in% (.data %>%
												attr("internal") %>% names())) &&
	(
		.data %>%
		attr("internal") %$%
		tt_columns %>%
		names() %>%
		grep("scaled", .) %>%
		length() %>%
		equals(1)
	)
}

get_assay_scaled_if_exists_SE = function(.data){
	if(counts_scaled_exist_SE(.data))
		.data %>%
		attr("internal") %$%
		tt_columns %$%
		.abundance_scaled %>%
		quo_name()
	else
		.data %>%
		assays() %>%
		names() %>%
		head(1)
}

filter_if_abundant_were_identified = function(.data){
	.data %>%

		# Filter abundant if performed
		when(
			".abundant" %in% (rowData(.data) %>% colnames()) ~ .data[rowData(.data)[,".abundant"],],
			~ {
				warning("tidybulk says: highly abundant transcripts were not identified (i.e. identify_abundant()) or filtered (i.e., keep_abundant), therefore this operation will be performed on unfiltered data. In rare occasions this could be wanted. In standard whole-transcriptome workflows is generally unwanted.")
				(.)
			}
		)
}

#' Identify variable genes for dimensionality reduction
#'
#' @keywords internal
#' @noRd
#'
#' @param .data A tibble
#' @param .sample A character name of the sample column
#' @param .transcript A character name of the transcript/gene column
#' @param .abundance A character name of the read count column
#' @param top An integer. How many top genes to select
#' @param transform A function that will tranform the counts, by default it is log1p for RNA sequencing data, but for avoinding tranformation you can use identity
#'
#' @return A tibble filtered genes
#'
keep_variable_transcripts_SE = function(.data,
																		 top = 500,
																		 transform = log1p) {


	# Manage Inf
	top = min(top, .data %>% nrow)

	message(sprintf("Getting the %s most variable genes", top))

	x =
		.data %>%

		# Check if log transform is needed
		transform()


	s <- rowMeans((x - rowMeans(x, na.rm=TRUE)) ^ 2, na.rm=TRUE)
	o <- order(s, decreasing = TRUE)
	x <- x[o[1L:top], , drop = FALSE]
	variable_trancripts = rownames(x)

	.data[variable_trancripts,]

}


#' Drop redundant elements (e.g., samples) for which feature (e.g., genes) aboundances are correlated
#'
#' @keywords internal
#' @noRd
#'
#' 
#' 
#' @import tibble
#' @importFrom rlang :=
#'
#' @param .data A tibble
#' @param correlation_threshold A real number between 0 and 1
#' @param top An integer. How many top genes to select
#' @param of_samples A boolean
#'
#' @return A tibble with redundant elements removed
#'
#'
remove_redundancy_elements_through_correlation_SE <- function(.data,
																													 correlation_threshold = 0.9,
																													 of_samples = TRUE) {
	# Comply with CRAN NOTES
	. = NULL

	# Check if package is installed, otherwise install
	if (find.package("widyr", quiet = TRUE) %>% length %>% equals(0)) {
		message("Installing widyr needed for correlation analyses")
		install.packages("widyr", repos = "https://cloud.r-project.org")
	}

	# Get the redundant data frame
	.data %>%

		# check that there are non-NA genes for enough samples
		when(# First condition
			(.) %>% nrow == 0 	~ stop(
				"tidybulk says: In calculating correlation there is no gene that have non NA values is all samples"
			),

			# Second condition
			(.) %>% nrow < 100 ~ {
				message(
					"tidybulk says: In calculating correlation there is < 100 genes (that have non NA values) is all samples.
The correlation calculation might not be reliable"
				)
				.x
			},
			~ (.)) %>%

		as_tibble(rownames="transcript") %>%

		# Prepare the data frame
		gather(sample,abundance,-transcript) %>%

		when(
			of_samples ~ 	dplyr::rename(., rc = abundance,
																	element = sample,
																	feature = transcript) ,
			~ dplyr::rename(., rc = abundance,
											element = transcript,
											feature = sample)
		) %>%

		# Is this necessary?
		mutate_if(is.factor, as.character) %>%

		# Run pairwise correlation and return a tibble
		widyr::pairwise_cor(
			element,
			feature,
			rc,
			sort = TRUE,
			diag = FALSE,
			upper = FALSE
		) %>%
		filter(correlation > correlation_threshold) %>%
		distinct(item1) %>%
		pull(item1)

}

#' Identifies the closest pairs in a MDS context and return one of them
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom stats setNames
#' @importFrom stats dist
#'
#' @param .data A tibble
#' @param of_samples A boolean
#'
#' @return A tibble with pairs dropped
#'
#'
remove_redundancy_elements_though_reduced_dimensions_SE <-
	function(.data) {
		# This function identifies the closest pairs and return one of them


		# Calculate distances
		.data %>%
			dist() %>%

			# Prepare matrix
			as.matrix() %>%
			as_tibble(rownames = "sample a") %>%
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
			pull(1)

	}



#' Get differential transcription information to a tibble using edgeR.
#'
#' @keywords internal
#' @noRd
#'
#' 
#' 
#' @import tibble
#' @importFrom magrittr set_colnames
#' @importFrom stats model.matrix
#' @importFrom utils install.packages
#' @importFrom purrr when
#' @importFrom magrittr extract2
#'
#'
#' @param .data A tibble
#' @param .formula a formula with no response variable, referring only to numeric variables
#' @param .contrasts A character vector. See edgeR makeContrasts specification for the parameter `contrasts`. If contrasts are not present the first covariate is the one the model is tested against (e.g., ~ factor_of_interest)
#' @param method A string character. Either "edgeR_quasi_likelihood" (i.e., QLF), "edgeR_likelihood_ratio" (i.e., LRT)
#' @param test_above_log2_fold_change A positive real value. This works for edgeR and limma_voom methods. It uses the `treat` function, which tests that the difference in abundance is bigger than this threshold rather than zero \url{https://pubmed.ncbi.nlm.nih.gov/19176553}.
#' @param scaling_method A character string. The scaling method passed to the backend function (i.e., edgeR::calcNormFactors; "TMM","TMMwsp","RLE","upperquartile")
#' @param omit_contrast_in_colnames If just one contrast is specified you can choose to omit the contrast label in the colnames.
#'
#' @return A tibble with edgeR results
#'
get_differential_transcript_abundance_bulk_SE <- function(.data,
																											 .formula,
																											 .abundance = NULL,
																											 sample_annotation,
																											 .contrasts = NULL,
																											 method = "edgeR_quasi_likelihood",
																											 test_above_log2_fold_change = NULL,
																											 scaling_method = "TMM",
																											 omit_contrast_in_colnames = FALSE,
																											 prefix = "",
																											 ...) {

  .abundance = enquo(.abundance)
  
	# Check if omit_contrast_in_colnames is correctly setup
	if(omit_contrast_in_colnames & length(.contrasts) > 1){
		warning("tidybulk says: you can omit contrasts in column names only when maximum one contrast is present")
		omit_contrast_in_colnames = FALSE
	}

	# Create design matrix
	design =
		model.matrix(
			object = .formula,
			data = sample_annotation
		)

	# Print the design column names in case I want contrasts
	message(
		sprintf(
			"tidybulk says: The design column names are \"%s\"",
			design %>% colnames %>% paste(collapse = ", ")
		)
	)

	my_contrasts =
		.contrasts %>%
		ifelse_pipe(length(.) > 0,
								~ limma::makeContrasts(contrasts = .x, levels = design),
								~ NULL)

	# Check if package is installed, otherwise install
	if (find.package("edgeR", quiet = TRUE) %>% length %>% equals(0)) {
		message("Installing edgeR needed for differential transcript abundance analyses")
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
	
	edgeR_object =
		.data |> 
	  assay(my_assay) |> 
		edgeR::DGEList() %>%

		# Scale data if method is not "none"
		when(
			scaling_method != "none" ~ (.) %>% edgeR::calcNormFactors(method = scaling_method),
			~ (.)
		) %>%

		# select method
		when(
			tolower(method) ==  "edger_likelihood_ratio" ~ (.) %>% 	edgeR::estimateDisp(design) %>% edgeR::glmFit(design),
			tolower(method) ==  "edger_quasi_likelihood" ~ (.) %>% 	edgeR::estimateDisp(design) %>% edgeR::glmQLFit(design),
			tolower(method) == "edger_robust_likelihood_ratio" ~ (.) %>% edgeR::estimateGLMRobustDisp(design) %>% edgeR::glmFit(design)
		)

	# Return
	list(
		result_raw = edgeR_object,
		result =
			edgeR_object %>%

			# If I have multiple .contrasts merge the results
			ifelse_pipe(
				my_contrasts %>% is.null | omit_contrast_in_colnames,

				# Simple comparison
				~ .x %>%

					# select method
					when(
						!is.null(test_above_log2_fold_change) ~ (.) %>% edgeR::glmTreat(coef = 2, contrast = my_contrasts, lfc=test_above_log2_fold_change),
						tolower(method) %in%  c("edger_likelihood_ratio", "edger_robust_likelihood_ratio") ~ (.) %>% edgeR::glmLRT(coef = 2, contrast = my_contrasts) ,
						tolower(method) ==  "edger_quasi_likelihood" ~ (.) %>% edgeR::glmQLFTest(coef = 2, contrast = my_contrasts)
					)	%>%

					# Convert to tibble
					edgeR::topTags(n = Inf) %$%
					table %>%
					as_tibble(rownames = "transcript") %>%

					# # Mark DE genes
					# mutate(significant = FDR < significance_threshold) 	%>%

					# Arrange
					arrange(FDR),

				# Multiple comparisons
				~ {
					edgeR_obj = .x

					1:ncol(my_contrasts) %>%
						map_dfr(
							~ edgeR_obj %>%

								# select method
								when(
								    !is.null(test_above_log2_fold_change) ~ (.) %>% edgeR::glmTreat(coef = 2, contrast = my_contrasts[, .x], lfc=test_above_log2_fold_change),
									tolower(method) %in%  c("edger_likelihood_ratio", "edger_robust_likelihood_ratio") ~ (.) %>% edgeR::glmLRT(coef = 2, contrast = my_contrasts[, .x]) ,
									tolower(method) ==  "edger_quasi_likelihood" ~ (.) %>% edgeR::glmQLFTest(coef = 2, contrast = my_contrasts[, .x])
								)	%>%

								# Convert to tibble
								edgeR::topTags(n = Inf) %$%
								table %>%
								as_tibble(rownames = "transcript") %>%
								mutate(constrast = colnames(my_contrasts)[.x])
							# %>%
							#
							# # Mark DE genes
							# mutate(significant = FDR < significance_threshold)
						) %>%
						pivot_wider(values_from = -c(transcript, constrast),
												names_from = constrast, names_sep = "___")
				}
			)	 %>%

			# Attach prefix
			setNames(c(
				colnames(.)[1],
				sprintf("%s%s", prefix, colnames(.)[2:ncol(.)])
			))
	)




}


#' Get differential transcription information to a tibble using voom.
#'
#' @keywords internal
#' @noRd
#'
#' 
#' 
#' @import tibble
#' @importFrom magrittr set_colnames
#' @importFrom stats model.matrix
#' @importFrom utils install.packages
#' @importFrom purrr when
#'
#'
#' @param .data A tibble
#' @param .formula a formula with no response variable, referring only to numeric variables
#' @param .contrasts A character vector. See voom makeContrasts specification for the parameter `contrasts`. If contrasts are not present the first covariate is the one the model is tested against (e.g., ~ factor_of_interest)
#' @param method A string character. Either "limma_voom", "limma_voom_sample_weights"
#' @param scaling_method A character string. The scaling method passed to the backend function (i.e., edgeR::calcNormFactors; "TMM","TMMwsp","RLE","upperquartile")
#' @param omit_contrast_in_colnames If just one contrast is specified you can choose to omit the contrast label in the colnames.
#'
#' @return A tibble with voom results
#'
get_differential_transcript_abundance_bulk_voom_SE <- function(.data,
																														.formula,
																														.abundance = NULL,
																														sample_annotation,
																														.contrasts = NULL,
																														method = NULL,
																														test_above_log2_fold_change = NULL,
																														scaling_method = "TMM",
																														omit_contrast_in_colnames = FALSE,
																														prefix = "") {

  .abundance = enquo(.abundance)

	# Check if omit_contrast_in_colnames is correctly setup
	if(omit_contrast_in_colnames & length(.contrasts) > 1){
		warning("tidybulk says: you can omit contrasts in column names only when maximum one contrast is present")
		omit_contrast_in_colnames = FALSE
	}


	# Create design matrix
	design =
		model.matrix(
			object = .formula,
			data = sample_annotation
		)

	# Print the design column names in case I want contrasts
	message(
		sprintf(
			"tidybulk says: The design column names are \"%s\"",
			design %>% colnames %>% paste(collapse = ", ")
		)
	)

	my_contrasts =
		.contrasts %>%
		ifelse_pipe(length(.) > 0,
								~ limma::makeContrasts(contrasts = .x, levels = design),
								~ NULL)

	# Check if package is installed, otherwise install
	if (find.package("limma", quiet = TRUE) %>% length %>% equals(0)) {
		message("Installing limma needed for differential transcript abundance analyses")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("limma", ask = FALSE)
	}

	# If no assay is specified take first
	my_assay = ifelse(
	  quo_is_symbol(.abundance), 
	  quo_name(.abundance), 
	  .data |>
	    assayNames() |>
	    extract2(1)
	)
	
	voom_object =
		.data %>%

	  assay(my_assay) |> 
	  edgeR::DGEList() %>%

		# Scale data if method is not "none"
		when(
			scaling_method != "none" ~ (.) %>% edgeR::calcNormFactors(method = scaling_method),
			~ (.)
		) %>%

		# select method
		when(
			tolower(method) == "limma_voom" ~ (.) %>% limma::voom(design, plot=FALSE),
			tolower(method) == "limma_voom_sample_weights" ~ (.) %>% limma::voomWithQualityWeights(design, plot=FALSE)
		) %>%

		limma::lmFit(design)

	# Return
	list(
		result_raw = voom_object,
		result =
			voom_object %>%

			# If I have multiple .contrasts merge the results
			ifelse_pipe(
				my_contrasts %>% is.null | omit_contrast_in_colnames,

				# Simple comparison
				~ .x %>%

					# Contrasts
					limma::contrasts.fit(contrasts=my_contrasts, coefficients =  when(my_contrasts, is.null(.) ~ 2)) %>%
					limma::eBayes() %>%

			        when(

				    	!is.null(test_above_log2_fold_change) ~ (.) %>%
				    	    limma::treat(lfc=test_above_log2_fold_change) %>%
				    	    limma::topTreat(n = Inf),

				    	~ (.) %>% limma::topTable(n = Inf)

				    ) %>%

			        # Convert to tibble
					as_tibble(rownames = "transcript") %>%

					# # Mark DE genes
					# mutate(significant = adj.P.Val < significance_threshold) 	%>%

					# Arrange
					arrange(adj.P.Val),

				# Multiple comparisons
				~ {
					voom_obj = .x

					1:ncol(my_contrasts) %>%
						map_dfr(
							~ voom_obj %>%

								# Contrasts
								limma::contrasts.fit(contrasts=my_contrasts[, .x]) %>%
								limma::eBayes() %>%
						        when(

						            !is.null(test_above_log2_fold_change) ~ (.) %>%
						            limma::treat(lfc=test_above_log2_fold_change) %>%
						            limma::topTreat(n = Inf),

						            ~ (.) %>% limma::topTable(n = Inf)
						        ) %>%

							    # Convert to tibble
								as_tibble(rownames = "transcript") %>%
								mutate(constrast = colnames(my_contrasts)[.x])
							# %>%
							#
							# # Mark DE genes
							# mutate(significant = adj.P.Val < significance_threshold)
						) %>%
						pivot_wider(values_from = -c(transcript, constrast),
												names_from = constrast, names_sep = "___")
				}
			)	 %>%

			# Attach prefix
			setNames(c(
				colnames(.)[1],
				sprintf("%s%s", prefix, colnames(.)[2:ncol(.)])
			))
	)


}


#' Get differential transcription information to a tibble using glmmSeq
#'
#' @keywords internal
#' @noRd
#'
#' 
#' 
#' @import tibble
#' @importFrom magrittr set_colnames
#' @importFrom stats model.matrix
#' @importFrom utils install.packages
#' @importFrom purrr when
#'
#'
#' @param .data A tibble
#' @param .formula a formula with no response variable, referring only to numeric variables
#' @param .contrasts A character vector. See edgeR makeContrasts specification for the parameter `contrasts`. If contrasts are not present the first covariate is the one the model is tested against (e.g., ~ factor_of_interest)
#' @param method A string character. Either "edgeR_quasi_likelihood" (i.e., QLF), "edgeR_likelihood_ratio" (i.e., LRT)
#' @param scaling_method A character string. The scaling method passed to the backend function (i.e., edgeR::calcNormFactors; "TMM","TMMwsp","RLE","upperquartile")
#' @param omit_contrast_in_colnames If just one contrast is specified you can choose to omit the contrast label in the colnames.
#' @param ... Additional arguments for glmmSeq
#'
#' @return A tibble with glmmSeq results
#'
get_differential_transcript_abundance_glmmSeq_SE <- function(.data,
                                                            .formula,
                                                            .abundance = NULL,
                                                            .contrasts = NULL,
                                                            sample_annotation ,
                                                            method = "deseq2",
                                                            
                                                            test_above_log2_fold_change = NULL,
                                                            
                                                            scaling_method = "TMM",
                                                            omit_contrast_in_colnames = FALSE,
                                                            prefix = "",
                                                            ...) {
  
  .abundance = enquo(.abundance)
  
  # Check if contrasts are of the same form
  if(
    .contrasts %>% is.null %>% not() &
    .contrasts %>% class %>% equals("list") %>% not()
  )
    stop("tidybulk says: for DESeq2 the list of constrasts should be given in the form list(c(\"condition_column\",\"condition1\",\"condition2\")) i.e. list(c(\"genotype\",\"knockout\",\"wildtype\"))")
  
  # Check if omit_contrast_in_colnames is correctly setup
  if(omit_contrast_in_colnames & length(.contrasts) > 1){
    warning("tidybulk says: you can omit contrasts in column names only when maximum one contrast is present")
    omit_contrast_in_colnames = FALSE
  }
  
  # Check if package is installed, otherwise install
  if (find.package("edgeR", quiet = TRUE) %>% length %>% equals(0)) {
    message("tidybulk says: Installing edgeR needed for differential transcript abundance analyses")
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install("edgeR", ask = FALSE)
  }
  
  # Check if package is installed, otherwise install
  if (find.package("glmmSeq", quiet = TRUE) %>% length %>% equals(0)) {
    message("tidybulk says: Installing glmmSeq needed for differential transcript abundance analyses")
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install("glmmSeq", ask = FALSE)
  }
  
  # If no assay is specified take first
  my_assay = ifelse(
    quo_is_symbol(.abundance), 
    quo_name(.abundance), 
    .data |>
      assayNames() |>
      extract2(1)
  )
  
  metadata = 
    .data |> 
    colData() 
  
  counts = 
    .data %>%
    assay(my_assay)
  
  glmmSeq_object = 
    glmmSeq::glmmSeq( .formula,
                      countdata = counts ,
                      metadata =   metadata |> as.data.frame(),
                      dispersion = setNames(edgeR::estimateDisp(counts)$tagwise.dispersion, rownames(counts)),
                      progress = TRUE, 
                      method = method |> str_remove("(?i)^glmmSeq_" ),
                      ...
    ) 
  
  glmmSeq_object |> 
    summary() |> 
    as_tibble(rownames = "transcript") |>
    mutate(across(starts_with("P_"), list(adjusted = function(x) p.adjust(x, method="BH")), .names = "{.col}_{.fn}")) |> 
    
    # Attach attributes
    reattach_internals(.data) %>%
    
    # select method
    memorise_methods_used("glmmSeq") %>% 
    
    # # Add raw object
    # attach_to_internals(glmmSeq_object, "glmmSeq") %>%
    
    # Communicate the attribute added
    {
      rlang::inform("tidybulk says: to access the raw results (fitted GLM) do `attr(..., \"internals\")$glmmSeq`", .frequency_id = "Access DE results glmmSeq",  .frequency = "once")
      (.)
    }  %>%
    
    # Attach prefix
    setNames(c(
      colnames(.)[1],
      sprintf("%s%s", prefix, colnames(.)[2:ncol(.)])
    )) |> 
    
    list() |> 
    setNames("result") |> 
    c(list(result_raw = glmmSeq_object))
  
  
}


#' Get differential transcription information to a tibble using DESeq2
#'
#' @keywords internal
#' @noRd
#'
#' 
#' 
#' @import tibble
#' @importFrom magrittr set_colnames
#' @importFrom stats model.matrix
#' @importFrom utils install.packages
#' @importFrom purrr when
#'
#'
#' @param .data A tibble
#' @param .formula a formula with no response variable, referring only to numeric variables
#' @param .contrasts A character vector. See edgeR makeContrasts specification for the parameter `contrasts`. If contrasts are not present the first covariate is the one the model is tested against (e.g., ~ factor_of_interest)
#' @param method A string character. Either "edgeR_quasi_likelihood" (i.e., QLF), "edgeR_likelihood_ratio" (i.e., LRT)
#' @param scaling_method A character string. The scaling method passed to the backend function (i.e., edgeR::calcNormFactors; "TMM","TMMwsp","RLE","upperquartile")
#' @param omit_contrast_in_colnames If just one contrast is specified you can choose to omit the contrast label in the colnames.
#' @param ... Additional arguments for DESeq2
#'
#' @return A tibble with DESeq2 results
#'
get_differential_transcript_abundance_deseq2_SE <- function(.data,
                                                            .formula,
                                                            .abundance = NULL,
                                                            .contrasts = NULL,
                                                            method = "deseq2",

                                                            test_above_log2_fold_change = NULL,

                                                            scaling_method = "TMM",
                                                            omit_contrast_in_colnames = FALSE,
                                                            prefix = "",
                                                            ...) {

  .abundance = enquo(.abundance)
  
	# Check if contrasts are of the same form
	if(
		.contrasts %>% is.null %>% not() &
		.contrasts %>% class %>% equals("list") %>% not()
	)
		stop("tidybulk says: for DESeq2 the list of constrasts should be given in the form list(c(\"condition_column\",\"condition1\",\"condition2\")) i.e. list(c(\"genotype\",\"knockout\",\"wildtype\"))")

	# Check if omit_contrast_in_colnames is correctly setup
	if(omit_contrast_in_colnames & length(.contrasts) > 1){
		warning("tidybulk says: you can omit contrasts in column names only when maximum one contrast is present")
		omit_contrast_in_colnames = FALSE
	}

	# Check if package is installed, otherwise install
	if (find.package("DESeq2", quiet = TRUE) %>% length %>% equals(0)) {
		message("Installing DESeq2 needed for differential transcript abundance analyses")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("DESeq2", ask = FALSE)
	}

        if (is.null(test_above_log2_fold_change)) {
          test_above_log2_fold_change <- 0
        }
  
	my_contrasts = .contrasts

	# If no assay is specified take first
	my_assay = ifelse(
	  quo_is_symbol(.abundance), 
	  quo_name(.abundance), 
	  .data |>
	    assayNames() |>
	    extract2(1)
	)
	
	deseq2_object =
	  
		# DESeq2
		DESeq2::DESeqDataSetFromMatrix(
		  countData = .data |> assay(my_assay),
		  colData = colData(.data),
		  design = .formula
		) %>%
		DESeq2::DESeq(...)

	# Return
	list(
		result_raw = deseq2_object,
		result =
			# Read ft object
			deseq2_object %>%

			# If I have multiple .contrasts merge the results
			when(

				# Simple comparison continuous
				(my_contrasts %>% is.null ) &
					(deseq2_object@colData[,parse_formula(.formula)[1]] %>%
					 	class %in% c("numeric", "integer", "double")) 	~
					(.) %>%
					DESeq2::results(lfcThreshold=test_above_log2_fold_change) %>%
					as_tibble(rownames = "transcript"),

				# Simple comparison discrete
				my_contrasts %>% is.null 	~
					(.) %>%
					DESeq2::results(contrast = c(
						parse_formula(.formula)[1],
						deseq2_object@colData[,parse_formula(.formula)[1]] %>% as.factor() %>% levels %>% .[2],
						deseq2_object@colData[,parse_formula(.formula)[1]] %>% as.factor() %>% levels %>% .[1]
					), lfcThreshold=test_above_log2_fold_change) %>%
					as_tibble(rownames = "transcript"),

				# Simple comparison discrete
				my_contrasts %>% is.null %>% not() & omit_contrast_in_colnames	~
					(.) %>%
					DESeq2::results(contrast = my_contrasts[[1]], lfcThreshold=test_above_log2_fold_change)%>%
					as_tibble(rownames = "transcript"),

				# Multiple comparisons NOT USED AT THE MOMENT
				~ {
					deseq2_obj = (.)

					1:length(my_contrasts) %>%
						map_dfr(
							~ 	deseq2_obj %>%

								# select method
								DESeq2::results(contrast = my_contrasts[[.x]], lfcThreshold=test_above_log2_fold_change)	%>%

								# Convert to tibble
								as_tibble(rownames = "transcript") %>%
								mutate(constrast = sprintf("%s %s-%s", my_contrasts[[.x]][1], my_contrasts[[.x]][2], my_contrasts[[.x]][3]) )

						) %>%
						pivot_wider(values_from = -c(transcript, constrast),
												names_from = constrast, names_sep = "___")
				}
			)	 %>%

			# Attach prefix
			setNames(c(
				colnames(.)[1],
				sprintf("%s%s", prefix, colnames(.)[2:ncol(.)])
			))
	)


}

#'
#' @keywords internal
#' @noRd
#'
#' @importFrom stringr str_remove
#' @importFrom stringr str_replace_all
#'
multivariable_differential_tissue_composition_SE = function(
	deconvoluted,
	method,
	.my_formula,
	min_detected_proportion
){
	results_regression =
		deconvoluted %>%
		as_tibble(rownames = "sample") %>%

		# Replace 0s - before
		mutate(across(starts_with(method), function(.x) if_else(.x==0, min_detected_proportion, .x))) %>%
		mutate(across(starts_with(method), boot::logit)) %>%

		# Rename columns - after
		setNames(
			str_remove(colnames(.), sprintf("%s:", method)) %>%
				str_replace_all("[ \\(\\)]", "___")
		) %>%

		# Beta or Cox
		when(
			grepl("Surv", .my_formula) %>% any ~ {
				# Check if package is installed, otherwise install
				if (find.package("survival", quiet = TRUE) %>% length %>% equals(0)) {
					message("Installing betareg needed for analyses")
					install.packages("survival", repos = "https://cloud.r-project.org")
				}

				if (find.package("boot", quiet = TRUE) %>% length %>% equals(0)) {
					message("Installing boot needed for analyses")
					install.packages("boot", repos = "https://cloud.r-project.org")
				}

				(.) %>%
					survival::coxph(.my_formula, .)	%>%
					broom::tidy()
			} ,
			~ {
				(.) %>%
					lm(.my_formula, .) %>%
					broom::tidy() %>%
					filter(term != "(Intercept)")
			}
		)

	# Join results
	deconvoluted %>%
		as_tibble(rownames = "sample") %>%
		pivot_longer(
			names_prefix = sprintf("%s: ", method),
			cols = starts_with(method),
			names_to = ".cell_type",
			values_to = ".proportion"
		) %>%
		tidyr::nest(cell_type_proportions = -.cell_type) %>%
		bind_cols(
			results_regression %>%
				select(-term)
		)
}

univariable_differential_tissue_composition_SE = function(
	deconvoluted,
	method,
	.my_formula,
	min_detected_proportion
){
	deconvoluted %>%
		as_tibble(rownames = "sample") %>%

		# Test
		pivot_longer(
			names_prefix = sprintf("%s: ", method),
			cols = starts_with(method),
			names_to = ".cell_type",
			values_to = ".proportion"
		) %>%

		# Replace 0s
		mutate(.proportion_0_corrected = if_else(.proportion==0, min_detected_proportion, .proportion)) %>%

		# Test survival
		tidyr::nest(cell_type_proportions = -.cell_type) %>%
		mutate(surv_test = map(
			cell_type_proportions,
			~ {
				if(pull(., .proportion_0_corrected) %>% unique %>% length %>%  `<=` (3)) return(NULL)

				# See if regression if censored or not
				.x %>%
					when(
						grepl("Surv", .my_formula) %>% any ~ {
							# Check if package is installed, otherwise install
							if (find.package("survival", quiet = TRUE) %>% length %>% equals(0)) {
								message("Installing betareg needed for analyses")
								install.packages("survival", repos = "https://cloud.r-project.org")
							}

							if (find.package("boot", quiet = TRUE) %>% length %>% equals(0)) {
								message("Installing boot needed for analyses")
								install.packages("boot", repos = "https://cloud.r-project.org")
							}

							(.) %>%
								mutate(.proportion_0_corrected = .proportion_0_corrected  %>% boot::logit()) %>%
								survival::coxph(.my_formula, .)	%>%
								broom::tidy() %>%
								select(-term)
						} ,
						~ {
							# Check if package is installed, otherwise install
							if (find.package("betareg", quiet = TRUE) %>% length %>% equals(0)) {
								message("Installing betareg needed for analyses")
								install.packages("betareg", repos = "https://cloud.r-project.org")
							}
							(.) %>%
								betareg::betareg(.my_formula, .) %>%
								broom::tidy() %>%
								filter(component != "precision") %>%
								pivot_wider(names_from = term, values_from = c(estimate, std.error, statistic,   p.value)) %>%
								select(-c(`std.error_(Intercept)`, `statistic_(Intercept)`, `p.value_(Intercept)`)) %>%
								select(-component)
						}
					)
			}
		)) %>%

		unnest(surv_test, keep_empty = TRUE)
}
