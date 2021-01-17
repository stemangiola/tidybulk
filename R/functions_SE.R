#' Get K-mean clusters to a tibble
#'
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom stats kmeans
#' @importFrom rlang :=
#'
#' @param .data A tibble
#' @param .abundance A column symbol with the value the clustering is based on (e.g., `count`)
#' @param .feature A column symbol. The column that is represents entities to cluster (i.e., normally samples)
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally genes)
#' @param of_samples A boolean
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#' @param ... Further parameters passed to the function kmeans
#'
#' @return A tibble with additional columns
#'
#'
get_clusters_kmeans_bulk_SE <-
	function(.data,
					 of_samples = TRUE,
					 log_transform = TRUE,
					 ...) {
		
		# Check if centers is in dots
		dots_args = rlang::dots_list(...)
		if ("centers" %in% names(dots_args) %>% not())
			stop("tidybulk says: for kmeans you need to provide the \"centers\" integer argument")
		
		.data %>%
			
			# Check if log transform is needed
			when(log_transform ~ log1p(.), ~ (.) ) %>%
			
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
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang :=
#' @importFrom utils install.packages
#'
#' @param .data A tibble
#' @param .abundance A column symbol with the value the clustering is based on (e.g., `count`)
#' @param .feature A column symbol. The column that is represents entities to cluster (i.e., normally samples)
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally genes)
#' @param of_samples A boolean
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#' @param ... Further parameters passed to the function kmeans
#'
#' @return A tibble with additional columns
#'
get_clusters_SNN_bulk_SE <-
	function(.data,
					 of_samples = TRUE,
					 log_transform = TRUE,
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
#'
#' @import dplyr
#' @import tidyr
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
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#'
#' @return A tibble with additional columns
#'
#'
get_reduced_dimensions_MDS_bulk_SE <-
	function(.data,
					 .dims = 2,
					 top = 500,
					 of_samples = TRUE,
					 log_transform = TRUE) {
		# Comply with CRAN NOTES
		. = NULL
		
		# Get components from dims
		components = 1:.dims
		
		mds_object = limma::plotMDS(.data, ndim = .dims, plot = FALSE, top = top)
		
		# Return
		list(
			raw_result = mds_object,
			result = 
				# Parse results
				mds_object %$%	cmdscale.out %>%
				as.data.frame %>%
				setNames(c(sprintf("Dim%s", 1:.dims))) 
		)
		
		
	}

#' Get principal component information to a tibble using PCA
#'
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
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
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
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
					 log_transform = TRUE,
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
				select(sample, sprintf("PC%s", components)) 
		)
		
		
	}

#' Get principal component information to a tibble using tSNE
#'
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
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
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#' @param ... Further parameters passed to the function Rtsne
#'
#' @return A tibble with additional columns
#'
get_reduced_dimensions_TSNE_bulk_SE <-
	function(.data,
					 .dims = 2,
					 top = 500,
					 of_samples = TRUE,
					 log_transform = TRUE,
					 ...) {
		# Comply with CRAN NOTES
		. = NULL
		
		
		# Evaluate ...
		arguments <- list(...)
		if (!"check_duplicates" %in% names(arguments))
			arguments = arguments %>% c(check_duplicates = TRUE)
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
				.data %>% distinct(!!.element) %>% nrow %>% sum(-1)
			) / 3 / 2) %>% floor() %>% min(30))
		
		# If not enough samples stop
		if (arguments$perplexity <= 2)
			stop("tidybulk says: You don't have enough samples to run tSNE")
		
		# Calculate the most variable genes, from plotMDS Limma
		do.call(Rtsne::Rtsne, c(list(.data), arguments)) %$%
			Y %>%
			as_tibble(.name_repair = "minimal") %>%
			setNames(c("tSNE1", "tSNE2")) %>%
			
			# add element name
			dplyr::mutate(!!.element := .data %>% rownames) %>%
			select(!!.element, everything()) 
		
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
#'
#' @param .data A tibble
#' @param .sample A character name of the sample column
#' @param .transcript A character name of the transcript/gene column
#' @param .abundance A character name of the read count column
#' @param top An integer. How many top genes to select
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#'
#' @return A tibble filtered genes
#'
keep_variable_transcripts_SE = function(.data,
																		 top = 500,
																		 log_transform = TRUE) {

	
	# Manage Inf
	top = min(top, .data %>% nrow)
	
	message(sprintf("Getting the %s most variable genes", top))
	
	x =
		.data %>%

		# Check if log transform is needed
		when(log_transform ~ log1p(.), ~ (.) ) 
	
	s <- rowMeans((x - rowMeans(x)) ^ 2)
	o <- order(s, decreasing = TRUE)
	x <- x[o[1L:top], , drop = FALSE]
	variable_trancripts = rownames(x)
	
	.data[variable_trancripts,]
		
}


#' Drop redundant elements (e.g., samples) for which feature (e.g., genes) aboundances are correlated
#'
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
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

#' Get adjusted count for some batch effect
#'
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr set_colnames
#' @importFrom stats model.matrix
#' @importFrom stats as.formula
#' @importFrom utils install.packages
#' @importFrom stats rnorm
#'
#' @param .data A tibble
#' @param .formula a formula with no response variable, of the kind ~ factor_of_interest + batch
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#' @param ... Further parameters passed to the function sva::ComBat
#'
#' @return A tibble with adjusted counts
#'
#'
get_adjusted_counts_for_unwanted_variation_bulk <- function(.data,
																														.formula,
																														sample_annotation,
																														...) {

	


	mat_for_combat <-
		.data 
	
	

	mat_for_combat 
	
}