


#' Create tt object from tibble
#'
#' @importFrom rlang enquo
#' @importFrom magrittr %>%
#'
#' @param .data A tibble
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param .abundance_scaled The name of the transcript/gene scaled abundance column
#'
#' @return A tibble with an additional column
#'
#'
create_tt_from_tibble_bulk = function(.data,
																			.sample,
																			.transcript,
																			.abundance,
																			.abundance_scaled = NULL) {
	# Make col names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	.abundance_scaled = enquo(.abundance_scaled)

	.data %>%

		# Add tt_columns attribute
		add_tt_columns(!!.sample,!!.transcript,!!.abundance,!!.abundance_scaled) %>%

		# Add class
		add_class("tt") %>%
		add_class("tidyBulk")
}



#' Convert bam/sam files to a tidy gene transcript counts data frame
#'
#' @importFrom purrr reduce
#'
#' @param file_names A character vector
#' @param genome A character string
#' @param ... Further parameters passed to the function Rsubread::featureCounts
#'
#' @return A tibble of gene counts
#'
create_tt_from_bam_sam_bulk <-
	function(file_names, genome = "hg38", ...) {
		# This function uses Subread to count the gene features,
		# annotate gene features with symbols, and
		# convert the data frame to tibble format

		n_cores <- system("nproc", intern = TRUE) %>%
			as.integer() %>%
			`-`(2)

		file_names %>%

			# Run subread
			Rsubread::featureCounts(annot.inbuilt = genome,
															nthreads = n_cores,
															...) %>%

			# Anonymous function
			# input: Subread results
			# output edgeR::DGEList object
			{
				edgeR::DGEList(
					counts = (.)$counts,
					genes = (.)$annotation[, c("GeneID", "Length")],
					samples = (.)$stat %>% as_tibble() %>% gather(sample, temp,-Status) %>% spread(Status, temp)
				)
			} %>%

			# Anonymous function
			# input: edgeR::DGEList object
			# output: edgeR::DGEList object with added transcript symbol
			{
				dge <- (.)
				dge$genes$transcript <-
					AnnotationDbi::mapIds(
						org.Hs.eg.db::org.Hs.eg.db,
						keys = as.character(dge$genes$GeneID),
						column = "SYMBOL",
						keytype = "ENTREZID",
						multiVals = "first"
					)

				dge
			} %>%

			# Anonymous function
			# input: annotated edgeR::DGEList object
			# output: tibble
			{
				reduce(
					list(
						(.) %$% counts %>% as_tibble(rownames = "GeneID") %>% mutate(GeneID = GeneID %>% as.integer()) %>% gather(sample, `count`,-GeneID),
						(.) %$% genes %>% select(GeneID, transcript) %>% as_tibble(),
						(.) %$% samples %>% as_tibble()
					),
					dplyr::left_join
				) %>%
					rename(entrez = GeneID) %>%
					mutate(entrez = entrez %>% as.character())
			}
	}


#' Get count per million for TMM scaling.
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang :=
#' @importFrom stats median
#'
#' @param .data A tibble
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param cpm_threshold A real positive number. Minimum counts per million required for a selected proportion of samples
#'
#' @return A tibble with an additional column
add_scaled_counts_bulk.get_cpm <- function(.data,
																					 .sample = `sample`,
																					 .transcript = `transcript`,
																					 .abundance = `count`,
																					 cpm_threshold = 0.5) {
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	if (cpm_threshold < 0)
		stop("The parameter cpm_threshold must be > 0")

	# Add cmp and cmp threshold to the data set, and return
	.data %>%
		dplyr::left_join(
			(.) %>%

				######################################################################
			# I don't know what this was for, but is dangerous. Delete after check
			# select(-contains("ct")) %>%
			######################################################################

			select(!!.transcript,!!.sample,!!.abundance) %>%
				spread(!!.sample,!!.abundance) %>%
				drop_na() %>%
				do(
					dplyr::bind_cols(
						!!.transcript := (.) %>% pull(!!.transcript),
						tibble::as_tibble((.) %>% select(-!!.transcript) %>% edgeR::cpm())
					)
				) %>%
				gather(!!.sample, cpm,-!!.transcript) %>%
				mutate(cpm_threshold = cpm_threshold),
			by = c(quo_name(.transcript), quo_name(.sample))
		) %>%

		# Attach attributes
		reattach_internals(.data)


}

filter_transcript_high_prop_cpm = function(.data,
																					 .sample,
																					 .transcript,
																					 cpm_threshold,
																					 prop_threshold) {
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)

	.data %>%
		select(!!.sample,!!.transcript, cpm) %>%
		spread(!!.sample, cpm) %>%
		as_matrix(rownames = !!.transcript) %>%
		`>` (cpm_threshold) %>%
		rowSums() %>%
		`<` ((max(.) * prop_threshold) %>% floor)  %>%
		which %>%
		names

}


#' Drop lowly tanscribed genes for TMM normalization
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang :=
#' @importFrom stats median
#'
#' @param .data A tibble
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param cpm_threshold A real positive number. Minimum counts per million required for a selected proportion of samples
#' @param prop_threshold A real positive number between 0 and 1. It is the threshold of proportion of samples for each transcripts/genes that have to be characterised by a cmp bigger than the threshold to be included for scaling procedure.
#'
#' @return A tibble filtered
add_scaled_counts_bulk.get_low_expressed <- function(.data,
																										 .sample = `sample`,
																										 .transcript = `transcript`,
																										 .abundance = `count`,
																										 cpm_threshold = 0.5,
																										 prop_threshold = 3 / 4) {
	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
	.sample = col_names$.sample
	.transcript = col_names$.transcript
	.abundance = col_names$.abundance

	if (cpm_threshold < 0)
		stop("The parameter cpm_threshold must be > 0")
	if (prop_threshold < 0 |
			prop_threshold > 1)
		stop("The parameter prop_threshold must be between 0 and 1")

	.data %>%

		# Prepare the data frame
		select(!!.transcript,!!.sample,!!.abundance) %>%

		# Calculate cpm
		add_scaled_counts_bulk.get_cpm(!!.sample,!!.transcript,!!.abundance) %>%
		# Filter based on how many samples have a gene below the threshold
		filter_transcript_high_prop_cpm(!!.sample,
																		!!.transcript,
																		cpm_threshold = cpm_threshold,
																		prop_threshold = prop_threshold) %>%
		# Attach attributes
		reattach_internals(.data)
}


#' Calculate the norm factor with calcNormFactor from limma
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang :=
#' @importFrom stats setNames
#'
#' @param .data A tibble
#' @param reference A reference matrix, not sure if used anymore
#' @param cpm_threshold A real positive number. Minimum counts per million required for a selected proportion of samples
#' @param prop_threshold A real positive number between 0 and 1. It is the threshold of proportion of samples for each transcripts/genes that have to be characterised by a cmp bigger than the threshold to be included for scaling procedure.
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param method A string character
#'
#'
#' @return A list including the filtered data frame and the normalization factors
add_scaled_counts_bulk.calcNormFactor <- function(.data,
																									reference = NULL,
																									cpm_threshold = 0.5,
																									prop_threshold = 3 / 4,
																									.sample = `sample`,
																									.transcript = `transcript`,
																									.abundance = `count`,
																									method) {
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	error_if_log_transformed(.data,!!.abundance)

	# Get list of low transcribed genes
	gene_to_exclude <-
		add_scaled_counts_bulk.get_low_expressed(
			.data %>%
				filter(!!.sample != "reference"),!!.sample,!!.transcript,!!.abundance,
			cpm_threshold = cpm_threshold,
			prop_threshold = prop_threshold
		)

	# Check if transcript after filtering is 0
	if (length(gene_to_exclude) == .data %>%
			dplyr::distinct(!!.transcript) %>%
			nrow()) {
		stop("The gene expression matrix has been filtered completely for lowly expressed genes")
	}

	# Get data frame for the higly transcribed transcripts
	df.filt <-
		.data %>%
		dplyr::filter(!(!!.transcript %in% gene_to_exclude)) %>%
		droplevels()



	# List of low abundant transcripts
	gene_to_exclude = gene_to_exclude

	# scaled data set
	nf =
		tibble::tibble(
			# Sample factor
			sample = factor(levels(df.filt %>% pull(!!.sample))),

			# scaled data frame
			nf = edgeR::calcNormFactors(
				df.filt %>%
					tidyr::spread(!!.sample,!!.abundance) %>%
					tidyr::drop_na() %>%
					dplyr::select(-!!.transcript),
				refColumn = which(reference == factor(levels(
					df.filt %>% pull(!!.sample)
				))),
				method = method
			)
		) %>%

		setNames(c(quo_name(.sample), "nf")) %>%

		# Add the statistics about the number of genes filtered
		dplyr::left_join(
			df.filt %>%
				dplyr::group_by(!!.sample) %>%
				dplyr::summarise(tot_filt = sum(!!.abundance, na.rm = TRUE)) %>%
				dplyr::mutate(!!.sample := as.factor(as.character(!!.sample))),
			by = quo_name(.sample)
		)

	# Return
	list(gene_to_exclude = gene_to_exclude, nf = nf) %>%

		# Attach attributes
		reattach_internals(.data)
}

#' Get a tibble with scaled counts using TMM
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr equals
#' @importFrom rlang :=
#' @importFrom stats median
#' @importFrom utils installed.packages
#' @importFrom utils install.packages
#'
#' @param .data A tibble
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param cpm_threshold A real positive number. Minimum counts per million required for a selected proportion of samples
#' @param prop_threshold A real positive number between 0 and 1. It is the threshold of proportion of samples for each transcripts/genes that have to be characterised by a cmp bigger than the threshold to be included for scaling procedure.
#' @param method A character string
#' @param reference_selection_function A function between median, mean and max
#'
#' @return A tibble including additional columns
#'
#'
get_scaled_counts_bulk <- function(.data,
																	 .sample = NULL,
																	 .transcript = NULL,
																	 .abundance = NULL,
																	 cpm_threshold = 0.5,
																	 prop_threshold = 3 / 4,
																	 method = "TMM",
																	 reference_selection_function = median) {
	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
	.sample = col_names$.sample
	.transcript = col_names$.transcript
	.abundance = col_names$.abundance

	# Check if package is installed, otherwise install
	if ("edgeR" %in% rownames(installed.packages()) == FALSE) {
		writeLines("Installing edgeR needed for differential transcript abundance analyses")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("edgeR")
	}

	# Set column name for value scaled
	value_scaled = as.symbol(sprintf("%s scaled",  quo_name(.abundance)))

	# Reformat input data set
	df <-
		.data %>%

		# # Check input types
		# error_if_wrong_input(
		#   as.list(environment())[-1],
		#   c("spec_tbl_df",  "quosure",  "quosure",  "quosure")
		# ) %>%

		# Stop if any counts is NA
		error_if_counts_is_na(!!.abundance) %>%

		# Stop if there are duplicated transcripts
		error_if_duplicated_genes(!!.sample,!!.transcript,!!.abundance) %>%

		# Rename
		dplyr::select(!!.sample,!!.transcript,!!.abundance) %>%
		#setNames(c("!!.sample", "gene", "count")) %>%

		# Set samples and genes as factors
		dplyr::mutate(!!.sample := factor(!!.sample),!!.transcript := factor(!!.transcript))


	# Get norm factor object
	reference <-
		df %>%
		group_by(!!.sample) %>%
		summarise(sum = sum(!!.abundance)) %>%
		mutate(med = reference_selection_function(sum)) %>%
		mutate(diff = abs(sum - med)) %>%
		arrange(diff) %>%
		head(n = 1) %>%
		pull(!!.sample) %>%
		as.character()

	nf_obj <-
		add_scaled_counts_bulk.calcNormFactor(
			df,
			reference,
			cpm_threshold,
			prop_threshold,
			.sample = !!.sample,
			.transcript = !!.transcript,
			.abundance = !!.abundance,
			method
		)

	# Calculate normalization factors

	nf <- nf_obj$nf %>%
		dplyr::left_join(
			df %>%
				group_by(!!.sample) %>%
				summarise(tot = sum(!!.abundance, na.rm = TRUE)) %>%
				ungroup() %>%
				dplyr::mutate(!!.sample := as.factor(as.character(!!.sample))),
			by = quo_name(.sample)
		) %>%
		mutate(multiplier =
					 	1 /
					 	(tot_filt * nf) *
					 	((.) %>% filter(!!.sample == reference) %>% pull(tot))) %>%

		# I have correct the strange behaviour of edgeR of reference
		# sample not being 1
		ifelse_pipe(
			"reference" %in% ((.) %>% pull(!!.sample)),
			~ .x %>%
				mutate(
					multiplier =
						multiplier /
						(.) %>%
						filter(!!.sample == "reference") %>%
						pull(multiplier)
				)
		)

	# Return
	df_norm =
		df %>%
		dplyr::mutate(!!.sample := as.factor(as.character(!!.sample))) %>%
		dplyr::left_join(nf, by = quo_name(.sample)) %>%

		# Calculate scaled values
		dplyr::mutate(!!value_scaled := !!.abundance * multiplier) %>%

		# Format df for join
		dplyr::select(!!.sample, !!.transcript, !!value_scaled,
									everything()) %>%
		dplyr::mutate(`filter out low counts` = !!.transcript %in% nf_obj$gene_to_exclude) %>%
		dplyr::select(-!!.abundance,-tot,-tot_filt) %>%
		dplyr::rename(TMM = nf) %>%
		arrange(!!.sample,!!.transcript)
	#dplyr::select(-!!.sample,-!!.transcript)

	# Attach attributes
	df_norm %>%
		add_tt_columns(!!.sample,!!.transcript,!!.abundance,!!(function(x, v)
			enquo(v))(x,!!value_scaled))

}

#' Add a tibble with scaled counts using TMM
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom stats median
#'
#' @param .data A tibble
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param cpm_threshold A real positive number. Minimum counts per million required for a selected proportion of samples
#' @param prop_threshold A real positive number between 0 and 1. It is the threshold of proportion of samples for each transcripts/genes that have to be characterised by a cmp bigger than the threshold to be included for scaling procedure.
#' @param method A character string
#' @param reference_selection_function A function between median, mean and max
#'
#' @return A tibble including additional columns
#'
#'
add_scaled_counts_bulk <- function(.data,
																	 .sample = NULL,
																	 .transcript = NULL,
																	 .abundance = NULL,
																	 cpm_threshold = 0.5,
																	 prop_threshold = 3 / 4,
																	 method = "TMM",
																	 reference_selection_function = median) {
	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
	.sample = col_names$.sample
	.transcript = col_names$.transcript
	.abundance = col_names$.abundance


	.data_norm =
		.data %>%
		get_scaled_counts_bulk(
			.sample = !!.sample,
			.transcript = !!.transcript,
			.abundance = !!.abundance,
			cpm_threshold = cpm_threshold,
			prop_threshold = prop_threshold,
			method = method,
			reference_selection_function = reference_selection_function
		) %>%
		arrange(!!.sample,!!.transcript)

	.data %>%
		arrange(!!.sample,!!.transcript) %>%

		# Add scaled data set
		bind_cols(.data_norm %>%
								select(-one_of(quo_name(.sample)), -one_of(quo_name(.transcript))))		%>%

		# Attach attributes
		reattach_internals(.data_norm)


}


#' Get differential transcription information to a tibble using edgeR.
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr set_colnames
#' @importFrom stats model.matrix
#' @importFrom utils installed.packages
#' @importFrom utils install.packages
#'
#'
#' @param .data A tibble
#' @param .formula a formula with no response variable, referring only to numeric variables
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param .coef An integer. See edgeR specifications
#' @param .contrasts A character vector. See edgeR makeContrasts specification for the parameter `contrasts`
#' @param significance_threshold A real between 0 and 1
#' @param cpm_threshold A real positive number. Minimum counts per million required for a selected proportion of samples
#' @param prop_threshold A real positive number between 0 and 1. It is the threshold of proportion of samples for each transcripts/genes that have to be characterised by a cmp bigger than the threshold to be included for scaling procedure.
#' @param fill_missing_values A boolean. Whether to fill missing sample/transcript values with the median of the transcript. This is rarely needed.
#' @param scaling_method A character string. The scaling method passed to the backend function (i.e., edgeR::calcNormFactors)
#'
#' @return A tibble with edgeR results
#'
get_differential_transcript_abundance_bulk <- function(.data,
																											 .formula,
																											 .sample = NULL,
																											 .transcript = NULL,
																											 .abundance = NULL,
																											 .coef = 2,
																											 .contrasts = NULL,
																											 significance_threshold = 0.05,
																											 cpm_threshold = 0.5,
																											 prop_threshold = 3 / 4,
																											 fill_missing_values = FALSE,
																											 scaling_method = "TMM") {
	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
	.sample = col_names$.sample
	.transcript = col_names$.transcript
	.abundance = col_names$.abundance

	# distinct_at is not released yet for dplyr, thus we have to use this trick
	df_for_edgeR <- .data %>%

		# Stop if any counts is NA
		error_if_counts_is_na(!!.abundance) %>%

		# Stop if there are duplicated transcripts
		error_if_duplicated_genes(!!.sample,!!.transcript,!!.abundance) %>%

		# Prepare the data frame
		select(!!.transcript,
					 !!.sample,
					 !!.abundance,
					 one_of(parse_formula(.formula))) %>%
		distinct() %>%

		# Check if data rectangular
		ifelse_pipe(
			(.) %>% check_if_data_rectangular(!!.sample,!!.transcript,!!.abundance, type = "soft") %>% `!` &
				!fill_missing_values,
			~ .x %>% eliminate_sparse_transcripts(!!.transcript)
		)

	# Check if at least two samples for each group
	if (# If I have some discrete covariates
		df_for_edgeR %>%
		select(one_of(parse_formula(.formula))) %>%
		select_if(function(col)
			is.character(col) | is.factor(col) | is.logical(col)) %>%
		ncol %>% `>` (0) &

		# If I have at least 2 samples per group
		df_for_edgeR %>%
		select(!!.sample, one_of(parse_formula(.formula))) %>%
		distinct %>%
		count(!!as.symbol(parse_formula(.formula))) %>%
		distinct(n) %>%
		pull(1) %>%
		min %>%
		`<` (2))
	stop("You need at least two replicated for each condition for edgeR to work")

	# Create design matrix
	design =
		model.matrix(
			object = .formula,
			data = df_for_edgeR %>% select(!!.sample, one_of(parse_formula(.formula))) %>% distinct %>% arrange(!!.sample)
		)

	# Print the design column names in case I want constrasts
	message(
		sprintf(
			"tidyBulk says: The design column names are \"%s\" in case you are interested in contrasts",
			design %>% colnames %>% paste(collapse = ", ")
		)
	)

	my_contrasts =
		.contrasts %>%
		ifelse_pipe(length(.) > 0,
								~ limma::makeContrasts(contrasts = .x, levels = design),
								~ NULL)

	# Check if package is installed, otherwise install
	if ("edgeR" %in% rownames(installed.packages()) == FALSE) {
		writeLines("Installing edgeR needed for differential transcript abundance analyses")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("edgeR")
	}

	df_for_edgeR.filt <-
		df_for_edgeR %>%
		select(!!.transcript,!!.sample,!!.abundance) %>%
		mutate(
			`filter out low counts` = !!.transcript %in% add_scaled_counts_bulk.get_low_expressed(
				.,
				!!.sample,
				!!.transcript,
				!!.abundance,
				cpm_threshold = cpm_threshold,
				prop_threshold = prop_threshold
			)
		)

	edgeR_object =
		df_for_edgeR.filt %>%
		filter(!`filter out low counts`) %>%
		select(!!.transcript,!!.sample,!!.abundance) %>%
		spread(!!.sample,!!.abundance) %>%
		as_matrix(rownames = !!.transcript) %>%

		# If fill missing values
		ifelse_pipe(fill_missing_values,
								~ .x %>% fill_NA_with_row_median) %>%

		edgeR::DGEList(counts = .) %>%
		edgeR::calcNormFactors(method = scaling_method) %>%
		edgeR::estimateGLMCommonDisp(design) %>%
		edgeR::estimateGLMTagwiseDisp(design) %>%
		edgeR::glmFit(design)

	edgeR_object %>%

		# If I have multiple .contrasts merge the results
		ifelse_pipe(
			my_contrasts %>% is.null | ncol(my_contrasts) < 2,

			# Simple comparison
			~ .x %>%
				edgeR::glmLRT(coef = .coef, contrast = my_contrasts) %>%
				edgeR::topTags(n = 999999) %$%
				table %>%
				as_tibble(rownames = quo_name(.transcript)) %>%

				# Mark DE genes
				mutate(is_de = FDR < significance_threshold) 	%>%

				# Arrange
				arrange(FDR),

			# Multiple comparisons
			~ {
				edgeR_obj = .x

				1:ncol(my_contrasts) %>%
					map_dfr(
						~ edgeR_obj %>%
							edgeR::glmLRT(coef = .coef, contrast = my_contrasts[, .x]) %>%
							edgeR::topTags(n = 999999) %$%
							table %>%
							as_tibble(rownames = quo_name(.transcript)) %>%
							mutate(constrast = colnames(my_contrasts)[.x]) %>%

							# Mark DE genes
							mutate(is_de = FDR < significance_threshold)
					) %>%
					pivot_wider(values_from = -c(!!.transcript, constrast),
											names_from = constrast)
			}
		)	 %>%

		# Add filtering info
		full_join(df_for_edgeR.filt %>%
								select(!!.transcript, `filter out low counts`) %>%
								distinct()) %>%


		# Attach attributes
		reattach_internals(.data) %>%

		# Add raw object
		attach_to_internals(edgeR_object, "edgeR") %>%
		# Communicate the attribute added
		{
			message(
				"tidyBulk says: to access the raw results (glmFit) do `attr(..., \"tt_internals\")$edgeR`"
			)
			(.)
		}
}

#' Add differential transcription information to a tibble using edgeR.
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr set_colnames
#'
#'
#' @param .data A tibble
#' @param .formula a formula with no response variable, referring only to numeric variables
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param .coef An integer. See edgeR specifications
#' @param .contrasts A character vector. See edgeR makeContrasts specification for the parameter `contrasts`
#' @param significance_threshold A real between 0 and 1
#' @param cpm_threshold A real positive number. Minimum counts per million required for a selected proportion of samples
#' @param prop_threshold A real positive number between 0 and 1. It is the threshold of proportion of samples for each transcripts/genes that have to be characterised by a cmp bigger than the threshold to be included for scaling procedure.
#' @param fill_missing_values A boolean. Whether to fill missing sample/transcript values with the median of the transcript. This is rarely needed.
#' @param scaling_method A character string. The scaling method passed to the backend function (i.e., edgeR::calcNormFactors)
#'
#' @return A tibble with differential_transcript_abundance results
#'
#'
add_differential_transcript_abundance_bulk <- function(.data,
																											 .formula,
																											 .sample = NULL,
																											 .transcript = NULL,
																											 .abundance = NULL,
																											 .coef = 2,
																											 .contrasts = NULL,
																											 significance_threshold = 0.05,
																											 cpm_threshold = 0.5,
																											 prop_threshold = 3 / 4,
																											 fill_missing_values = FALSE,
																											 scaling_method = "TMM") {
	# Comply with CRAN NOTES
	. = NULL
	FDR = NULL

	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
	.sample = col_names$.sample
	.transcript = col_names$.transcript
	.abundance = col_names$.abundance

	.data_processed =
		.data %>%
		get_differential_transcript_abundance_bulk(
			.formula,
			.sample = !!.sample,
			.transcript = !!.transcript,
			.abundance = !!.abundance,
			.coef = .coef,
			.contrasts = .contrasts,
			significance_threshold = significance_threshold,
			cpm_threshold = cpm_threshold,
			prop_threshold = prop_threshold,
			fill_missing_values = fill_missing_values,
			scaling_method = scaling_method
		)

	.data %>%
		dplyr::left_join(.data_processed) %>%

		# Arrange
		ifelse_pipe(.contrasts %>% is.null,
								~ .x %>% arrange(FDR))	%>%

		# Attach attributes
		reattach_internals(.data_processed)
}


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
#' symbol_to_entrez(tidyBulk::counts_mini, .transcript = transcript, .sample = sample)
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
	if ("org.Hs.eg.db" %in% rownames(installed.packages()) == FALSE) {
		writeLines("Installing org.Hs.eg.db needed for differential transcript abundance analyses")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("org.Hs.eg.db")
	}

	.data %>%
		dplyr::left_join(
			# Get entrez mapping 1:1
			AnnotationDbi::mapIds(
				org.Hs.eg.db::org.Hs.eg.db,
				.data %>% distinct(!!.transcript) %>% pull(1),
				'ENTREZID',
				'SYMBOL'
			) %>%
				enframe(name = quo_name(.transcript), value = "entrez") %>%
				filter(entrez %>% is.na %>% `!`) %>%
				group_by(!!.transcript) %>%
				slice(1) %>%
				ungroup()
		)

}

#' Get gene enrichment analyses using EGSEA
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr set_colnames
#' @importFrom purrr map2_dfr
#' @importFrom stats model.matrix
#' @importFrom utils installed.packages
#' @importFrom utils install.packages
#'
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
#' @return A tibble with edgeR results
#'
analyse_gene_enrichment_bulk_EGSEA <- function(.data,
																							 .formula,
																							 .sample = NULL,
																							 .entrez,
																							 .abundance = NULL,
																							 .contrasts = NULL,
																							 species,
																							 cores = 10) {
	# Comply with CRAN NOTES
	. = NULL

	# Get column names
	.sample = enquo(.sample)
	.abundance = enquo(.abundance)
	col_names = get_sample_counts(.data, .sample, .abundance)
	.sample = col_names$.sample
	.abundance = col_names$.abundance

	.entrez = enquo(.entrez)

	# distinct_at is not released yet for dplyr, thus we have to use this trick
	df_for_edgeR <- .data %>%

		# Stop if any counts is NA
		error_if_counts_is_na(!!.abundance) %>%

		# Stop if there are duplicated transcripts
		error_if_duplicated_genes(!!.sample,!!.entrez,!!.abundance) %>%

		# Prepare the data frame
		select(!!.entrez, !!.sample, !!.abundance,
					 one_of(parse_formula(.formula))) %>%
		distinct() %>%

		# Add entrez from symbol
		filter(!!.entrez %>% is.na %>% `!`)

	# Check if at least two samples for each group
	if (df_for_edgeR %>%
			select(!!.sample, one_of(parse_formula(.formula))) %>%
			distinct %>%
			count(!!as.symbol(parse_formula(.formula))) %>%
			distinct(n) %>%
			pull(1) %>%
			min %>%
			`<` (2))
		stop("You need at least two replicated for each condition for EGSEA to work")

	# Create design matrix
	design =
		model.matrix(
			object = .formula,
			data = df_for_edgeR %>% select(!!.sample, one_of(parse_formula(.formula))) %>% distinct %>% arrange(!!.sample)
		) %>%
		magrittr::set_colnames(c("(Intercept)",
														 (.) %>% colnames %>% `[` (-1)))

	# # Check if package is installed, otherwise install
	# if ("EGSEA" %in% rownames(installed.packages()) == FALSE) {
	# 	writeLines("Installing EGSEA needed for differential transcript abundance analyses")
	# 	if (!requireNamespace("BiocManager", quietly = TRUE))
	# 		install.packages("BiocManager", repos = "https://cloud.r-project.org")
	# 	BiocManager::install("EGSEA")
	# }

	# Check if package is installed, otherwise install
	if ("EGSEA" %in% rownames(installed.packages()) == FALSE) {
		writeLines("EGSEA not installed. Please install it with.")
		writeLines("BiocManager::install(\"EGSEA\")")
	}
	if (!"EGSEA" %in% (.packages())) {
		writeLines("EGSEA package not loaded. Please run library(\"EGSEA\")")
	}

	df_for_edgeR.filt <-
		df_for_edgeR %>%
		select(!!.entrez,!!.sample,!!.abundance) %>%
		mutate(
			`filter out low counts` = !!.entrez %in% add_scaled_counts_bulk.get_low_expressed(.,!!.sample,!!.entrez,!!.abundance)
		) %>%
		filter(!`filter out low counts`) %>%

		# Make sure transcrpt names are adjacent
		arrange(!!.entrez)

	dge =
		df_for_edgeR.filt %>%
		select(!!.entrez,!!.sample,!!.abundance) %>%
		spread(!!.sample,!!.abundance) %>%
		as_matrix(rownames = !!.entrez) %>%
		edgeR::DGEList(counts = .)

	idx =  buildIdx(entrezIDs = rownames(dge), species = species)

	res =
		dge %>%

		# Calculate weights
		limma::voom(design, plot = FALSE) %>%

		# Execute EGSEA
		egsea(
			contrasts = .contrasts,
			gs.annots = idx,
			# symbolsMap=
			# 	df_for_edgeR.filt %>%
			# 	select(entrez, !!.transcript) %>%
			# 	distinct() %>%
			# 	arrange(match(entrez, rownames(dge))) %>%
			# 	setNames(c("FeatureID", "Symbols")),
			baseGSEAs = egsea.base()[-c(6, 7, 8, 9, 12)],
			sort.by = "med.rank",
			num.threads = cores
		)

	res@results %>%
		map2_dfr(
			res@results %>% names,
			~ .x[[1]][[1]] %>%
				as_tibble(rownames = "pathway") %>%
				mutate(data_base = .y)
		) %>%
		arrange(med.rank) %>%
		select(data_base, pathway, everything())


}

#' Get K-mean clusters to a tibble
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
get_clusters_kmeans_bulk <-
	function(.data,
					 .element = NULL,
					 .feature = NULL,
					 .abundance = NULL,
					 of_samples = TRUE,
					 log_transform = TRUE,
					 ...) {
		# Check if centers is in dots
		dots_args = rlang::dots_list(...)
		if ("centers" %in% names(dots_args) %>% `!`)
			stop("tidyBulk says: for kmeans you need to provide the \"centers\" integer argument")

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

		.data %>%

			# Through error if some counts are NA
			error_if_counts_is_na(!!.abundance) %>%

			# Prepare data frame
			distinct(!!.feature,!!.element,!!.abundance) %>%

			# Check if log tranfrom is needed
			ifelse_pipe(log_transform,
									~ .x %>% dplyr::mutate(!!.abundance := !!.abundance %>%  `+`(1) %>%  log())) %>%

			# Prepare data frame for return
			spread(!!.feature,!!.abundance) %>%
			as_matrix(rownames = !!.element) %>%

			# Wrap the do.call because of the centers check
			{
				do.call(kmeans, list(x = (.), iter.max = 1000) %>% c(dots_args))
			}	 %$%
			cluster %>%
			as.list() %>%
			as_tibble() %>%
			gather(!!.element, `cluster kmeans`) %>%
			mutate(`cluster kmeans` = `cluster kmeans` %>% as.factor()) %>%

			# Attach attributes
			reattach_internals(.data)
	}

#' Add K-mean clusters to a tibble
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom stats kmeans
#'
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
add_clusters_kmeans_bulk <-
	function(.data,
					 .element = NULL,
					 .feature = NULL,
					 .abundance,
					 of_samples = TRUE,
					 log_transform = TRUE,
					 ...) {
		# Comply with CRAN NOTES
		. = NULL

		# Get column names
		.element = enquo(.element)
		.feature = enquo(.feature)
		col_names = get_elements_features(.data, .element, .feature, of_samples)
		.element = col_names$.element
		.feature = col_names$.feature

		.abundance = enquo(.abundance)

		.data %>%
			dplyr::left_join(
				(.) %>%
					get_clusters_kmeans_bulk(
						.abundance = !!.abundance,
						.element = !!.element,
						.feature = !!.feature,
						log_transform = log_transform,
						...
					)
			) %>%

			# Attach attributes
			reattach_internals(.data)
	}

#' Get SNN shared nearest neighbour clusters to a tibble
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang :=
#' @importFrom utils installed.packages
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
get_clusters_SNN_bulk <-
	function(.data,
					 .element = NULL,
					 .feature = NULL,
					 .abundance,
					 of_samples = TRUE,
					 log_transform = TRUE,
					 ...) {
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

		# Check if package is installed, otherwise install
		if ("Seurat" %in% rownames(installed.packages()) == FALSE) {
			writeLines("Installing Seurat")
			install.packages("Seurat", repos = "https://cloud.r-project.org")
		}
		if ("KernSmooth" %in% rownames(installed.packages()) == FALSE) {
			writeLines("Installing KernSmooth")
			install.packages("KernSmooth", repos = "https://cloud.r-project.org")
		}

		my_df =
			.data %>%

			# Through error if some counts are NA
			error_if_counts_is_na(!!.abundance) %>%

			# Prepare data frame
			distinct(!!.element,!!.feature,!!.abundance) %>%

			# Check if log tranfrom is needed
			#ifelse_pipe(log_transform, ~ .x %>% dplyr::mutate(!!.abundance := !!.abundance %>%  `+`(1) %>%  log())) %>%

			# Prepare data frame for return
			spread(!!.element,!!.abundance)

		my_df %>%
			data.frame(row.names = quo_name(.feature)) %>%
			Seurat::CreateSeuratObject() %>%
			Seurat::ScaleData(display.progress = TRUE,
												num.cores = 4,
												do.par = TRUE) %>%
			Seurat::FindVariableFeatures(selection.method = "vst") %>%
			Seurat::RunPCA(npcs = 30) %>%
			Seurat::FindNeighbors() %>%
			Seurat::FindClusters(method = "igraph", ...) %>%
			`[[` ("seurat_clusters") %>%
			as_tibble(rownames = quo_name(.element)) %>%
			rename(`cluster SNN` = seurat_clusters) %>%
			dplyr::mutate(!!.element := gsub("\\.", "-",!!.element)) %>%

			# Attach attributes
			reattach_internals(.data)
	}

#' Add SNN shared nearest neighbour clusters to a tibble
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom stats kmeans
#'
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
add_clusters_SNN_bulk <-
	function(.data,
					 .element = NULL,
					 .feature = NULL,
					 .abundance,
					 of_samples = TRUE,
					 log_transform = TRUE,
					 ...) {
		# Get column names
		.element = enquo(.element)
		.feature = enquo(.feature)
		col_names = get_elements_features(.data, .element, .feature, of_samples)
		.element = col_names$.element
		.feature = col_names$.feature

		.abundance = enquo(.abundance)

		.data %>%
			dplyr::left_join(
				(.) %>%
					get_clusters_SNN_bulk(
						.abundance = !!.abundance,
						.element = !!.element,
						.feature = !!.feature,
						log_transform = log_transform,
						...
					)
			) %>%

			# Attach attributes
			reattach_internals(.data)
	}

#' Get dimensionality information to a tibble using MDS
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
get_reduced_dimensions_MDS_bulk <-
	function(.data,
					 .element = NULL,
					 .feature = NULL,
					 .abundance = NULL,
					 .dims = 2,
					 top = 500,
					 of_samples = TRUE,
					 log_transform = TRUE) {
		# Comply with CRAN NOTES
		. = NULL

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

		# Get components from dims
		components = 1:.dims

		mds_object =
			.data %>%

			# Through error if some counts are NA
			error_if_counts_is_na(!!.abundance) %>%

			# Filter lowly transcribed (I have to avoid the use of scaling function)
			filter_abundant(!!.element, !!.feature,!!.abundance) %>%
			distinct(!!.feature,!!.element,!!.abundance) %>%

			# Check if logtansform is needed
			ifelse_pipe(log_transform,
									~ .x %>% dplyr::mutate(!!.abundance := !!.abundance %>% `+`(1) %>%  log())) %>%

			# Stop any column is not if not numeric or integer
			ifelse_pipe(
				(.) %>% select(!!.abundance) %>% summarise_all(class) %>% `%in%`(c("numeric", "integer")) %>% `!`() %>% any(),
				~ stop(".abundance must be numerical or integer")
			) %>%
			spread(!!.element,!!.abundance) %>%
			as_matrix(rownames = !!.feature, do_check = FALSE) %>%
			limma::plotMDS(ndim = .dims, plot = FALSE, top = top)

		# Pase results
		mds_object %$%	cmdscale.out %>%
			as.data.frame %>%
			as_tibble(rownames = quo_name(.element)) %>%
			setNames(c(quo_name(.element), sprintf("Dim%s", 1:.dims))) %>%


			# Attach attributes
			reattach_internals(.data) %>%

			# Add raw object
			attach_to_internals(mds_object, "MDS") %>%
			# Communicate the attribute added
			{
				message("tidyBulk says: to access the raw results do `attr(..., \"tt_internals\")$MDS`")
				(.)
			}

	}

#' Add dimensionality information to a tibble using MDS
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
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
add_reduced_dimensions_MDS_bulk <-
	function(.data,
					 .element = NULL,
					 .feature = NULL,
					 .abundance = NULL,
					 .dims = 2,
					 top = 500,
					 of_samples = TRUE,
					 log_transform = TRUE) {
		# Comply with CRAN NOTES
		. = NULL

		# Get column names
		.element = enquo(.element)
		.feature = enquo(.feature)
		.abundance = enquo(.abundance)
		col_names = get_elements_features_abundance(.data, .element, .feature, .abundance, of_samples)
		.element = col_names$.element
		.feature = col_names$.feature
		.abundance = col_names$.abundance

		.data_processed =
			.data %>%
			get_reduced_dimensions_MDS_bulk(
				.abundance = !!.abundance,
				.dims = .dims,
				.element = !!.element,
				.feature = !!.feature,
				top = top,
				of_samples = of_samples,
				log_transform = log_transform
			)

		.data %>%	dplyr::left_join(.data_processed,	by = quo_name(.element)) %>%

			# Attach attributes
			reattach_internals(.data_processed)
	}

#' Get principal component information to a tibble using PCA
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang :=
#' @importFrom stats prcomp
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
get_reduced_dimensions_PCA_bulk <-
	function(.data,
					 .element = NULL,
					 .feature = NULL,

					 .abundance  = NULL,
					 .dims = 2,
					 top = 500,
					 of_samples = TRUE,
					 log_transform = TRUE,
					 scale = TRUE,
					 ...) {
		# Comply with CRAN NOTES
		. = NULL

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

		# Get components from dims
		components = 1:.dims

		prcomp_obj =
			.data %>%

			# Through error if some counts are NA
			error_if_counts_is_na(!!.abundance) %>%

			# Prepare data frame
			distinct(!!.feature,!!.element,!!.abundance) %>%

			# Check if logtansform is needed
			ifelse_pipe(log_transform,
									~ .x %>% dplyr::mutate(!!.abundance := !!.abundance %>% `+`(1) %>%  log())) %>%

			# Stop any column is not if not numeric or integer
			ifelse_pipe(
				(.) %>% select(!!.abundance) %>% summarise_all(class) %>% `%in%`(c("numeric", "integer")) %>% `!`() %>% any(),
				~ stop(".abundance must be numerical or integer")
			) %>%

			# Filter most variable genes
			filter_variable_transcripts(!!.element,!!.feature,!!.abundance, top) %>%

			spread(!!.element,!!.abundance) %>%

			drop_na %>% # Is this necessary?

			# check that there are non-NA genes for enough samples
			ifelse2_pipe(# First condition
				(.) %>% nrow == 0,

				# Second condition
				(.) %>% nrow < 100,

				# First function
				~ stop(
					"In calculating PCA there is no gene that have non NA values is all samples"
				),

				# Second function
				~ {
					warning(
						"
						In PCA correlation there is < 100 genes that have non NA values is all samples.
						The correlation calculation would not be reliable,
						we suggest to partition the dataset for sample clusters.
						"
					)
					.x
				}) %>%

			# Transform to matrix
			as_matrix(rownames = !!.feature, do_check = FALSE) %>%

			# Calculate principal components
			prcomp(scale = scale, ...)

		prcomp_obj %>%

			# Anonymous function - Prints fraction of variance
			# input: PCA object
			# output: PCA object
			{
				writeLines("Fraction of variance explained by the selected principal components")

				(.) %$% sdev %>% `^` (2) %>% # Eigen value
					`/` (sum(.)) %>%
					`[` (components) %>%
					enframe() %>%
					select(-name) %>%
					rename(`Fraction of variance` = value) %>%
					mutate(PC = components) %>%
					print(n = 9999999)

				(.)

			} %$%

			# Parse the PCA results to a tibble
			rotation %>%
			as_tibble(rownames = quo_name(.element)) %>%
			select(!!.element, sprintf("PC%s", components)) %>%

			# Attach attributes
			reattach_internals(.data) %>%

			# Add raw object
			attach_to_internals(prcomp_obj, "PCA") %>%
			# Communicate the attribute added
			{
				message("tidyBulk says: to access the raw results do `attr(..., \"tt_internals\")$PCA`")
				(.)
			}

	}

#' Add principal component information to a tibble using PCA
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
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
add_reduced_dimensions_PCA_bulk <-
	function(.data,
					 .element = NULL,
					 .feature = NULL,

					 .abundance = NULL,
					 .dims = 2,
					 top = 500,
					 of_samples = TRUE,
					 log_transform = TRUE,
					 scale = TRUE,
					 ...) {
		# Comply with CRAN NOTES
		. = NULL

		# Get column names
		.element = enquo(.element)
		.feature = enquo(.feature)
		.abundance = enquo(.abundance)
		col_names = get_elements_features_abundance(.data, .element, .feature, .abundance, of_samples)
		.element = col_names$.element
		.feature = col_names$.feature
		.abundance = col_names$.abundance

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

		.data %>%
			dplyr::left_join(.data_processed,	by = quo_name(.element)) %>%

			# Attach attributes
			reattach_internals(.data_processed)
	}

#' Get principal component information to a tibble using tSNE
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang :=
#' @importFrom stats setNames
#' @importFrom utils installed.packages
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
#' @param scale A boolean
#' @param ... Further parameters passed to the function Rtsne
#'
#' @return A tibble with additional columns
#'
get_reduced_dimensions_TSNE_bulk <-
	function(.data,
					 .element = NULL,
					 .feature = NULL,

					 .abundance = NULL,
					 .dims = 2,
					 top = 500,
					 of_samples = TRUE,
					 log_transform = TRUE,
					 scale = TRUE,
					 ...) {
		# Comply with CRAN NOTES
		. = NULL

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

		# Evaluate ...
		arguments <- list(...)
		if (!"check_duplicates" %in% names(arguments))
			arguments = arguments %>% c(check_duplicates = TRUE)
		if (!"verbose" %in% names(arguments))
			arguments = arguments %>% c(verbose = TRUE)
		if (!"dims" %in% names(arguments))
			arguments = arguments %>% c(dims = .dims)


		# Check if package is installed, otherwise install
		if ("Rtsne" %in% rownames(installed.packages()) == FALSE) {
			writeLines("Installing Rtsne")
			install.packages("Rtsne", repos = "https://cloud.r-project.org")
		}

		# Set perprexity to not be too high
		if (!"perplexity" %in% names(arguments))
			arguments = arguments %>% c(perplexity = ((
				.data %>% distinct(!!.element) %>% nrow %>% sum(-1)
			) / 3 / 2) %>% floor() %>% min(30))

		# If not enough samples stop
		if (arguments$perplexity <= 2)
			stop("You don't have enough samples to run tSNE")

		# Calculate the most variable genes, from plotMDS Limma


		df_tsne =
			.data %>%

			# Check if duplicates
			error_if_duplicated_genes(!!.element,!!.feature,!!.abundance)  %>%

			# Filter NA symbol
			filter(!!.feature %>% is.na %>% `!`) %>%

			# Prepare data frame
			distinct(!!.feature,!!.element,!!.abundance) %>%

			# Check if data rectangular
			ifelse_pipe(
				(.) %>% check_if_data_rectangular(!!.element,!!.feature,!!.abundance, type = "soft"),
				~ .x %>% eliminate_sparse_transcripts(!!.feature)
			) %>%

			# Check if logtansform is needed
			ifelse_pipe(log_transform,
									~ .x %>% dplyr::mutate(!!.abundance := !!.abundance %>% `+`(1) %>%  log())) %>%

			# Filter most variable genes
			filter_variable_transcripts(!!.element,!!.feature,!!.abundance, top) %>%

			spread(!!.feature,!!.abundance) %>%
			# select(-sample) %>%
			# distinct %>%
			as_matrix(rownames = quo_name(.element))

		do.call(Rtsne::Rtsne, c(list(df_tsne), arguments)) %$%
			Y %>%
			as_tibble(.name_repair = "minimal") %>%
			setNames(c("tSNE1", "tSNE2")) %>%

			# add element name
			dplyr::mutate(!!.element := df_tsne %>% rownames) %>%
			select(!!.element, everything()) %>%

			# Attach attributes
			reattach_internals(.data)

	}

#' Add principal component information to a tibble using tSNE
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
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
#' @param ... Further parameters passed to the function tSNE
#'
#' @return A tibble with additional columns
#'
#'
add_reduced_dimensions_TSNE_bulk <-
	function(.data,
					 .element = NULL,
					 .feature = NULL,

					 .abundance = NULL,
					 .dims = 2,
					 top = 500,
					 of_samples = TRUE,
					 log_transform = TRUE,
					 scale = TRUE,
					 ...) {
		# Comply with CRAN NOTES
		. = NULL

		# Get column names
		.element = enquo(.element)
		.feature = enquo(.feature)
		.abundance = enquo(.abundance)
		col_names = get_elements_features_abundance(.data, .element, .feature, .abundance, of_samples)
		.element = col_names$.element
		.feature = col_names$.feature
		.abundance = col_names$.abundance

		.data %>%
			dplyr::left_join(
				(.) %>%
					get_reduced_dimensions_TSNE_bulk(
						.abundance = !!.abundance,
						.dims = .dims,
						.element = !!.element,
						.feature = !!.feature,
						top = top,
						log_transform = log_transform,
						scale = TRUE,
						...
					),
				by = quo_name(.element)
			) %>%

			# Attach attributes
			reattach_internals(.data)
	}


#' Get rotated dimensions of two principal components or MDS dimension of choice, of an angle
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang quo_is_null
#'
#'
#' @param .data A tibble
#' @param dimension_1_column A column symbol. The column of the dimension 1
#' @param dimension_2_column   A column symbol. The column of the dimension 2
#' @param rotation_degrees A real number between 0 and 360
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally samples)
#' @param of_samples A boolean
#' @param dimension_1_column_rotated A column symbol. The column of the dimension 1 rotated
#' @param dimension_2_column_rotated   A column symbol. The column of the dimension 2 rotated
#'
#' @return A tibble with additional rotated columns
#'
#'
get_rotated_dimensions =
	function(.data,
					 dimension_1_column,
					 dimension_2_column,
					 rotation_degrees,
					 .element = NULL,
					 of_samples = TRUE,
					 dimension_1_column_rotated = NULL,
					 dimension_2_column_rotated = NULL) {
		# Get column names
		.element = enquo(.element)
		col_names = get_elements(.data, .element)
		.element = col_names$.element

		# Parse other colnames
		dimension_1_column = enquo(dimension_1_column)
		dimension_2_column = enquo(dimension_2_column)
		dimension_1_column_rotated = enquo(dimension_1_column_rotated)
		dimension_2_column_rotated = enquo(dimension_2_column_rotated)

		if (.data %>%
				distinct(!!.element, !!dimension_1_column, !!dimension_2_column) %>%
				count(!!.element, !!dimension_1_column, !!dimension_2_column) %>%
				pull(n) %>%
				max %>%
				`>` (1))
			stop(sprintf(
				"%s must be unique for each row for the calculation of rotation",
				quo_name(.element)
			))

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

		# Function that rotates a 2D space of a arbitrary angle
		rotation = function(m, d) {
			r = d * pi / 180
			((dplyr::bind_rows(
				c(`1` = cos(r), `2` = -sin(r)),
				c(`1` = sin(r), `2` = cos(r))
			) %>% as_matrix) %*% m)
		}

		# Sanity check of the angle selected
		if (rotation_degrees %>% between(-360, 360) %>% `!`)
			stop("rotation_degrees must be between -360 and 360")

		# Return
		.data %>%
			distinct(!!.element, !!dimension_1_column, !!dimension_2_column) %>%
			as_matrix(rownames = !!.element) %>% t %>%
			rotation(rotation_degrees) %>%
			as_tibble() %>%
			mutate(`rotated dimensions` =
						 	c(
						 		quo_name(dimension_1_column_rotated),
						 		quo_name(dimension_2_column_rotated)
						 	)) %>%
			gather(!!.element, value,-`rotated dimensions`) %>%
			spread(`rotated dimensions`, value) %>%

			# Attach attributes
			reattach_internals(.data)

	}

#' Add Rotated dimensions of two principal components or MDS dimensions, of an angle
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang quo_is_null
#'
#'
#' @param .data A tibble
#' @param dimension_1_column A column symbol. The column of the dimension 1
#' @param dimension_2_column   A column symbol. The column of the dimension 2
#' @param rotation_degrees A real number between 0 and 360
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally samples)
#' @param of_samples A boolean
#' @param dimension_1_column_rotated A column symbol. The column of the dimension 1 rotated
#' @param dimension_2_column_rotated   A column symbol. The column of the dimension 2 rotated
#'
#' @return A tibble with additional rotated columns
#'
#'
add_rotated_dimensions =
	function(.data,
					 dimension_1_column,
					 dimension_2_column,
					 rotation_degrees,
					 .element = NULL,
					 of_samples = TRUE,
					 dimension_1_column_rotated = NULL,
					 dimension_2_column_rotated = NULL) {
		# Comply with CRAN NOTES
		. = NULL

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

		.data %>%
			dplyr::left_join(
				(.) %>%
					get_rotated_dimensions(
						dimension_1_column = !!dimension_1_column,
						dimension_2_column = !!dimension_2_column,
						rotation_degrees = rotation_degrees,
						.element = !!.element,
						of_samples = of_samples,
						dimension_1_column_rotated = !!dimension_1_column_rotated,
						dimension_2_column_rotated = !!dimension_2_column_rotated
					),
				by = quo_name(.element)
			) %>%

			# Attach attributes
			reattach_internals(.data)
	}

#' Aggregates multiple counts from the same samples (e.g., from isoforms)
#' This function aggregates counts over samples, concatenates other character columns, and averages other numeric columns
#'
#' @importFrom dplyr summarise_all
#' @importFrom dplyr bind_rows
#' @importFrom magrittr %$%
#' @importFrom rlang :=
#'
#' @param .data A tibble
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param aggregation_function A function for counts aggregation (e.g., sum)
#' @param keep_integer A boolean
#'
#' @return A tibble with aggregated genes and annotation
#'
#'
aggregate_duplicated_transcripts_bulk =
	function(.data,
					 .sample = NULL,
					 .transcript = NULL,
					 .abundance = NULL,
					 aggregation_function = sum,

					 keep_integer = TRUE) {
		# Comply with CRAN NOTES
		. = NULL

		# Get column names
		.sample = enquo(.sample)
		.transcript = enquo(.transcript)
		.abundance = enquo(.abundance)
		col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
		.sample = col_names$.sample
		.transcript = col_names$.transcript
		.abundance = col_names$.abundance

		# Robust paste function that preserves NAs
		paste3 <- function(..., sep = ", ") {
			L <- list(...)
			L <- lapply(L, function(x) {
				x[is.na(x)] <- ""
				x
			})
			ret <- gsub(paste0("(^", sep, "|", sep, "$)"),
									"",
									gsub(paste0(sep, sep), sep,
											 do.call(paste, c(
											 	L, list(sep = sep)
											 ))))
			is.na(ret) <- ret == ""
			ret
		}

		# Through warning if there are logicals of factor in the data frame
		# because they cannot be merged if they are not unique
		if ((lapply(.data, class) %>% unlist %in% c("logical", "factor")) %>% any) {
			warning("for aggregation fctors and logical columns were converted to character")
			writeLines("Converted to characters")
			lapply(.data, class) %>% unlist %>% `[` (. %in% c("logical", "factor") %>% which) %>% print
		}

		# Select which are the numerical columns
		numerical_columns =
			.data %>%
			ungroup() %>%
			select_if(is.numeric) %>%
			select(-!!.abundance) %>%

			# If scaled add the column to the exclusion
			ifelse_pipe((
				".abundance_scaled" %in% (.data %>% get_tt_columns() %>% names) &&
					# .data %>% get_tt_columns() %$% .abundance_scaled %>% is.null %>% `!` &&
					quo_name(.data %>% get_tt_columns() %$% .abundance_scaled) %in% (.data %>% colnames)
			),
			~ .x %>% select(-!!(
				.data %>% get_tt_columns() %$% .abundance_scaled
			)))	%>%
			colnames() %>%
			c("n_aggr")

		# ggregates read .data over samples, concatenates other character columns, and averages other numeric columns
		.data %>%

			# Through error if some counts are NA
			error_if_counts_is_na(!!.abundance) %>%

			# transform logials and factors
			mutate_if(is.factor, as.character) %>%
			mutate_if(is.logical, as.character) %>%

			# Add the nuber of duplicates for each gene
			dplyr::left_join((.) %>% count(!!.sample,!!.transcript, name = "n_aggr"),
											 by = c(quo_name(.sample), quo_name(.transcript))) %>%

			# Anonymous function - binds the unique and the reduced genes,
			# in the way we have to reduce redundancy just for the duplicated genes
			# input: tibble
			# output tibble distinct
			{
				dplyr::bind_rows(
					# Unique symbols
					(.) %>%
						filter(n_aggr == 1),

					# Duplicated symbols
					(.) %>%
						filter(n_aggr > 1) %>%
						group_by(!!.sample,!!.transcript) %>%
						dplyr::mutate(!!.abundance := !!.abundance %>% aggregation_function()) %>%

						# If scaled abundance exists aggragate that as well
						ifelse_pipe((
							".abundance_scaled" %in% (.data %>% get_tt_columns() %>% names) &&
								# .data %>% get_tt_columns() %$% .abundance_scaled %>% is.null %>% `!` &&
								quo_name(.data %>% get_tt_columns() %$% .abundance_scaled) %in% (.data %>% colnames)
						),
						~ {
							.abundance_scaled = .data %>% get_tt_columns() %$% .abundance_scaled
							.x %>% dplyr::mutate(!!.abundance_scaled := !!.abundance_scaled %>% aggregation_function())
						}) %>%

						mutate_at(vars(numerical_columns), mean) %>%
						mutate_at(
							vars(
								-group_cols(),
								-contains(quo_name(.abundance)),
								-!!numerical_columns
							),
							list( ~ paste3(unique(.), collapse = ", "))
						) %>%
						distinct()
				)
			} %>%

			# Rename column of number of duplicates for each gene
			rename(`merged transcripts` = n_aggr) %>%

			# Attach attributes
			reattach_internals(.data)

	}

#' Drop redundant elements (e.g., samples) for which feature (e.g., genes) aboundances are correlated
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang :=
#'
#' @param .data A tibble
#' @param .abundance A column symbol with the value the clustering is based on (e.g., `count`)
#' @param correlation_threshold A real number between 0 and 1
#' @param top An integer. How many top genes to select
#' @param .feature A column symbol. The column that is represents entities to cluster (i.e., normally genes)
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally samples)
#' @param of_samples A boolean
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#'
#' @return A tibble with redundant elemens removed
#'
#'
remove_redundancy_elements_through_correlation <- function(.data,
																													 .element = NULL,
																													 .feature = NULL,
																													 .abundance = NULL,
																													 correlation_threshold = 0.9,
																													 top = Inf,
																													 of_samples = TRUE,
																													 log_transform = FALSE) {
	# Comply with CRAN NOTES
	. = NULL

	# Get column names
	.element = enquo(.element)
	.feature = enquo(.feature)
	.abundance = enquo(.abundance)
	col_names = get_elements_features_abundance(.data, .element, .feature, .abundance, of_samples)
	.element = col_names$.element
	.feature = col_names$.feature
	.abundance = col_names$.abundance

	# Check if package is installed, otherwise install
	if ("widyr" %in% rownames(installed.packages()) == FALSE) {
		writeLines("Installing widyr needed for correlation analyses")
		install.packages("widyr")
	}

	# Get the redundant data frame
	.data.correlated =
		.data %>%

		# Stop if any counts is NA
		error_if_counts_is_na(!!.abundance) %>%

		# Stop if there are duplicated transcripts
		error_if_duplicated_genes(!!.element,!!.feature,!!.abundance) %>%

		# Prepare the data frame
		select(!!.feature,!!.element,!!.abundance) %>%

		# Filter variable genes
		filter_variable_transcripts(!!.element,!!.feature,!!.abundance, top = top) %>%

		# Check if logtansform is needed
		ifelse_pipe(log_transform,
								~ .x %>% dplyr::mutate(!!.abundance := !!.abundance %>% `+`(1) %>%  log())) %>%
		distinct() %>%
		spread(!!.element,!!.abundance) %>%
		drop_na() %>%

		# check that there are non-NA genes for enough samples
		ifelse2_pipe(# First condition
			(.) %>% nrow == 0,

			# Second condition
			(.) %>% nrow < 100,

			# First function
			~ stop(
				"In calculating correlation there is no gene that have non NA values is all samples"
			),

			# Second function
			~ {
				warning(
					"
					In calculating correlation there is < 100 genes that have non NA values is all samples.
					The correlation calculation would not be reliable,
					we suggest to partition the dataset for sample clusters.
					"
				)
				.x
			}) %>%

		# Prepare the data frame
		gather(!!.element,!!.abundance,-!!.feature) %>%
		dplyr::rename(rc := !!.abundance,
									sample := !!.element,
									transcript := !!.feature) %>% # Is rename necessary?
		mutate_if(is.factor, as.character) %>%

		# Run pairwise correlation and return a tibble
		widyr::pairwise_cor(
			sample,
			transcript,
			rc,
			sort = TRUE,
			diag = FALSE,
			upper = FALSE
		) %>%
		filter(correlation > correlation_threshold) %>%
		distinct(item1) %>%
		dplyr::rename(!!.element := item1)

	# Return non redudant data frame
	.data %>% anti_join(.data.correlated) %>%

		# Attach attributes
		reattach_internals(.data)
}

#' Identifies the closest pairs in a MDS contaxt and return one of them
#'
#' @importFrom stats setNames
#' @importFrom stats dist
#'
#' @param .data A tibble
#' @param Dim_a_column A column symbol. The column of one principal component
#' @param Dim_b_column A column symbol. The column of another principal component
#' @param .element A column symbol. The column that is represents entities to cluster (i.e., normally samples)
#' @param of_samples A boolean
#'
#' @return A tibble with pairs dropped
#'
#'
remove_redundancy_elements_though_reduced_dimensions <-
	function(.data,
					 Dim_a_column,
					 Dim_b_column,
					 .element = NULL,
					 of_samples = TRUE) {
		# This function identifies the closest pairs and return one of them

		# Get column names
		.element = enquo(.element)
		col_names = get_elements(.data, .element)
		.element = col_names$.element

		Dim_a_column = enquo(Dim_a_column)
		Dim_b_column = enquo(Dim_b_column)

		# Find redundant samples
		.data.redundant =

			# Calculate distances
			.data %>%
			select(!!.element,!!Dim_a_column,!!Dim_b_column) %>%
			distinct() %>%
			as_matrix(rownames = !!.element) %>%
			dist() %>%

			# Prepare matrix
			as.matrix() %>% as_tibble(rownames = "sample a") %>%
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
			select(1) %>%

			# Set the column names
			setNames(quo_name(.element))

		# Drop samples that are correlated with others and return
		.data %>% anti_join(.data.redundant) %>%

			# Attach attributes
			reattach_internals(.data)
	}

# #' after wget, this function merges hg37 and hg38 mapping data bases - Do not execute!
# #'
# #' @return A tibble with ensembl-transcript mapping
# #'
# get_ensembl_symbol_mapping <- function() {
#   # wget -O mapping_38.txt 'http://www.ensembl.org/biomart/martservice?query=  <Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6">  <Dataset name="hsapiens_gene_ensembl" interface="default"> <Attribute name="ensembl_transcript_id"/> <Attribute name="ensembl_gene_id"/><Attribute name="transcript"/> </Dataset> </Query>'
#   # wget -O mapping_37.txt 'http://grch37.ensembl.org/biomart/martservice?query=<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6"><Dataset name="hsapiens_gene_ensembl" interface="default"><Attribute name="ensembl_transcript_id"/><Attribute name="ensembl_gene_id"/><Attribute name="transcript"/></Dataset></Query>'
#   all =
#   	read_table2("~/third_party_sofware/ensembl_mapping/mapping_37.txt",
#               col_names = FALSE) %>%
#     setNames(c("ensembl_transcript_id", "ensembl_gene_id", "transcript")) %>%
#     mutate(hg = "hg37") %>%
#     bind_rows(
#       read_table2(
#         "~/third_party_sofware/ensembl_mapping/mapping_38.txt",
#         col_names = FALSE
#       ) %>%
#         setNames(
#           c("ensembl_transcript_id", "ensembl_gene_id", "transcript")
#         ) %>%
#         mutate(hg = "hg38")
#     ) %>%
#     drop_na()
#
#
#   bind_rows(
#   	# Gene annotation
#   	all %>%
# 		select(-ensembl_transcript_id) %>%
# 		group_by(ensembl_gene_id) %>%
# 		arrange(hg %>% desc()) %>%
# 		slice(1) %>%
# 		ungroup() %>%
# 		rename(ensembl_id = ensembl_gene_id),
#
# 		# Transcript annotation
# 		all %>%
# 		select(-ensembl_gene_id) %>%
# 		group_by(ensembl_transcript_id) %>%
# 		arrange(hg %>% desc()) %>%
# 		slice(1) %>%
# 		ungroup() %>%
# 		rename(ensembl_id = ensembl_transcript_id)
#   ) %>%
#
#   # Write to file and return
#   {
#     (.) %>% write_csv("~/third_party_sofware/ensembl_mapping/ensembl_symbol_mapping.csv")
#     (.)
#   }
# }

#' get_symbol_from_ensembl
#'
#' @description Get transcript column from ensembl gene id
#'
#' @param .data A tibble
#' @param .ensembl A column symbol. The column that is represents ensembl gene id
#'
#' @return A tibble with added annotation
#'
#'
get_symbol_from_ensembl <-
	function(.data, .ensembl) {
		.ensembl = enquo(.ensembl)

		.data %>%
			select(!!.ensembl) %>%
			distinct() %>%

			# Add name information
			dplyr::left_join(
				tidyBulk::ensembl_symbol_mapping %>%
					distinct(ensembl_id, transcript, hg) %>%
					dplyr::rename(!!.ensembl := ensembl_id) %>%
					rename(transcript = transcript),
				by = quo_name(.ensembl)
			)

	}

#' Add transcript column from ensembl gene id
#'
#' @param .data A tibble
#' @param .ensembl A column symbol. The column that is represents ensembl gene id
#'
#' @return A tibble with added annotation
#'
#'
add_symbol_from_ensembl <-
	function(.data, .ensembl) {
		# Comply with CRAN NOTES
		. = NULL

		.ensembl = enquo(.ensembl)

		# Add new symbols column
		.data %>%
			dplyr::left_join((.) %>%
											 	get_symbol_from_ensembl(!!.ensembl)) %>%

			# Attach attributes
			reattach_internals(.data)
	}

#' Perform linear equation system analysis through llsr
#'
#' @importFrom stats lsfit
#'
#' @param mix A data frame
#' @param reference A data frame
#'
#' @return A data frame
#'
#'
run_llsr = function(mix, reference) {
	# Get common markers
	markers = intersect(rownames(mix), rownames(reference))

	X <- (reference[markers, , drop = FALSE])
	Y <- (mix[markers, , drop = FALSE])

	results <- t(data.frame(lsfit(X, Y)$coefficients)[-1, , drop = FALSE])
	results[results < 0] <- 0
	results <- results / apply(results, 1, sum)
	rownames(results) = colnames(Y)

	results
}


#' Get cell type proportions from cibersort
#'
#' @import parallel
#' @import preprocessCore
#' @importFrom stats setNames
#' @importFrom rlang dots_list
#' @importFrom magrittr equals
#' @importFrom utils installed.packages
#' @importFrom utils install.packages
#'
#' @param .data A tibble
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param reference A data frame. The transcript/cell_type data frame of integer transcript abundance
#' @param method A character string. The method to be used. At the moment Cibersort (default) and llsr (linear least squares regression) are available.
#' @param ... Further parameters passed to the function Cibersort
#'
#' @return A tibble including additional columns
#'
#'
get_cell_type_proportions = function(.data,
																		 .sample = NULL,
																		 .transcript = NULL,
																		 .abundance = NULL,
																		 reference = X_cibersort,
																		 method = "cibersort",
																		 ...) {
	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
	.sample = col_names$.sample
	.transcript = col_names$.transcript
	.abundance = col_names$.abundance

	# Check if package is installed, otherwise install
	if ("e1071" %in% rownames(installed.packages()) == FALSE) {
		writeLines("Installing e1071 needed for Cibersort")
		install.packages("e1071", repos = "https://cloud.r-project.org")
	}

	# Check if package is installed, otherwise install
	if ("preprocessCore" %in% rownames(installed.packages()) == FALSE) {
		writeLines("Installing preprocessCore needed for Cibersort")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("preprocessCore")

	}

	# Load library which is optional for the whole package
	#library(preprocessCore)

	# Check if there are enough genes for the signature
	if ((.data %>%
			 pull(!!.transcript) %in% (reference %>% rownames)) %>%
			which %>%
			length %>%
			`<` (50))
		stop(
			"You have less than 50 genes that overlap the Cibersort signature. Please check again your input dataframe"
		)

	# Check if rownames exist
	if (reference %>% sapply(class) %in% c("numeric", "double", "integer") %>% `!` %>% any)
		stop("tidyBulk says: your reference has non-numeric/integer columns.")

	# Get the dots arguments
	dots_args = rlang::dots_list(...)

	.data %>%

		# Check if some transcripts are duplicated
		error_if_duplicated_genes(!!.sample,!!.transcript,!!.abundance) %>%

		# Prepare data frame
		distinct(!!.sample,!!.transcript,!!.abundance) %>%
		spread(!!.sample,!!.abundance) %>%
		data.frame(row.names = 1, check.names = FALSE) %>%

		# Run Cibersort or llsr through custom function, depending on method choice
		ifelse2_pipe(
			method %>% tolower %>% equals("cibersort"),
			method %>% tolower %>% equals("llsr"),

			# Execute do.call because I have to deal with ...
			~ do.call(my_CIBERSORT, list(Y = .x, X = reference) %>% c(dots_args)) %$%
				proportions %>%
				as_tibble(rownames = quo_name(.sample)) %>%
				select(-`P-value`,-Correlation,-RMSE),

			# Don't need to execute do.call
			~ .x %>%
				run_llsr(reference) %>%
				as_tibble(rownames = quo_name(.sample)),

			~ stop(
				"tidyBulk syas: please choose between cibersort and llsr methods"
			)
		)	 %>%

		# Parse results and return
		setNames(c(
			quo_name(.sample),
			(.) %>% select(-1) %>% colnames() %>% sprintf("%s: %s", method, .)
		)) %>%
		#%>%
		#gather(`Cell type`, proportion,-!!.sample) %>%

		# Attach attributes
		reattach_internals(.data)


}

#' Add cell type proportions from cibersort
#'
#' @import parallel
#'
#'
#' @param .data A tibble
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param reference A data frame. The transcript/cell_type data frame of integer transcript abundance
#' @param method A character string. The method to be used. At the moment Cibersort (default) and llsr (linear least squares regression) are available.
#' @param ... Further parameters passed to the function Cibersort
#'
#' @return A tibble including additional columns
#'
#'
add_cell_type_proportions = function(.data,
																		 .sample = NULL,
																		 .transcript = NULL,
																		 .abundance = NULL,
																		 reference = X_cibersort,
																		 method = "cibersort",
																		 ...) {
	# Comply with CRAN NOTES
	. = NULL

	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
	.sample = col_names$.sample
	.transcript = col_names$.transcript
	.abundance = col_names$.abundance

	.data %>%

		# Add new annotation
		dplyr::left_join(
			(.) %>%
				get_cell_type_proportions(
					.sample = !!.sample,
					.transcript = !!.transcript,
					.abundance = !!.abundance,
					reference = reference,
					method = method,
					...
				),
			by = quo_name(.sample)
		) %>%

		# Attach attributes
		reattach_internals(.data)

}

#' Get adjusted count for some batch effect
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr set_colnames
#' @importFrom stats model.matrix
#' @importFrom stats as.formula
#' @importFrom utils installed.packages
#' @importFrom utils install.packages
#'
#' @param .data A tibble
#' @param .formula a formula with no response variable, of the kind ~ factor_of_intrest + batch
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#' @param ... Further parameters passed to the function sva::ComBat
#'
#' @return A tibble with adjusted counts
#'
#'
get_adjusted_counts_for_unwanted_variation_bulk <- function(.data,
																														.formula,
																														.sample = NULL,
																														.transcript = NULL,
																														.abundance = NULL,
																														log_transform = TRUE,
																														...) {
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

	# Check that .formula includes at least two covariates
	if (parse_formula(.formula) %>% length %>% `<` (2))
		stop(
			"The .formula must contain two covariates, the first being the factor of interest, the second being the factor of unwanted variation"
		)

	# Check that .formula includes no more than two covariates at the moment
	if (parse_formula(.formula) %>% length %>% `>` (3))
		warning("Only the second covariate in the .formula is adjusted for, at the moment")

	# New column name
	value_adjusted = as.symbol(sprintf("%s adjusted",  quo_name(.abundance)))

	# Stop is any counts are NAs
	.data %>% error_if_counts_is_na(!!.abundance)

	df_for_combat <-
		.data %>%

		# Filter low counts
		filter_abundant(!!.sample, !!.transcript, !!.abundance) %>%
		{
			# Give warning of filtering
			message(
				"Combat is applied excluding `filter out low counts`, as it performs on non sparse matrices. Therefore NAs will be used for those lowly abundant transcripts "
			)
			(.)
		} %>%

		select(!!.transcript,
					 !!.sample,
					 !!.abundance,
					 one_of(parse_formula(.formula))) %>%
		distinct() %>%

		# Check if logtansform is needed
		ifelse_pipe(log_transform,
								~ .x %>% dplyr::mutate(!!.abundance := !!.abundance %>% `+`(1) %>%  log()))


	# Create design matrix
	design =
		model.matrix(
			object = as.formula("~" %>% paste0(parse_formula(.formula)[1])),
			# get first argument of the .formula
			data = df_for_combat %>% select(!!.sample, one_of(parse_formula(.formula))) %>% distinct %>% arrange(!!.sample)
		) %>%
		set_colnames(c("(Intercept)", parse_formula(.formula)[1]))

	# Check if package is installed, otherwise install
	if ("sva" %in% rownames(installed.packages()) == FALSE) {
		writeLines("Installing sva - Combat needed for adjustment for unwanted variation")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("sva")
	}

	my_batch =
		df_for_combat %>%
		distinct(!!.sample,!!as.symbol(parse_formula(.formula)[2])) %>%
		arrange(!!.sample) %>%
		pull(2)

	mat = df_for_combat %>%
		# Select relevant info
		distinct(!!.transcript,!!.sample,!!.abundance) %>%

		# Stop any column is not if not numeric or integer
		ifelse_pipe(
			(.) %>% select(!!.abundance) %>% summarise_all(class) %>% `%in%`(c("numeric", "integer")) %>% `!`() %>% any(),
			~ stop(".abundance must be numerical or integer")
		) %>%

		spread(!!.sample,!!.abundance) %>%
		as_matrix(rownames = !!.transcript,
							do_check = FALSE)
	mat %>%

		# Run combat
		sva::ComBat(batch = my_batch,
								mod = design,
								prior.plots = FALSE,
								...) %>%

		as_tibble(rownames = quo_name(.transcript)) %>%
		gather(!!.sample,!!.abundance,-!!.transcript) %>%

		# ReverseLog transform if tranformed in the first place
		ifelse_pipe(
			log_transform,
			~ .x %>%
				dplyr::mutate(!!.abundance := !!.abundance %>% exp %>% `-`(1)) %>%
				dplyr::mutate(!!.abundance := ifelse(!!.abundance < 0, 0,!!.abundance)) %>%
				dplyr::mutate(!!.abundance := !!.abundance %>% as.integer)
		) %>%

		# Reset column names
		dplyr::rename(!!value_adjusted := !!.abundance)  %>%

		# # Add filtering info
		# right_join(
		# 	df_for_combat %>%
		# 		distinct(!!.transcript,!!.sample,
		# 						 `filter out low counts`),
		# 	by = c(quo_name(.transcript), quo_name(.sample))
		# )%>%

		# Attach attributes
		reattach_internals(.data)
}

#' Add adjusted count for some batch effect
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr set_colnames
#'
#' @param .data A tibble
#' @param .formula a formula with no response variable, of the kind ~ factor_of_intrest + batch
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#' @param ... Further parameters passed to the function sva::ComBat
#'
#' @return A tibble with adjusted counts
#'
#'
add_adjusted_counts_for_unwanted_variation_bulk <- function(.data,
																														.formula,
																														.sample = NULL,
																														.transcript = NULL,
																														.abundance = NULL,
																														log_transform = TRUE,
																														...) {
	# Comply with CRAN NOTES
	. = NULL

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

	.data %>%

		# Add adjsted column
		dplyr::left_join(
			(.) %>%
				get_adjusted_counts_for_unwanted_variation_bulk(
					.formula,
					.sample = !!.sample,
					.transcript = !!.transcript,
					.abundance = !!.abundance,
					log_transform = log_transform,
					...
				) ,
			by = c(quo_name(.transcript), quo_name(.sample))
		) %>%

		# Attach attributes
		reattach_internals(.data)

}

#' Identify variable genes for dimensionality reduction
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
filter_variable_transcripts = function(.data,
																			 .sample = NULL,
																			 .transcript = NULL,
																			 .abundance = NULL,
																			 top = 500,
																			 log_transform = TRUE) {
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

	# Manage Inf
	top = min(top, .data %>% distinct(!!.transcript) %>% nrow)

	writeLines(sprintf("Getting the %s most variable genes", top))

	x =
		.data %>%
		distinct(!!.sample,!!.transcript,!!.abundance) %>%

		# Check if logtansform is needed
		ifelse_pipe(log_transform,
								~ .x %>% dplyr::mutate(!!.abundance := !!.abundance %>% `+`(1) %>%  log())) %>%

		spread(!!.sample,!!.abundance) %>%
		as_matrix(rownames = quo_name(.transcript))

	s <- rowMeans((x - rowMeans(x)) ^ 2)
	o <- order(s, decreasing = TRUE)
	x <- x[o[1L:top], , drop = FALSE]
	variable_trancripts = rownames(x)

	.data %>% filter(!!.transcript %in% variable_trancripts)
}

#'ttBulk_to_SummarizedExperiment
#'
#' @importFrom utils data
#'
#' @param .data A tibble
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#'
#' @return A SummarizedExperiment
#'
ttBulk_to_SummarizedExperiment = function(.data,
																					.sample = NULL,
																					.transcript = NULL,
																					.abundance = NULL) {
	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
	.sample = col_names$.sample
	.transcript = col_names$.transcript
	.abundance = col_names$.abundance

	# Check if package is installed, otherwise install
	if ("SummarizedExperiment" %in% rownames(installed.packages()) == FALSE) {
		writeLines("Installing SummarizedExperiment")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("SummarizedExperiment")
	}
	if ("S4Vectors" %in% rownames(installed.packages()) == FALSE) {
		writeLines("Installing S4Vectors")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("S4Vectors")
	}
	# If present get the scaled abundance
	.abundance_scaled =
		.data %>%
		ifelse_pipe(
			".abundance_scaled" %in% ((.) %>% get_tt_columns() %>% names) &&
				# .data %>% get_tt_columns() %$% .abundance_scaled %>% is.null %>% `!` &&
				quo_name((.) %>% get_tt_columns() %$% .abundance_scaled) %in% ((.) %>% colnames),
			~ .x %>% get_tt_columns() %$% .abundance_scaled,
			~ NULL
		)

	# Get which columns are sample wise and which are feature wise
	col_direction = get_x_y_annotation_columns(.data,
																						 !!.sample,
																						 !!.transcript,
																						 !!.abundance,
																						 !!.abundance_scaled)
	sample_cols = col_direction$horizontal_cols
	feature_cols = col_direction$vertical_cols
	counts_cols = col_direction$counts_cols

	colData = .data %>% select(!!.sample, sample_cols) %>% distinct %>% arrange(!!.sample) %>% {
		S4Vectors::DataFrame((.) %>% select(-!!.sample),
												 row.names = (.) %>% pull(!!.sample))
	}

	rowData = .data %>% select(!!.transcript, feature_cols) %>% distinct %>% arrange(!!.transcript) %>% {
		S4Vectors::DataFrame((.) %>% select(-!!.transcript),
												 row.names = (.) %>% pull(!!.transcript))
	}

	my_assays =
		.data %>%
		select(!!.sample,
					 !!.transcript,
					 !!.abundance,
					 !!.abundance_scaled,
					 counts_cols) %>%
		distinct() %>%
		gather(`assay`, .a,-!!.transcript,-!!.sample) %>%
		nest(`data` = -`assay`) %>%
		mutate(`data` = `data` %>%  map(
			~ .x %>% spread(!!.sample, .a) %>% as_matrix(rownames = quo_name(.transcript))
		))

	# Build the object
	SummarizedExperiment::SummarizedExperiment(
		assays = my_assays %>% pull(`data`) %>% setNames(my_assays$assay),
		rowData = rowData,
		colData = colData
	)

}

#' Get matrix from tibble
#'
#' @import dplyr
#' @import tidyr
#' @importFrom magrittr set_rownames
#' @importFrom rlang quo_is_null
#'
#' @param tbl A tibble
#' @param rownames A character string of the rownames
#' @param do_check A boolean
#'
#' @return A matrix
#'
#' @examples
#'
#' as_matrix(head(dplyr::select(tidyBulk::counts_mini, transcript, count)), rownames=transcript)
#'
#' @export
as_matrix <- function(tbl,
											rownames = NULL,
											do_check = TRUE) {
	rownames = enquo(rownames)
	tbl %>%

		# Through warning if data frame is not numerical beside the rownames column (if present)
		ifelse_pipe(
			do_check &&
				tbl %>%
				# If rownames defined eliminate it from the data frame
				ifelse_pipe(!quo_is_null(rownames), ~ .x[,-1], ~ .x) %>%
				dplyr::summarise_all(class) %>%
				tidyr::gather(variable, class) %>%
				pull(class) %>%
				unique() %>%
				`%in%`(c("numeric", "integer")) %>% `!`() %>% any(),
			~ {
				warning("to_matrix says: there are NON-numerical columns, the matrix will NOT be numerical")
				.x
			}
		) %>%
		as.data.frame() %>%

		# Deal with rownames column if present
		ifelse_pipe(
			!quo_is_null(rownames),
			~ .x %>%
				magrittr::set_rownames(tbl %>% pull(!!rownames)) %>%
				select(-1)
		) %>%

		# Convert to matrix
		as.matrix()
}
