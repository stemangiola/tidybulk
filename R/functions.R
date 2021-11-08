#' Create tt object from tibble
#'
#' @keywords internal
#' @noRd
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
		memorise_methods_used("tidyverse") %>%

		# Add class
		add_class("tt") %>%
		add_class("tidybulk")
}



#' Convert bam/sam files to a tidy gene transcript counts data frame
#'
#' @keywords internal
#' @noRd
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
			when(
				"annot.ext" %in% (rlang::dots_list(...) %>% names) %>% not() ~ {
					dge <- (.)
					dge$genes$symbol <-
						AnnotationDbi::mapIds(
							org.Hs.eg.db::org.Hs.eg.db,
							keys = as.character(dge$genes$GeneID),
							column = "SYMBOL",
							keytype = "ENTREZID",
							multiVals = "first"
						)

					dge
				},
				~ (.)
			) %>%

			# Anonymous function
			# input: annotated edgeR::DGEList object
			# output: tibble
			{
				reduce(
					list(

						# Counts
						(.) %$%
							counts %>%
							as_tibble(rownames = "GeneID") %>%
							gather(sample, count,-GeneID),

						# Genes
						(.) %$%
							genes %>%
							select(
								suppressWarnings(
									one_of("GeneID", "symbol")
								)
								) %>%
							as_tibble() %>%
							mutate(GeneID = GeneID %>% as.character()),

						# Sample
						(.) %$%
							samples %>%
							as_tibble()
					),
					dplyr::left_join
				) %>%
					rename(transcript = GeneID)
			} %>%

			# Add tt_columns attribute

			tidybulk_to_SummarizedExperiment(sample, transcript, count)
			#
			# add_tt_columns(sample,symbol,count) %>%
			# memorise_methods_used("featurecounts") %>%
			#
			# # Add class
			# add_class("tt") %>%
			# add_class("tidybulk")
	}

#' Calculate the norm factor with calcNormFactor from limma
#'
#' @keywords internal
#' @noRd
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang :=
#' @importFrom stats setNames
#'
#' @param .data A tibble
#' @param reference A reference matrix, not sure if used anymore
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param method A string character. The scaling method passed to the backend function (i.e., edgeR::calcNormFactors; "TMM","TMMwsp","RLE","upperquartile")
#'
#'
#' @return A list including the filtered data frame and the normalization factors
add_scaled_counts_bulk.calcNormFactor <- function(.data,
                                                  reference = NULL,
                                                  .sample = `sample`,
                                                  .transcript = `transcript`,
                                                  .abundance = `count`,
                                                  method) {
 	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	error_if_log_transformed(.data,!!.abundance)

	# Get data frame for the highly transcribed transcripts
	df.filt <-
		.data %>%
		# dplyr::filter(!(!!.transcript %in% gene_to_exclude)) %>%
		droplevels() %>%
		select(!!.sample, !!.transcript, !!.abundance)

	df.filt.spread =
	  df.filt %>%
	  tidyr::spread(!!.sample,!!.abundance) %>%
	  tidyr::drop_na() %>%
	  dplyr::select(-!!.transcript)

	# If not enough genes, warning
	if(nrow(df.filt.spread)<100) warning(warning_for_scaling_with_few_genes)

	# scaled data set
	nf =
		tibble::tibble(
			# Sample factor
			sample = factor(levels(df.filt %>% pull(!!.sample))),

			# scaled data frame
			nf = edgeR::calcNormFactors(
			  df.filt.spread,
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
  list(
    # gene_to_exclude = gene_to_exclude,
    nf = nf
  ) %>%

    # Attach attributes
    reattach_internals(.data)
}

#' Get a tibble with scaled counts using TMM
#'
#' @keywords internal
#' @noRd
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr equals
#' @importFrom rlang :=
#' @importFrom stats median
#' @importFrom utils install.packages
#'
#' @param .data A tibble
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param method A character string. The scaling method passed to the backend function (i.e., edgeR::calcNormFactors; "TMM","TMMwsp","RLE","upperquartile")
#' @param reference_sample A character string. The name of the reference sample. If NULL the sample with highest total read count will be selected as reference.
#'
#'
#' @return A tibble including additional columns
#'
#'
get_scaled_counts_bulk <- function(.data,
																	 .sample = NULL,
																	 .transcript = NULL,
																	 .abundance = NULL,
																	 method = "TMM",
																	 reference_sample = NULL) {
	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	# Check if package is installed, otherwise install
	if (find.package("edgeR", quiet = TRUE) %>% length %>% equals(0)) {
		message("tidybulk says: Installing edgeR needed for analyses")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("edgeR", ask = FALSE)
	}

	# Reformat input data set
	df <-
		.data %>%

		# Rename
		dplyr::select(!!.sample,!!.transcript,!!.abundance) %>%

		# Set samples and genes as factors
		dplyr::mutate(!!.sample := factor(!!.sample),!!.transcript := factor(!!.transcript))


	# Get reference
	reference <-
		reference_sample %>%
		when(
			!is.null(.) ~ (.),

			# If not specified take most abundance sample
			df %>%
				group_by(!!.sample) %>%
				summarise(sum = median(!!.abundance)) %>%
				mutate(med = max(sum)) %>%
				mutate(diff = abs(sum - med)) %>%
				arrange(diff) %>%
				head(n = 1) %>%
				pull(!!.sample) %>%
				as.character()
		)

	# Communicate the reference if chosen by default
  if(is.null(reference_sample)) message(sprintf("tidybulk says: the sample with largest library size %s was chosen as reference for scaling", reference))

	nf_obj <-
		add_scaled_counts_bulk.calcNormFactor(
			df,
			reference,
			.sample = !!.sample,
			.transcript = !!.transcript,
			.abundance = !!.abundance,
			method
		)

	# Calculate normalization factors
	nf_obj$nf %>%
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

					 	# Put everything to the reference sample scale
					 	((.) %>% filter(!!.sample == reference) %>% pull(tot))) %>%

		# I have correct the strange behaviour of edgeR of reference
		# sample not being 1
		# I HAD TO COMMENT BECAUSE TEST FAILING
		# {
		# 	mult_ref = (.) %>%  filter(!!.sample == reference) %>% pull(multiplier)
		# 	(.) %>%  mutate(
		# 		multiplier =
		# 			multiplier /mult_ref
		# 	)
		# } %>%

		dplyr::select(-tot,-tot_filt) %>%
		dplyr::rename(TMM = nf) %>%

		# # Attach internals
		# add_tt_columns(!!.sample,!!.transcript,!!.abundance) %>%

		# Add methods
		memorise_methods_used(c("edger", "tmm"))

}



#' Get differential transcription information to a tibble using edgeR.
#'
#' @keywords internal
#' @noRd
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr set_colnames
#' @importFrom stats model.matrix
#' @importFrom utils install.packages
#' @importFrom purrr when
#' @importFrom rlang inform
#'
#'
#' @param .data A tibble
#' @param .formula a formula with no response variable, referring only to numeric variables
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param .contrasts A character vector. See edgeR makeContrasts specification for the parameter `contrasts`. If contrasts are not present the first covariate is the one the model is tested against (e.g., ~ factor_of_interest)
#' @param method A string character. Either "edgeR_quasi_likelihood" (i.e., QLF), "edgeR_likelihood_ratio" (i.e., LRT)
#' @param test_above_log2_fold_change A positive real value. This works for edgeR and limma_voom methods. It uses the `treat` function, which tests that the difference in abundance is bigger than this threshold rather than zero \url{https://pubmed.ncbi.nlm.nih.gov/19176553}.
#' @param scaling_method A character string. The scaling method passed to the backend function (i.e., edgeR::calcNormFactors; "TMM","TMMwsp","RLE","upperquartile")
#' @param omit_contrast_in_colnames If just one contrast is specified you can choose to omit the contrast label in the colnames.
#' @param .sample_total_read_count
#'
#' @return A tibble with edgeR results
#'
get_differential_transcript_abundance_bulk <- function(.data,
																											 .formula,
																											 .sample = NULL,
																											 .transcript = NULL,
																											 .abundance = NULL,
																											 .contrasts = NULL,
																											 method = "edgeR_quasi_likelihood",
																											 test_above_log2_fold_change = NULL,
																											 scaling_method = "TMM",
																											 omit_contrast_in_colnames = FALSE,
																											 prefix = "",
																											 .sample_total_read_count = NULL) {
	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	.sample_total_read_count = enquo(.sample_total_read_count)

	# Check if omit_contrast_in_colnames is correctly setup
	if(omit_contrast_in_colnames & length(.contrasts) > 1){
		warning("tidybulk says: you can omit contrasts in column names only when maximum one contrast is present")
		omit_contrast_in_colnames = FALSE
	}

	# distinct_at is not released yet for dplyr, thus we have to use this trick
	df_for_edgeR <- .data %>%

		# Prepare the data frame
		select(!!.transcript,
					 !!.sample,
					 !!.abundance,
					 one_of(parse_formula(.formula))) %>%
		distinct() %>%

		# drop factors as it can affect design matrix
		droplevels()
	#%>%

		# # Check if data rectangular
		# ifelse2_pipe(
		# 	(.) %>% check_if_data_rectangular(!!.sample,!!.transcript,!!.abundance) %>% not() & fill_missing_values,
		# 	(.) %>% check_if_data_rectangular(!!.sample,!!.transcript,!!.abundance) %>% not() & !fill_missing_values,
		# 	~ .x %>% fill_NA_using_formula(.formula,!!.sample, !!.transcript, !!.abundance),
		# 	~ .x %>% eliminate_sparse_transcripts(!!.transcript)
		# )


	# # Check if at least two samples for each group
	# if (
	# 	# If I have some discrete covariates
	# 	df_for_edgeR %>%
	# 	select(one_of(parse_formula(.formula))) %>%
	# 	select_if(function(col)
	# 		is.character(col) | is.factor(col) | is.logical(col)) %>%
	# 	ncol %>% gt(0) &
	#
	# 	# If I have at least 2 samples per group
	# 	df_for_edgeR %>%
	# 	select(!!.sample, one_of(parse_formula(.formula))) %>%
	# 	select_if(function(col) !is.numeric(col) & !is.integer(col) & !is.double(col) ) %>%
	# 	distinct %>%
	# 	group_by_at(vars(-!!.sample)) %>%
	# 	count() %>%
	# 	ungroup() %>%
	# 	{
	# 		(.) %>% nrow %>% st(2) |
	# 		(.) %>% distinct(n) %>%	pull(n) %>%	min %>% st(2)
	# 	}
	# )
	# message("tidybulk says: Just so you know. You have less than two replicates for each factorial combination")

	# Create design matrix
	design =
		model.matrix(
			object = .formula,
			data = df_for_edgeR %>% select(!!.sample, one_of(parse_formula(.formula))) %>% distinct %>% arrange(!!.sample)
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
		message("tidybulk says: Installing edgeR needed for differential transcript abundance analyses")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("edgeR", ask = FALSE)
	}



	edgeR_object =
		df_for_edgeR %>%
		select(!!.transcript,!!.sample,!!.abundance) %>%
		spread(!!.sample,!!.abundance) %>%
		as_matrix(rownames = !!.transcript) %>%

		edgeR::DGEList(counts = .) %>%

		# Override lib.size if imposed
		# This is useful in case you are analysing a small amount of genes,
		# for which the overall lib.size cannot be calculated
		when(
			!quo_is_null(.sample_total_read_count) ~ {

				# New library size dataset
				new_lib_size = .data %>% pivot_sample(!!.sample) %>% select(!!.sample, !!.sample_total_read_count)

				x = (.)
				x$samples$lib.size =
					new_lib_size %>%
					slice(match(rownames(x$samples), !!.sample)) %>%
					pull(!!.sample_total_read_count)

				x
			},
			~ (.)
		) %>%

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
				as_tibble(rownames = quo_name(.transcript)) %>%

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
							as_tibble(rownames = quo_name(.transcript)) %>%
							mutate(constrast = colnames(my_contrasts)[.x])
							# %>%
							#
							# # Mark DE genes
							# mutate(significant = FDR < significance_threshold)
					) %>%
					pivot_wider(values_from = -c(!!.transcript, constrast),
											names_from = constrast, names_sep = "___")
			}
		)	 %>%

		# Attach prefix
		setNames(c(
			colnames(.)[1],
			sprintf("%s%s", prefix, colnames(.)[2:ncol(.)])
		)) %>%

		# Attach attributes
		reattach_internals(.data) %>%

	    # select method
			when(
				tolower(method) ==  "edger_likelihood_ratio" ~ (.) %>% memorise_methods_used(c("edger", "edgeR_likelihood_ratio")),
				tolower(method) ==  "edger_quasi_likelihood" ~ (.) %>% memorise_methods_used(c("edger", "edgeR_quasi_likelihood")),
				tolower(method) ==  "edger_robust_likelihood_ratio" ~ (.) %>% memorise_methods_used(c("edger", "edger_robust_likelihood_ratio"))
			)	%>%
        	when(
        	    !is.null(test_above_log2_fold_change) ~ (.) %>% memorise_methods_used("treat"),
        	    ~ (.)
        	) %>%

		# Add raw object
		attach_to_internals(edgeR_object, "edgeR") %>%
		# Communicate the attribute added
		{

		  rlang::inform("tidybulk says: to access the raw results (fitted GLM) do `attr(..., \"internals\")$edgeR`", .frequency_id = "Access DE results edgeR",  .frequency = "once")

			(.)
		}
}


#' Get differential transcription information to a tibble using voom.
#'
#' @keywords internal
#' @noRd
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr set_colnames
#' @importFrom stats model.matrix
#' @importFrom utils install.packages
#' @importFrom purrr when
#' @importFrom rlang inform
#'
#' @param .data A tibble
#' @param .formula a formula with no response variable, referring only to numeric variables
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param .contrasts A character vector. See voom makeContrasts specification for the parameter `contrasts`. If contrasts are not present the first covariate is the one the model is tested against (e.g., ~ factor_of_interest)
#' @param method A string character. Either "limma_voom", "limma_voom_sample_weights"
#' @param scaling_method A character string. The scaling method passed to the backend function (i.e., edgeR::calcNormFactors; "TMM","TMMwsp","RLE","upperquartile")
#' @param omit_contrast_in_colnames If just one contrast is specified you can choose to omit the contrast label in the colnames.
#'
#' @return A tibble with voom results
#'
get_differential_transcript_abundance_bulk_voom <- function(.data,
																											 .formula,
																											 .sample = NULL,
																											 .transcript = NULL,
																											 .abundance = NULL,
																											 .contrasts = NULL,
                                                                                                             method = NULL,
																											 test_above_log2_fold_change = NULL,
																											 scaling_method = "TMM",
																											 omit_contrast_in_colnames = FALSE,
																											 prefix = "") {
	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	# Check if omit_contrast_in_colnames is correctly setup
	if(omit_contrast_in_colnames & length(.contrasts) > 1){
		warning("tidybulk says: you can omit contrasts in column names only when maximum one contrast is present")
		omit_contrast_in_colnames = FALSE
	}

	# distinct_at is not released yet for dplyr, thus we have to use this trick
	df_for_voom <- .data %>%

		# Prepare the data frame
		select(!!.transcript,
					 !!.sample,
					 !!.abundance,
					 one_of(parse_formula(.formula))) %>%
		distinct() %>%

		# drop factors as it can affect design matrix
		droplevels()


	# Create design matrix
	design =
		model.matrix(
			object = .formula,
			data = df_for_voom %>% select(!!.sample, one_of(parse_formula(.formula))) %>% distinct %>% arrange(!!.sample)
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
		message("tidybulk says: Installing limma needed for differential transcript abundance analyses")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("limma", ask = FALSE)
	}

	voom_object =
		df_for_voom %>%
		select(!!.transcript,!!.sample,!!.abundance) %>%
		spread(!!.sample,!!.abundance) %>%
		as_matrix(rownames = !!.transcript) %>%

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
				as_tibble(rownames = quo_name(.transcript)) %>%

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
							as_tibble(rownames = quo_name(.transcript)) %>%
							mutate(constrast = colnames(my_contrasts)[.x])
							# %>%
							#
							# # Mark DE genes
							# mutate(significant = adj.P.Val < significance_threshold)
					) %>%
					pivot_wider(values_from = -c(!!.transcript, constrast),
											names_from = constrast, names_sep = "___")
			}
		) %>%

		# Attach prefix
		setNames(c(
			colnames(.)[1],
			sprintf("%s%s", prefix, colnames(.)[2:ncol(.)])
		)) %>%

		# Attach attributes
		reattach_internals(.data) %>%

	    # select method
	    when(
			tolower(method) == "limma_voom" ~ (.) %>% memorise_methods_used("voom"),
			tolower(method) == "limma_voom_sample_weights" ~ (.) %>% memorise_methods_used("voom_sample_weights")
		) %>%
	    when(
			!is.null(test_above_log2_fold_change) ~ (.) %>% memorise_methods_used("treat"),
			~ (.)
		) %>%

		# Add raw object
		attach_to_internals(voom_object, "voom") %>%
		# Communicate the attribute added
		{
		  rlang::inform("tidybulk says: to access the raw results (fitted GLM) do `attr(..., \"internals\")$voom`", .frequency_id = "Access DE results voom",  .frequency = "once")

			(.)
		}
}


#' Get differential transcription information to a tibble using DESeq2
#'
#' @keywords internal
#' @noRd
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr set_colnames
#' @importFrom stats model.matrix
#' @importFrom utils install.packages
#' @importFrom purrr when
#' @importFrom rlang inform
#'
#' @param .data A tibble
#' @param .formula a formula with no response variable, referring only to numeric variables
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param .contrasts A character vector. See edgeR makeContrasts specification for the parameter `contrasts`. If contrasts are not present the first covariate is the one the model is tested against (e.g., ~ factor_of_interest)
#' @param method A string character. Either "edgeR_quasi_likelihood" (i.e., QLF), "edgeR_likelihood_ratio" (i.e., LRT)
#' @param scaling_method A character string. The scaling method passed to the backend function (i.e., edgeR::calcNormFactors; "TMM","TMMwsp","RLE","upperquartile")
#' @param omit_contrast_in_colnames If just one contrast is specified you can choose to omit the contrast label in the colnames.
#' @param ... Additional arguments for DESeq2
#'
#' @return A tibble with edgeR results
#'
get_differential_transcript_abundance_deseq2 <- function(.data,
																											 .formula,
																											 .sample = NULL,
																											 .transcript = NULL,
																											 .abundance = NULL,
																											 .contrasts = NULL,
																											 method = "edgeR_quasi_likelihood",
																											 scaling_method = "TMM",
																											 omit_contrast_in_colnames = FALSE,
																											 prefix = "",
																											 ...) {

	# Check if contrasts are of the same form
	if(
		.contrasts %>% is.null %>% not() &
		.contrasts %>% class %>% equals("list") %>% not()
	)
		stop("tidybulk says: for DESeq2 the list of constrasts should be given in the form list(c(\"condition_column\",\"condition1\",\"condition2\")) i.e. list(c(\"genotype\",\"knockout\",\"wildtype\"))")

	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	# Check if omit_contrast_in_colnames is correctly setup
	if(omit_contrast_in_colnames & length(.contrasts) > 1){
		warning("tidybulk says: you can omit contrasts in column names only when maximum one contrast is present")
		omit_contrast_in_colnames = FALSE
	}

	if (find.package("acepack", quiet = TRUE) %>% length %>% equals(0)) {
		message("tidybulk says: Installing acepack needed for analyses")
		install.packages("acepack", repos = "https://cloud.r-project.org")
	}

	# Check if package is installed, otherwise install
	if (find.package("DESeq2", quiet = TRUE) %>% length %>% equals(0)) {
		message("tidybulk says: Installing DESeq2 needed for differential transcript abundance analyses")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("DESeq2", ask = FALSE)
	}

	# # Print the design column names in case I want contrasts
	# message(
	# 	sprintf(
	# 		"tidybulk says: The design column names are \"%s\"",
	# 		design %>% colnames %>% paste(collapse = ", ")
	# 	)
	# )

	my_contrasts = .contrasts


	deseq2_object =
		.data %>%

		# Prepare the data frame
		select(!!.transcript,
					 !!.sample,
					 !!.abundance,
					 one_of(parse_formula(.formula))) %>%
		distinct() %>%

		# drop factors as it can affect design matrix
		droplevels() %>%

		# Needed for DESeq2
		mutate(!!.abundance := as.integer(!!.abundance)) %>%
		rename(counts = !!.abundance) %>%
		mutate_if(is.logical , as.factor) %>%
		mutate_if(is.character , as.factor) %>%

		# Filter
		tidybulk_to_SummarizedExperiment(!!.sample, !!.transcript, counts) %>%

		# DESeq2
		DESeq2::DESeqDataSet( design = .formula) %>%
		DESeq2::DESeq(...)

	# Read ft object
	deseq2_object %>%

		# If I have multiple .contrasts merge the results
		when(

			# Simple comparison continuous
			(my_contrasts %>% is.null ) &
				(deseq2_object@colData[,parse_formula(.formula)[1]] %>%
				 	class %in% c("numeric", "integer", "double")) 	~
				(.) %>%
				DESeq2::results() %>%
				as_tibble(rownames = quo_name(.transcript)),

			# Simple comparison discrete
			my_contrasts %>% is.null 	~
				(.) %>%
				DESeq2::results(contrast = c(
					parse_formula(.formula)[1],
					deseq2_object@colData[,parse_formula(.formula)[1]] %>% levels %>% .[2],
					deseq2_object@colData[,parse_formula(.formula)[1]] %>% levels %>% .[1]
				)) %>%
				as_tibble(rownames = quo_name(.transcript)),

			# Simple comparison discrete
			my_contrasts %>% is.null %>% not() & omit_contrast_in_colnames	~
				(.) %>%
				DESeq2::results(contrast = my_contrasts[[1]])%>%
				as_tibble(rownames = quo_name(.transcript)),

			# Multiple comparisons NOT USED AT THE MOMENT
			~ {
				deseq2_obj = (.)

				1:length(my_contrasts) %>%
					map_dfr(
						~ 	deseq2_obj %>%

							# select method
							DESeq2::results(contrast = my_contrasts[[.x]])	%>%

							# Convert to tibble
							as_tibble(rownames = quo_name(.transcript)) %>%
							mutate(constrast = sprintf("%s %s-%s", my_contrasts[[.x]][1], my_contrasts[[.x]][2], my_contrasts[[.x]][3]) )

					) %>%
					pivot_wider(values_from = -c(!!.transcript, constrast),
											names_from = constrast, names_sep = "___")
			}
		)	 %>%

		# Attach prefix
		setNames(c(
			colnames(.)[1],
			sprintf("%s%s", prefix, colnames(.)[2:ncol(.)])
		)) %>%

		# Attach attributes
		reattach_internals(.data) %>%
		memorise_methods_used("deseq2") %>%

		# Add raw object
		attach_to_internals(deseq2_object, "DESeq2") %>%

		# Communicate the attribute added
		{

		  rlang::inform("tidybulk says: to access the raw results (fitted GLM) do `attr(..., \"internals\")$DESeq2`", .frequency_id = "Access DE results deseq2",  .frequency = "once")

			(.)
		}
}

#' Get differential composition information to a tibble using edgeR.
#'
#' @keywords internal
#' @noRd
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr set_colnames
#' @importFrom stats model.matrix
#' @importFrom utils install.packages
#' @importFrom purrr when
#' @importFrom purrr map_lgl
#' @importFrom stringr str_replace
#' @importFrom stringr str_split
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_remove
#'
#'
#' @param .data A tibble
#' @param .formula a formula with no response variable, referring only to numeric variables
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param method A string character. Either "edgeR_quasi_likelihood" (i.e., QLF), "edgeR_likelihood_ratio" (i.e., LRT)
#' @param reference A data frame. The transcript/cell_type data frame of integer transcript abundance
#' @param significance_threshold A real between 0 and 1
#'
#' @return A tibble with edgeR results
#'
test_differential_cellularity_ <- function(.data,
																					 .formula,
																					 .sample = NULL,
																					 .transcript = NULL,
																					 .abundance = NULL,
																					 method = "cibersort",
																					 reference = NULL,
																					 significance_threshold = 0.05,
																					 ...
) {

	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)


	if (find.package("broom", quiet = TRUE) %>% length %>% equals(0)) {
		message("tidybulk says: Installing broom needed for analyses")
		install.packages("broom", repos = "https://cloud.r-project.org")
	}

	deconvoluted =
		.data %>%

		# Deconvolution
		deconvolve_cellularity(
			!!.sample, !!.transcript, !!.abundance,
			method=method,
			prefix = sprintf("%s:", method),
			reference = reference,
			action="get",
			...
		)

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
					when(grepl("Surv", .my_formula) ~ (.) %>% memorise_methods_used(c("survival", "boot"), object_containing_methods = .data),
							 ~ (.))

			}) %>%

		# Eliminate prefix
		mutate(.cell_type = str_remove(.cell_type, sprintf("%s:", method)))


}


#' Get differential composition information to a tibble using edgeR.
#'
#' @keywords internal
#' @noRd
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr set_colnames
#' @importFrom stats model.matrix
#' @importFrom utils install.packages
#' @importFrom purrr when
#' @importFrom purrr map_lgl
#' @importFrom stringr str_replace
#' @importFrom stringr str_split
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_remove
#'
#'
#' @param .data A tibble
#' @param .formula a formula with no response variable, referring only to numeric variables
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param method A string character. Either "edgeR_quasi_likelihood" (i.e., QLF), "edgeR_likelihood_ratio" (i.e., LRT)
#' @param reference A data frame. The transcript/cell_type data frame of integer transcript abundance
#' @param significance_threshold A real between 0 and 1
#'
#' @return A tibble with edgeR results
#'
test_stratification_cellularity_ <- function(.data,
																						 .formula,
																						 .sample = NULL,
																						 .transcript = NULL,
																						 .abundance = NULL,
																						 method = "cibersort",
																						 reference = NULL,
																						 ...
) {

	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)


	deconvoluted =
		.data %>%

		# Deconvolution
		deconvolve_cellularity(
			!!.sample, !!.transcript, !!.abundance,
			method=method,
			prefix = sprintf("%s:", method),
			reference = reference,
			action="get",
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
			univariable_differential_tissue_stratification(deconvoluted,
																										 method,
																										 .my_formula) %>%

				# Attach attributes
				reattach_internals(.data) %>%

				# Add methods used
				memorise_methods_used(c("survival", "boot", "survminer"), object_containing_methods = .data)
		} %>%

		# Eliminate prefix
		mutate(.cell_type = str_remove(.cell_type, sprintf("%s:", method)))


}



#' Get gene enrichment analyses using EGSEA
#'
#' @keywords internal
#' @noRd
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr set_colnames
#' @importFrom purrr map2_dfr
#' @importFrom stats model.matrix
#' @importFrom utils install.packages
#'
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .formula A formula with no response variable, representing the desired linear model
#' @param .sample The name of the sample column
#' @param .entrez The ENTREZ code of the transcripts/genes
#' @param .abundance The name of the transcript/gene abundance column
#' @param .contrasts A character vector. See edgeR makeContrasts specification for the parameter `contrasts`. If contrasts are not present the first covariate is the one the model is tested against (e.g., ~ factor_of_interest)
#' @param methods A character vector. One or 3 or more methods to use in the testing (currently EGSEA errors if 2 are used). Type EGSEA::egsea.base() to see the supported GSE methods.
#' @param gene_sets A character vector or a list. It can take one or more of the following built-in collections as a character vector: c("h", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "kegg_disease", "kegg_metabolism", "kegg_signaling"), to be used with EGSEA buildIdx. c1 is human specific. Alternatively, a list of user-supplied gene sets can be provided, to be used with EGSEA buildCustomIdx. In that case, each gene set is a character vector of Entrez IDs and the names of the list are the gene set names.
#' @param species A character. It can be human, mouse or rat.
#' @param cores An integer. The number of cores available
#'
#' @return A tibble with edgeR results
#'
test_gene_enrichment_bulk_EGSEA <- function(.data,
																							 .formula,
																							 .sample = NULL,
																							 .entrez,
																							 .abundance = NULL,
																							 .contrasts = NULL,
																						     methods,
																							 gene_sets,
																							 species,
																							 cores = 10) {
	# Comply with CRAN NOTES
	. = NULL

	# Get column names
	.sample = enquo(.sample)
	.abundance = enquo(.abundance)

	.entrez = enquo(.entrez)

	# distinct_at is not released yet for dplyr, thus we have to use this trick
	df_for_edgeR <- .data %>%

		# Prepare the data frame
		select(!!.entrez, !!.sample, !!.abundance,
					 one_of(parse_formula(.formula))) %>%
		distinct() %>%

		# Add entrez from symbol
		filter(!!.entrez %>% is.na %>% not())

	# Check if at least two samples for each group
	if (df_for_edgeR %>%
			select(!!.sample, one_of(parse_formula(.formula))) %>%
			distinct %>%
			count(!!as.symbol(parse_formula(.formula))) %>%
			distinct(n) %>%
			pull(n) %>%
			min %>%
			st(2))
		stop("tidybulk says: You need at least two replicates for each condition for EGSEA to work")


	# Create design matrix
	design =
		model.matrix(
			object = .formula,
			data = df_for_edgeR %>% select(!!.sample, one_of(parse_formula(.formula))) %>% distinct %>% arrange(!!.sample)
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
	if (find.package("EGSEA", quiet = TRUE) %>% length %>% equals(0)) {
		stop("
				 EGSEA not installed. Please install it. EGSEA requires manual installation to not overwhelm the user in case it is not needed.
				 BiocManager::install(\"EGSEA\", ask = FALSE)
				 ")
	}
	if (!"EGSEA" %in% (.packages())) {
		stop("EGSEA package not loaded. Please run library(\"EGSEA\"). With this setup, EGSEA require manual loading, for technical reasons.")
	}

	dge =
		df_for_edgeR %>%

		# Make sure transcript names are adjacent
		arrange(!!.entrez) %>%

		select(!!.sample, !!.entrez, !!.abundance) %>%
		spread(!!.sample,!!.abundance) %>%
		as_matrix(rownames = !!.entrez) %>%
		edgeR::DGEList(counts = .)

	# Add gene ids for Interpret Results tables in report
    dge$genes = rownames(dge$counts)


	if (is.list(gene_sets)) {

	    idx =  buildCustomIdx(geneIDs = rownames(dge), species = species, gsets=gene_sets)
	    nonkegg_genesets = idx
	    kegg_genesets = NULL

	} else {

    	# Check gene sets to include
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


    	idx =  buildIdx(entrezIDs = rownames(dge), species = species,  msigdb.gsets = msigdb.gsets,
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
    	message("tidybulk says: due to a bug in the call to KEGG database (http://supportupgrade.bioconductor.org/p/122172/#122218), the analysis for this database is run without report production.")

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
	    out %>% memorise_methods_used(c("egsea", collections_bib, methods), object_containing_methods = .data)
	}
}

#' Get K-mean clusters to a tibble
#'
#' @keywords internal
#' @noRd
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
		if ("centers" %in% names(dots_args) %>% not())
			stop("tidybulk says: for kmeans you need to provide the \"centers\" integer argument")

		# Get column names
		.element = enquo(.element)
		.feature = enquo(.feature)
		.abundance = enquo(.abundance)

		.data %>%

			# Prepare data frame
			distinct(!!.feature,!!.element,!!.abundance) %>%

			# Check if log transform is needed
			ifelse_pipe(log_transform,
									~ .x %>% dplyr::mutate(!!.abundance := !!.abundance %>%  `+`(1) %>%  log())) %>%

			# Prepare data frame for return
			spread(!!.feature,!!.abundance) %>%
			as_matrix(rownames = !!.element) %>%

			# Wrap the do.call because of the centrers check
			{
				do.call(kmeans, list(x = (.), iter.max = 1000) %>% c(dots_args))
			}	 %$%
			cluster %>%
			as.list() %>%
			as_tibble() %>%
			gather(!!.element, `cluster kmeans`) %>%
			mutate(`cluster kmeans` = `cluster kmeans` %>% as.factor()) %>%

			# Attach attributes
			reattach_internals(.data) %>%
			memorise_methods_used("stats")
	}

#' Get SNN shared nearest neighbour clusters to a tibble
#'
#' @keywords internal
#' @noRd
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
		.abundance = enquo(.abundance)

		# Check if package is installed, otherwise install
		if (find.package("cluster", quiet = TRUE) %>% length %>% equals(0)) {
			message("tidybulk says: Installing cluster")
			install.packages("cluster", repos = "https://cloud.r-project.org")
		}
		if (find.package("Seurat", quiet = TRUE) %>% length %>% equals(0)) {
			message("tidybulk says: Installing Seurat")
			install.packages("Seurat", repos = "https://cloud.r-project.org")
		}
		if (find.package("KernSmooth", quiet = TRUE) %>% length %>% equals(0)) {
			message("tidybulk says: Installing KernSmooth")
			install.packages("KernSmooth", repos = "https://cloud.r-project.org")
		}

		my_df =
			.data %>%

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
			Seurat::RunPCA(npcs = 10) %>%
			Seurat::FindNeighbors() %>%
			Seurat::FindClusters(method = "igraph", ...) %>%
			`[[` ("seurat_clusters") %>%
			as_tibble(rownames = quo_name(.element)) %>%
			rename(`cluster SNN` = seurat_clusters) %>%
			dplyr::mutate(!!.element := gsub("\\.", "-",!!.element)) %>%

			# Attach attributes
			reattach_internals(.data) %>%
			memorise_methods_used("seurat")
	}

#' Get dimensionality information to a tibble using MDS
#'
#' @keywords internal
#' @noRd
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom purrr map_dfr
#' @importFrom rlang :=
#' @importFrom stats setNames
#' @importFrom rlang inform
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
		.abundance = enquo(.abundance)

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

					distinct(!!.feature,!!.element,!!.abundance) %>%

					# Check if log transform is needed
					ifelse_pipe(log_transform,
											~ .x %>% dplyr::mutate(!!.abundance := !!.abundance %>% log1p())) %>%

					# Stop any column is not if not numeric or integer
					ifelse_pipe(
						(.) %>% select(!!.abundance) %>% summarise_all(class) %>% `%in%`(c("numeric", "integer")) %>% `!`() %>% any(),
						~ stop(".abundance must be numerical or integer")
					) %>%
					spread(!!.element,!!.abundance) %>%
					as_matrix(rownames = !!.feature, do_check = FALSE) %>%
					limma::plotMDS(dim.plot = .x, plot = FALSE, top = top)
			)

		map2_dfr(
			mds_object, components_list,
			~ 	{
				# Change of function from Bioconductor 3_13 of plotMDS
				my_rownames = .x %>% when(
					"distance.matrix.squared" %in% names(.x) ~ .x$distance.matrix.squared,
					~ .x$distance.matrix
				) %>%
					rownames()

				tibble(my_rownames, .x$x, .x$y) %>%
					rename(
						!!.element := my_rownames,
						!!as.symbol(.y[1]) := `.x$x`,
						!!as.symbol(.y[2]) := `.x$y`
					) %>%
					gather(Component, `Component value`,-!!.element)

				}

		)  %>%
			distinct() %>%
			spread(Component, `Component value`) %>%
			setNames(c((.) %>% select(1) %>% colnames(),
								 paste0("Dim", (.) %>% select(-1) %>% colnames())
			)) %>%


			# Attach attributes
			reattach_internals(.data) %>%
			memorise_methods_used("limma") %>%

			# Add raw object
			attach_to_internals(mds_object, "MDS") %>%
			# Communicate the attribute added
			{

			  rlang::inform("tidybulk says: to access the raw results do `attr(..., \"internals\")$MDS`", .frequency_id = "Access MDS results",  .frequency = "once")

				(.)
			}
	}

#' Get principal component information to a tibble using PCA
#'
#' @keywords internal
#' @noRd
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang :=
#' @importFrom stats prcomp
#' @importFrom utils capture.output
#' @importFrom magrittr divide_by
#' @importFrom rlang inform
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
					 scale = FALSE,
					 ...) {
		# Comply with CRAN NOTES
		. = NULL

		# Get column names
		.element = enquo(.element)
		.feature = enquo(.feature)
		.abundance = enquo(.abundance)

		# Get components from dims
		components = 1:.dims

		prcomp_obj =
			.data %>%

			# Filter most variable genes
			keep_variable_transcripts(!!.element,!!.feature,!!.abundance, top) %>%

			# Prepare data frame
			distinct(!!.feature,!!.element,!!.abundance) %>%

			# Check if log transform is needed
			ifelse_pipe(log_transform,
									~ .x %>% dplyr::mutate(!!.abundance := !!.abundance %>% log1p())) %>%

			# Stop any column is not if not numeric or integer
			ifelse_pipe(
				(.) %>% select(!!.abundance) %>% summarise_all(class) %>% `%in%`(c("numeric", "integer", "bouble")) %>% not() %>% any(),
				~ stop("tidybulk says: .abundance must be numerical or integer")
			) %>%

			spread(!!.element,!!.abundance) %>%

			drop_na %>% # Is this necessary?

			# check that there are non-NA genes for enough samples
			ifelse2_pipe(# First condition
				(.) %>% nrow == 0,

				# Second condition
				(.) %>% nrow < 100,

				# First function
				~ stop(
					"tidybulk says: In calculating PCA there is no gene that have non NA values is all samples"
				),

				# Second function
				~ {
					warning(
						"
						tidybulk says: In PCA correlation there is < 100 genes that have non NA values is all samples.
The correlation calculation would not be reliable,
we suggest to partition the dataset for sample clusters.
						"
					)
					.x
				}) %>%

			# Transform to matrix
			as_matrix(rownames = !!.feature, do_check = FALSE) %>%
			t() %>%

			# Calculate principal components
			prcomp(scale = scale, ...)

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
			as_tibble(rownames = quo_name(.element)) %>%
			select(!!.element, sprintf("PC%s", components)) %>%

			# Attach attributes
			reattach_internals(.data) %>%
			memorise_methods_used("stats") %>%

			# Add raw object
			attach_to_internals(prcomp_obj, "PCA") %>%
			# Communicate the attribute added
			{
			  rlang::inform("tidybulk says: to access the raw results do `attr(..., \"internals\")$PCA`", .frequency_id = "Access PCA results",  .frequency = "once")

				(.)
			}

	}

#' Get tSNE
#'
#' @keywords internal
#' @noRd
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
get_reduced_dimensions_TSNE_bulk <-
	function(.data,
					 .element = NULL,
					 .feature = NULL,

					 .abundance = NULL,
					 .dims = 2,
					 top = 500,
					 of_samples = TRUE,
					 log_transform = TRUE,
					 ...) {
		# Comply with CRAN NOTES
		. = NULL

		# Get column names
		.element = enquo(.element)
		.feature = enquo(.feature)
		.abundance = enquo(.abundance)

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
			message("tidybulk says: Installing Rtsne")
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


		df_tsne =
			.data %>%

			# Filter NA symbol
			filter(!!.feature %>% is.na %>% not()) %>%

			# Prepare data frame
			distinct(!!.feature,!!.element,!!.abundance) %>%

			# Check if data rectangular
			ifelse_pipe(
				(.) %>% check_if_data_rectangular(!!.element,!!.feature,!!.abundance),
				~ .x %>% eliminate_sparse_transcripts(!!.feature)
			) %>%

			# Check if log transform is needed
			ifelse_pipe(log_transform,
									~ .x %>% dplyr::mutate(!!.abundance := !!.abundance %>% log1p)) %>%

			# Filter most variable genes
			keep_variable_transcripts(!!.element,!!.feature,!!.abundance, top) %>%

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
			reattach_internals(.data) %>%
			memorise_methods_used("rtsne")

	}

#' Get UMAP
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
#' @param calculate_for_pca_dimensions An integer of length one. The number of PCA dimensions to based the UMAP calculatio on. If NULL all variable features are considered
#' @param ... Further parameters passed to the function uwot
#'
#' @return A tibble with additional columns
#'
get_reduced_dimensions_UMAP_bulk <-
  function(.data,
           .element = NULL,
           .feature = NULL,

           .abundance = NULL,
           .dims = 2,
           top = 500,
           of_samples = TRUE,
           log_transform = TRUE,
           scale = TRUE,
           calculate_for_pca_dimensions = 20,
           ...) {

    if(!is.null(calculate_for_pca_dimensions) & (
      !is(calculate_for_pca_dimensions, "numeric") |
      length(calculate_for_pca_dimensions) > 1
    ))
      stop("tidybulk says: the argument calculate_for_pca_dimensions should be NULL or an integer of size 1")

    # Comply with CRAN NOTES
    . = NULL

    # Get column names
    .element = enquo(.element)
    .feature = enquo(.feature)
    .abundance = enquo(.abundance)

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

    df_source =
      .data %>%

      # Filter NA symbol
      filter(!!.feature %>% is.na %>% not()) %>%

      # Prepare data frame
      distinct(!!.feature,!!.element,!!.abundance) %>%

      # Check if data rectangular
      when(
        check_if_data_rectangular(., !!.element,!!.feature,!!.abundance)  ~
          eliminate_sparse_transcripts(., !!.feature),
        ~ (.)
      ) %>%

      # Check if log transform is needed
      when(log_transform    ~ dplyr::mutate(., !!.abundance := !!.abundance %>% log1p), ~ (.)) %>%

      # Filter most variable genes
      keep_variable_transcripts(!!.element,!!.feature,!!.abundance, top)

    # Calculate based on PCA
    if(!is.null(calculate_for_pca_dimensions))
      df_UMAP =
      df_source %>%
      reduce_dimensions(
        !!.element,!!.feature,!!.abundance,
        method="PCA",
        .dims = calculate_for_pca_dimensions,
        action="get",
        scale = scale
      ) %>%
      suppressMessages() %>%
      as_matrix(rownames = quo_name(.element))

    # Calculate based on all features
    else
      df_UMAP =
      df_source %>%
      spread(!!.feature,!!.abundance) %>%
      as_matrix(rownames = quo_name(.element))

    do.call(uwot::tumap, c(list(df_UMAP), arguments)) %>%
      as_tibble(.name_repair = "minimal") %>%
      setNames(c("UMAP1", "UMAP2")) %>%

      # add element name
      dplyr::mutate(!!.element := df_UMAP %>% rownames) %>%
      select(!!.element, everything()) %>%

      # Attach attributes
      reattach_internals(.data) %>%
      memorise_methods_used("uwot")

  }

#' Get rotated dimensions of two principal components or MDS dimension of choice, of an angle
#'
#' @keywords internal
#' @noRd
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
		dimension_1_column = enquo(dimension_1_column)
		dimension_2_column = enquo(dimension_2_column)
		dimension_1_column_rotated = enquo(dimension_1_column_rotated)
		dimension_2_column_rotated = enquo(dimension_2_column_rotated)

		if (.data %>%
				distinct(!!.element, !!dimension_1_column, !!dimension_2_column) %>%
				count(!!.element, !!dimension_1_column, !!dimension_2_column) %>%
				pull(n) %>%
				max %>%
				gt(1))
			stop(sprintf(
				"tidybulk says: %s must be unique for each row for the calculation of rotation",
				quo_name(.element)
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
		if (rotation_degrees %>% between(-360, 360) %>% not())
			stop("tidybulk says: rotation_degrees must be between -360 and 360")

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

#' Aggregates multiple counts from the same samples (e.g., from isoforms)
#' This function aggregates counts over samples, concatenates other character columns, and averages other numeric columns
#'
#' @keywords internal
#' @noRd
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
			warning("tidybulk says: for aggregation, factors and logical columns were converted to character")
			message("Converted to characters")
			message(lapply(.data, class) %>% unlist %>% `[` (. %in% c("logical", "factor") %>% which))
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
					# .data %>% get_tt_columns() %$% .abundance_scaled %>% is.null %>% not() &&
					quo_name(.data %>% get_tt_columns() %$% .abundance_scaled) %in% (.data %>% colnames)
			),
			~ .x %>% select(-!!(
				.data %>% get_tt_columns() %$% .abundance_scaled
			)))	%>%
			colnames() %>%
			c("n_aggr")

		# aggregates read .data over samples, concatenates other character columns, and averages other numeric columns
		.data %>%

			# transform logicals and factors
			mutate_if(is.factor, as.character) %>%
			mutate_if(is.logical, as.character) %>%

			# Add the number of duplicates for each gene
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

						# If scaled abundance exists aggregate that as well
						ifelse_pipe((
							".abundance_scaled" %in% (.data %>% get_tt_columns() %>% names) &&
								# .data %>% get_tt_columns() %$% .abundance_scaled %>% is.null %>% not() &&
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
			rename(merged_transcripts = n_aggr) %>%

			# Attach attributes
			reattach_internals(.data)

	}

#' Aggregates multiple counts from the same samples (e.g., from isoforms)
#' This function aggregates counts over samples, concatenates other character columns, and averages other numeric columns
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom dplyr summarise_all
#' @importFrom dplyr bind_rows
#' @importFrom magrittr %$%
#' @importFrom rlang :=
#' @importFrom stringi stri_c
#' @importFrom tidyr replace_na
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
aggregate_duplicated_transcripts_DT =
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

    #
    #         data.table::setDT(.data)
    #
    #         input_df[,.(c = sum(c)), by=c("a", "b")] %>% #Sum the raw_counts of duplicated rows
    #       tidybulk::scale_abundance(.sample = sample, .abundance = abundance, .transcript = transcript)

    if(.data %>% filter(is.na(!!.transcript)) %>% nrow() %>% gt(0)){
      warning(sprintf("tidybulk says: some of your %s are NAs. Those will be eliminated to correctly aggregate the duplicates", quo_name(.transcript)))
      .data = .data %>% filter(!is.na(!!.transcript))
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
          # .data %>% get_tt_columns() %$% .abundance_scaled %>% is.null %>% not() &&
          quo_name(.data %>% get_tt_columns() %$% .abundance_scaled) %in% (.data %>% colnames)
      ),
      ~ .x %>% select(-!!(
        .data %>% get_tt_columns() %$% .abundance_scaled
      )))	%>%
      colnames()

    aggregate_count_columns =
      quo_name(.abundance) %>%
      when(
        ".abundance_scaled" %in% (.data %>% get_tt_columns() %>% names) &&
          quo_name(.data %>% get_tt_columns() %$% .abundance_scaled) %in% (.data %>% colnames)  ~
          (.) %>% c(.data %>% get_tt_columns() %$% .abundance_scaled),
        ~ (.)
      )

    pasted_strings___ = stringi::stri_c(pull(.data,quo_name(.transcript)), pull(.data,quo_name(.sample)), sep = "_")
    #.data = .data %>% mutate(pasted_strings___ = pasted_strings___)
    duplicates = pasted_strings___%in%pasted_strings___[which(duplicated(pasted_strings___))]

    dup = .data %>% filter(duplicates)
    dup_pasted_strings___ = pasted_strings___[duplicates]

    dup_counts =
      dup%>%
      group_by(!!.sample,!!.transcript) %>%
      dplyr::summarise(
        across(aggregate_count_columns, ~ .x %>% aggregation_function()),
        merged_transcripts = n()
      ) %>%
      ungroup()

    dup =
      dup_counts %>%

      # Bind cols with dropped duplicated
      left_join(
        dup[!duplicated(dup_pasted_strings___),] %>%
          select(-aggregate_count_columns)
        )

    .data %>%
      filter(!duplicates) %>%
      bind_rows(dup) %>%
      replace_na(list(merged_transcripts = 1)) %>%
      select(.data %>% colnames() %>% c("merged_transcripts"))

  }

#' Drop redundant elements (e.g., samples) for which feature (e.g., genes) aboundances are correlated
#'
#' @keywords internal
#' @noRd
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
#' @return A tibble with redundant elements removed
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

	# Check if .data has more than one element
	if(.data %>% distinct(!!.element) %>% nrow() <= 1 )
		stop("tidybulk says: You must have more than one element (trancripts if of_samples == FALSE) to perform remove_redundancy")

	# Check if package is installed, otherwise install
	if (find.package("widyr", quiet = TRUE) %>% length %>% equals(0)) {
		message("tidybulk says: Installing widyr needed for correlation analyses")
		install.packages("widyr", repos = "https://cloud.r-project.org")
	}

	# Get the redundant data frame
	.data.correlated =
		.data %>%

		# Prepare the data frame
		select(!!.feature,!!.element,!!.abundance) %>%

		# Filter variable genes
		keep_variable_transcripts(!!.element,!!.feature,!!.abundance, top = top) %>%

		# Check if log transform is needed
		ifelse_pipe(log_transform,
								~ .x %>% dplyr::mutate(!!.abundance := !!.abundance %>% `+`(1) %>%  log())) %>%
		distinct() %>%

# NO NEED OF RECTANGULAR
# 		spread(!!.element,!!.abundance) %>%
# 		drop_na() %>%
#
# 		# check that there are non-NA genes for enough samples
# 		ifelse2_pipe(# First condition
# 			(.) %>% nrow == 0,
#
# 			# Second condition
# 			(.) %>% nrow < 100,
#
# 			# First function
# 			~ stop(
# 				"tidybulk says: In calculating correlation there is no gene that have non NA values is all samples"
# 			),
#
# 			# Second function
# 			~ {
# 				message(
# 					"tidybulk says: In calculating correlation there is < 100 genes (that have non NA values) is all samples.
# The correlation calculation might not be reliable"
# 				)
# 				.x
# 			}) %>%
#
# 		# Prepare the data frame
# 		gather(!!.element,!!.abundance,-!!.feature) %>%

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

	# Return non redundant data frame
	.data %>% anti_join(.data.correlated, by = quo_name(.element)) %>%

		# Attach attributes
		reattach_internals(.data) %>%
		memorise_methods_used("widyr")
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

		# Check if .data has more than one element
		if(.data %>% distinct(!!.element) %>% nrow() <= 1 )
			stop("tidybulk says: You must have more than one element (trancripts if of_samples == FALSE) to perform remove_redundancy")

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
		.data %>% anti_join(.data.redundant, by = quo_name(.element)) %>%

			# Attach attributes
			reattach_internals(.data)
	}

##' after wget, this function merges hg37 and hg38 mapping data bases - Do not execute!
##'
##' @return A tibble with ensembl-transcript mapping
##'
# get_ensembl_symbol_mapping <- function() {
# 	# wget -O mapping_38.txt 'http://www.ensembl.org/biomart/martservice?query=  <Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6">  <Dataset name="hsapiens_gene_ensembl" interface="default"> <Attribute name="ensembl_transcript_id"/> <Attribute name="ensembl_gene_id"/><Attribute name="transcript"/> </Dataset> </Query>'
# 	# wget -O mapping_37.txt 'http://grch37.ensembl.org/biomart/martservice?query=<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6"><Dataset name="hsapiens_gene_ensembl" interface="default"><Attribute name="ensembl_transcript_id"/><Attribute name="ensembl_gene_id"/><Attribute name="transcript"/></Dataset></Query>'
# 	all =
# 		read_table2("~/third_party_sofware/ensembl_mapping/mapping_37.txt",
# 								col_names = FALSE) %>%
# 		setNames(c("ensembl_transcript_id", "ensembl_gene_id", "transcript")) %>%
# 		mutate(ref_genome = "hg37") %>%
# 		bind_rows(
# 			read_table2(
# 				"~/third_party_sofware/ensembl_mapping/mapping_38.txt",
# 				col_names = FALSE
# 			) %>%
# 				setNames(
# 					c("ensembl_transcript_id", "ensembl_gene_id", "transcript")
# 				) %>%
# 				mutate(ref_genome = "hg38")
# 		) %>%
# 		drop_na()
#
#
# 	bind_rows(
# 		# Gene annotation
# 		all %>%
# 			select(-ensembl_transcript_id) %>%
# 			group_by(ensembl_gene_id) %>%
# 			arrange(ref_genome %>% desc()) %>%
# 			slice(1) %>%
# 			ungroup() %>%
# 			rename(ensembl_id = ensembl_gene_id),
#
# 		# Transcript annotation
# 		all %>%
# 			select(-ensembl_gene_id) %>%
# 			group_by(ensembl_transcript_id) %>%
# 			arrange(ref_genome %>% desc()) %>%
# 			slice(1) %>%
# 			ungroup() %>%
# 			rename(ensembl_id = ensembl_transcript_id)
# 	) %>%
#
# 		# Write to file and return
# 		{
# 			(.) %>% write_csv("~/third_party_sofware/ensembl_mapping/ensembl_symbol_mapping.csv")
# 			(.)
# 		}
# }

##' after wget, this function merges hg37 and hg38 mapping data bases - Do not execute!
##'
##' @return A tibble with ensembl-transcript mapping
##'
# get_ensembl_symbol_mapping_mouse <- function() {
#
#
# 	left_join(
# 		AnnotationDbi::toTable(org.Mm.eg.db::org.Mm.egENSEMBL),
# 		AnnotationDbi::toTable(org.Mm.eg.db::org.Mm.egSYMBOL)
# 	) %>%
# 		as_tibble %>%
# 		select(-gene_id) %>%
# 		rename(ensembl_id = ensembl_id, transcript = symbol) %>%
# 		drop_na() %>%
#
# 		mutate(ref_genome = "mm10") %>%
#
# 		# Write to file and return
# 		{
# 			(.) %>% write_csv("~/third_party_sofware/ensembl_mapping/ensembl_symbol_mapping_mouse.csv")
# 			(.)
# 		}
# }

#' get_symbol_from_ensembl
#'
#' @keywords internal
#' @noRd
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

		# # Creating file
		# ensembl_symbol_mapping =
		# 	bind_rows(
		# 		read_csv("~/third_party_sofware/ensembl_mapping/ensembl_symbol_mapping.csv"),
		# 		read_csv("~/third_party_sofware/ensembl_mapping/ensembl_symbol_mapping_mouse.csv")
		# 	)
		#
		# save(ensembl_symbol_mapping, file="data/ensembl_symbol_mapping.rda", compress = "xz")


		.ensembl = enquo(.ensembl)

		.data %>%
			select(!!.ensembl) %>%
			distinct() %>%

			# Add name information
			dplyr::left_join(
				tidybulk::ensembl_symbol_mapping %>%
					distinct(ensembl_id, transcript, ref_genome) %>%
					dplyr::rename(!!.ensembl := ensembl_id) %>%
					rename(transcript = transcript),
				by = quo_name(.ensembl)
			)

	}

#' Perform linear equation system analysis through llsr
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom stats lsfit
#'
#' @param mix A data frame
#' @param reference A data frame
#'
#' @return A data frame
#'
#'
run_llsr = function(mix, reference = X_cibersort,  intercept= TRUE) {
	# Get common markers
	markers = intersect(rownames(mix), rownames(reference))

	X <- (reference[markers, , drop = FALSE])
	Y <- (mix[markers, , drop = FALSE])

	X = as.matrix(X)
	Y = as.matrix(Y)

	X <- (X - mean(X)) / sd(X)
	Y <- apply(Y, 2, function(mc) (mc - mean(mc)) / sd(mc)  )
	# Y <- (Y - mean(y)) / sd(Y)

	if(intercept)
	  results <- t(data.frame(lsfit(X, Y)$coefficients)[-1, , drop = FALSE])
	else
	  results <- t(data.frame(lsfit(X, Y, intercept=FALSE)$coefficients))
	results[results < 0] <- 0
	results <- results / apply(results, 1, sum)
	rownames(results) = colnames(Y)

	results
}

#' Perform linear equation system analysis through llsr
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom stats lsfit
#'
#' @param mix A data frame
#' @param reference A data frame
#'
#' @return A data frame
#'
#'
run_epic = function(mix, reference = NULL) {

	# Check if package is installed, otherwise install
	if (find.package("devtools", quiet = TRUE) %>% length %>% equals(0)) {
		message("tidybulk says: Installing class needed for EPIC")
		install.packages("devtools", repos = "https://cloud.r-project.org", dependencies = c("Depends", "Imports"))
	}

	# Check if package is installed, otherwise install
	if (find.package("EPIC", quiet = TRUE) %>% length %>% equals(0)) {
		message("tidybulk says: Installing class needed for EPIC")
		devtools::install_github("GfellerLab/EPIC")
	}

	if("EPIC" %in% .packages() %>% not) stop("tidybulk says: Please install and then load the package EPIC manually (i.e. library(EPIC)). This is because EPIC is not in Bioconductor or CRAN so it is not possible to seamlessly make EPIC part of the dependencies.")

	# Get common markers
	if( reference %>% class %>% equals("data.frame")){
		markers = intersect(rownames(mix), rownames(reference))

		X <- (reference[markers, , drop = FALSE])
		Y <- (mix[markers, , drop = FALSE])

		if(!is.null(reference))
			reference = list(
				refProfiles = X,
				sigGenes = rownames(X)
			)
	} else { Y <- mix }


	results <- EPIC(Y, reference = reference)$cellFractions %>% data.frame()
	#results[results < 0] <- 0
	#results <- results / apply(results, 1, sum)
	rownames(results) = colnames(Y)

	results
}


#' Get cell type proportions from cibersort
#'
#' @keywords internal
#' @noRd
#'
#' @import parallel
#' @import preprocessCore
#' @importFrom stats setNames
#' @importFrom rlang dots_list
#' @importFrom magrittr equals
#' @importFrom utils install.packages
#'
#' @param .data A tibble
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param reference A data frame. The transcript/cell_type data frame of integer transcript abundance
#' @param method A character string. The method to be used. At the moment Cibersort (default) and llsr (linear least squares regression) are available.
#' @param prefix A character string. The prefix you would like to add to the result columns. It is useful if you want to reshape data.
#' @param ... Further parameters passed to the function Cibersort
#'
#' @return A tibble including additional columns
#'
#'
get_cell_type_proportions = function(.data,
																		 .sample = NULL,
																		 .transcript = NULL,
																		 .abundance = NULL,
																		 reference = NULL,
																		 method = "cibersort",
																		 prefix = "",
																		 ...) {
	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	# Get the dots arguments
	dots_args = rlang::dots_list(...)

	.data %>%

		# Prepare data frame
		distinct(!!.sample,!!.transcript,!!.abundance) %>%
		spread(!!.sample,!!.abundance) %>%

		# Eliminate NA transcripts
		filter(!!.transcript %>% is.na %>% not()) %>%

		# Convert
		data.frame(row.names = 1, check.names = FALSE) %>%

		# Run Cibersort or llsr through custom function, depending on method choice
		when(

			# Execute do.call because I have to deal with ...
			method %>% tolower %>% equals("cibersort") 	~ {

				# Check if package is installed, otherwise install
				if (find.package("class", quiet = TRUE) %>% length %>% equals(0)) {
					message("tidybulk says: Installing class needed for Cibersort")
					install.packages("class", repos = "https://cloud.r-project.org", dependencies = c("Depends", "Imports"))
				}

				# Check if package is installed, otherwise install
				if (find.package("e1071", quiet = TRUE) %>% length %>% equals(0)) {
					message("tidybulk says: Installing e1071 needed for Cibersort")
					install.packages("e1071", repos = "https://cloud.r-project.org", dependencies = c("Depends", "Imports"))
				}

				# Check if package is installed, otherwise install
				if (find.package("preprocessCore", quiet = TRUE) %>% length %>% equals(0)) {
					message("tidybulk says: Installing preprocessCore needed for Cibersort")
					if (!requireNamespace("BiocManager", quietly = TRUE))
						install.packages("BiocManager", repos = "https://cloud.r-project.org")
					BiocManager::install("preprocessCore", ask = FALSE)

				}

				# Choose reference
				reference = reference %>% when(is.null(.) ~ X_cibersort, ~ .)

				# Validate reference
				validate_signature(.data, reference, !!.transcript)

				do.call(my_CIBERSORT, list(Y = ., X = reference, QN=FALSE) %>% c(dots_args)) %$%
				proportions %>%
				as_tibble(rownames = quo_name(.sample)) %>%
				select(-`P-value`,-Correlation,-RMSE)
			},

			# Don't need to execute do.call
			method %>% tolower %>% equals("llsr") ~ {

				# Choose reference
				reference = reference %>% when(is.null(.) ~ X_cibersort, ~ .)

				# Validate reference
				validate_signature(.data, reference, !!.transcript)

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

				# # Check if package is installed, otherwise install
				if (find.package("immunedeconv", quiet = TRUE) %>% length %>% equals(0)) {
					message("tidybulk says: Installing immunedeconv")
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
			quo_name(.sample),
			(.) %>% select(-1) %>% colnames() %>% sprintf("%s%s", prefix, .)

		)) %>%

		# Attach attributes
		reattach_internals(.data) %>%
		memorise_methods_used(tolower(method))


}

#' Get adjusted count for some batch effect
#'
#' @keywords internal
#' @noRd
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
	.abundance = enquo(.abundance)

	# Check that .formula includes at least two covariates
	if (parse_formula(.formula) %>% length %>% st(2))
		stop(
			"The .formula must contain two covariates, the first being the factor of interest, the second being the factor of unwanted variation"
		)

	# Check that .formula includes no more than two covariates at the moment
	if (parse_formula(.formula) %>% length %>% gt(3))
		warning("tidybulk says: Only the second covariate in the .formula is adjusted for, at the moment")

	# New column name
	value_adjusted = as.symbol(sprintf("%s%s",  quo_name(.abundance), adjusted_string))

	df_for_combat <-
		.data %>%

		select(!!.transcript,
					 !!.sample,
					 !!.abundance,
					 one_of(parse_formula(.formula))) %>%
		distinct() %>%

		# Check if log transform is needed
		ifelse_pipe(log_transform,
								~ .x %>% dplyr::mutate(!!.abundance := !!.abundance %>% log1p()))


	# Create design matrix
	design =
		model.matrix(
			object = as.formula("~" %>% paste0(parse_formula(.formula)[1])),
			# get first argument of the .formula
			data = df_for_combat %>% select(!!.sample, one_of(parse_formula(.formula))) %>% distinct %>% arrange(!!.sample)
		)

	# Maybe not needed and causing trouble if more columns that in the formula
	  # %>%
		#set_colnames(c("(Intercept)", parse_formula(.formula)[1]))

	# Check if package is installed, otherwise install
	if (find.package("sva", quiet = TRUE) %>% length %>% equals(0)) {
		message("tidybulk says: Installing sva - Combat needed for adjustment for unwanted variation")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("sva", ask = FALSE)
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
			(.) %>% select(!!.abundance) %>% summarise_all(class) %>% `%in%`(c("numeric", "integer")) %>% not() %>% any(),
			~ stop(".abundance must be numerical or integer")
		) %>%

		spread(!!.sample,!!.abundance) %>%
		as_matrix(rownames = !!.transcript,
							do_check = FALSE)
	mat %>%

		# Add little noise to avoid all 0s for a covariate that would error combat code (not statistics that would be fine)
		`+` (rnorm(length(mat), 0, 0.000001)) %>%

		# Run combat
		sva::ComBat(batch = my_batch,
								mod = design,
								prior.plots = FALSE,
								...) %>%

		as_tibble(rownames = quo_name(.transcript)) %>%
		gather(!!.sample,!!.abundance,-!!.transcript) %>%

		# Reverse-Log transform if transformed in the first place
		ifelse_pipe(
			log_transform,
			~ .x %>%
				dplyr::mutate(!!.abundance := !!.abundance %>% exp %>% `-`(1)) %>%
				dplyr::mutate(!!.abundance := ifelse(!!.abundance < 0, 0,!!.abundance)) %>%
				dplyr::mutate(!!.abundance := !!.abundance %>% as.integer)
		) %>%

		# Reset column names
		dplyr::rename(!!value_adjusted := !!.abundance)  %>%

		# Attach attributes
		reattach_internals(.data) %>%
		memorise_methods_used("sva")
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
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#'
#' @return A tibble filtered genes
#'
keep_variable_transcripts = function(.data,
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

	message(sprintf("Getting the %s most variable genes", top))

	x =
		.data %>%
		distinct(!!.sample,!!.transcript,!!.abundance) %>%

		# Check if logtansform is needed
		ifelse_pipe(log_transform,
								~ .x %>% dplyr::mutate(!!.abundance := !!.abundance %>% `+`(1) %>%  log())) %>%

		spread(!!.sample,!!.abundance) %>%
		as_matrix(rownames = quo_name(.transcript))

	s <- rowMeans((x - rowMeans(x, na.rm = TRUE)) ^ 2, na.rm = TRUE)
	o <- order(s, decreasing = TRUE)
	x <- x[o[1L:top], , drop = FALSE]
	variable_trancripts = rownames(x)

	.data %>%
		filter(!!.transcript %in% variable_trancripts) %>%

		# Add methods used. The correlation code comes from there
		memorise_methods_used(c("edger"))
}

#' tidybulk_to_SummarizedExperiment
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom utils data
#' @importFrom tidyr pivot_longer
#'
#' @param .data A tibble
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#'
#' @return A SummarizedExperiment
#'
tidybulk_to_SummarizedExperiment = function(.data,
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
	if (find.package("SummarizedExperiment", quiet = TRUE) %>% length %>% equals(0)) {
		message("tidybulk says: Installing SummarizedExperiment")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("SummarizedExperiment", ask = FALSE)
	}
	if (find.package("S4Vectors", quiet = TRUE) %>% length %>% equals(0)) {
		message("tidybulk says: Installing S4Vectors")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("S4Vectors", ask = FALSE)
	}
	# If present get the scaled abundance
	.abundance_scaled =
		.data %>%
		ifelse_pipe(
			".abundance_scaled" %in% ((.) %>% get_tt_columns() %>% names) &&
				# .data %>% get_tt_columns() %$% .abundance_scaled %>% is.null %>% not() &&
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
		pivot_longer( cols=-c(!!.transcript,!!.sample), names_to="assay", values_to= ".a") %>%
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

#' This function is needed for DE in case the matrix is not rectangular, but includes NA
#'
#' @keywords internal
#' @noRd
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr set_colnames
#' @importFrom stats model.matrix
#' @importFrom stats as.formula
#' @importFrom utils install.packages
#' @importFrom tidyr complete
#' @importFrom rlang quo_is_symbol
#'
#' @param .data A tibble
#' @param .formula a formula with no response variable, of the kind ~ factor_of_interest + batch
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param .abundance_scaled The name of the transcript/gene scaled abundance column
#'
#'
#' @return A tibble with adjusted counts
#'
#'
fill_NA_using_formula = function(.data,
																 .formula,
																 .sample = NULL,
																 .transcript = NULL,
																 .abundance = NULL,
																 .abundance_scaled = NUL,
																 suffix = "_imputed"){

	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	.abundance_scaled = enquo(.abundance_scaled)

	col_formula =
		.data %>%
		select(parse_formula(.formula)) %>%
		distinct() %>%
		select_if(function(x) is.character(x) | is.logical(x) | is.factor(x)) %>%
		colnames

	# Sample-wise columns
	sample_col = .data %>% get_specific_annotation_columns(!!.sample) %>% outersect(col_formula)

	need_log = .data %>% pull(!!.abundance) %>%  max(na.rm=TRUE) > 50

	# Create NAs for missing sample/transcript pair
 .data_completed =
		.data %>%

		# Add missing pairs
 		nest(ct_data = -c(col_formula)) %>%
 		mutate(ct_data = map(ct_data, ~ .x %>% droplevels() %>% complete(!!as.symbol(quo_name(.sample)), !!.transcript) )) %>%
 		unnest(ct_data)

 # For non scaled counts create a pseudo scale based on library size, then calculate imputed and scale back
 abundance_is_int = .data %>% slice(1) %>% pull(!!.abundance) %>% class() %>% equals("integer")
 .data =
   .data %>%
   group_by(!!.sample) %>%
   mutate(library_size__ = sum(!!.abundance)) %>%
   ungroup() %>%
   mutate(!!.abundance := !!.abundance / library_size__)

 imputed_column = sprintf("%s%s", quo_name(.abundance), suffix )
 imputed_column_scaled = sprintf("%s%s", quo_name(.abundance_scaled), suffix )

 # Divide the dataset
 .data_OK =
 	.data %>%
 	anti_join(.data_completed %>% filter(!!.abundance %>% is.na) %>% select( !!.transcript, col_formula) %>% distinct(), by = c(quo_name(.transcript), col_formula)) %>%

   # Add the imputed column
   mutate(!!as.symbol(imputed_column) := !!.abundance) %>%
   when(
     quo_is_symbol(.abundance_scaled) ~ .x %>%
       mutate(!!as.symbol(imputed_column_scaled) := !!.abundance_scaled),
     ~ (.)
   )

 .data_FIXED =
 .data %>%
 	inner_join(.data_completed %>% filter(!!.abundance %>% is.na) %>% select( !!.transcript, col_formula) %>% distinct(), by = c(quo_name(.transcript), col_formula)) %>%

 	# attach NAs
 	bind_rows(
	.data_completed %>%
		filter(!!.abundance %>% is.na) %>%
		select(!!.sample, !!.transcript) %>%
		left_join(.data %>% pivot_sample(!!.sample), by = quo_name(.sample)) %>%
		left_join(.data %>% pivot_transcript(!!.transcript), by = quo_name(.transcript))
	)


 # Clean environment
 rm(.data_completed)

 ~ {

   # Pseudo-scale if not scaled
   if(!grepl("_scaled", .y)) library_size = colSums(.x, na.rm = TRUE)
   if(!grepl("_scaled", .y)) .x = .x / library_size

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
   if(!grepl("_scaled", .y)) .x = .x * library_size

   # Return
   .x
 }

 .data_FIXED =
   .data_FIXED %>%

   when( need_log ~ mutate(., !!.abundance := log1p(!!.abundance)), ~ (.)   ) %>%
   when( need_log & quo_is_symbol(.abundance_scaled) ~ mutate(., !!.abundance_scaled := log1p(!!.abundance_scaled)), ~ (.)   ) %>%


	# Group by covariate
	nest(cov_data = -c(col_formula, !!.transcript)) %>%
	mutate(cov_data = map(cov_data, ~
											.x %>%
											mutate(
											  !!as.symbol(imputed_column) := ifelse(
													!!.abundance %>% is.na,
													median(!!.abundance, na.rm = TRUE),
													!!.abundance
												)
											) %>%

											# Impute scaled if exist
											ifelse_pipe(
												quo_is_symbol(.abundance_scaled),
												~ .x %>% mutate(
												  !!as.symbol(imputed_column_scaled) := ifelse(
														!!.abundance_scaled %>% is.na,
														median(!!.abundance_scaled, na.rm = TRUE),
														!!.abundance_scaled
													)
												)
											) %>%

											# Through warning if group of size 1
											ifelse_pipe((.) %>% nrow %>% `<` (2), warning("tidybulk says: According to your design matrix, u have sample groups of size < 2, so you your dataset could still be sparse."))
	)) %>%
	unnest(cov_data) %>%
   when( need_log ~ mutate(., !!.abundance := exp(!!.abundance)-1), ~ (.)   ) %>%
   when( need_log & quo_is_symbol(.abundance_scaled) ~ mutate(., !!.abundance_scaled := exp(!!.abundance_scaled)-1), ~ (.)   )


	.data_OK %>%
		bind_rows(.data_FIXED) %>%

	  # Scale back the pseudoscaling
	  mutate(!!.abundance := !!.abundance * library_size__) %>%
	  select(-library_size__) %>%
	  when(abundance_is_int ~ mutate(., !!.abundance := as.integer(!!.abundance)), ~ (.)) %>%

		# Reattach internals
		reattach_internals(.data)

}

#' This function is needed for DE in case the matrix is not rectangular, but includes NA
#'
#' @keywords internal
#' @noRd
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr set_colnames
#' @importFrom stats model.matrix
#' @importFrom stats as.formula
#'
#' @param .data A `tbl` formatted as | <element> | <feature> | <value> | <...> |
#' @param .sample The name of the element column
#' @param .transcript The name of the feature/gene column
#' @param .abundance The name of the feature/gene value column
#' @param fill_with A numerical value with which fill the missing data points
#'
#'
#' @return A tibble with adjusted counts
#'
#'
fill_NA_using_value = function(.data,
															 .sample,
															 .transcript,
															 .abundance,
															 fill_with){

	# Comply with CRAN NOTES
	. = NULL

	# Get column names
	.element = enquo(.sample)
	.feature = enquo(.transcript)
	.value = enquo(.abundance)

	# Scale based on library size


	# Create NAs for missing element/feature pair
	df_to_impute =
		.data %>%
		select(!!.element, !!.feature, !!.value) %>%
		distinct %>%
		pivot_wider(
			names_from = !!.feature,
			values_from = !!.value,
			names_sep = "___",
			names_prefix = "fill_miss_"
		) %>%
		pivot_longer(
			names_to = .data %>% select(!!.feature) %>% names,
			values_to = quo_names(.value),
			names_sep = purrr::when(quo_names(.feature), length(.) > 1 ~ "___", ~ NULL),
			names_prefix = "fill_miss_",
			cols = contains("fill_miss_")
		)

	# Select just features/covariates that have missing
	combo_to_impute = df_to_impute %>% anti_join(.data, by=c(quo_names(.element), quo_names(.feature))) %>% select(!!.feature, !!.element) %>% distinct()

	# Impute using median
	df_to_impute %>%
		inner_join(combo_to_impute, by = c(quo_names(.element), quo_names(.feature))) %>%

		# Fill
		mutate(!!.value := if_else(!!.value %>% is.na, fill_with, !!.value)) %>%
		# when(
		# 	quo_is_symbol(.value_scaled) ~ mutate(., !!.value_scaled := !!.value)) ,
		# 	~ (.)
		# ) %>%

		# In next command avoid error if no data to impute
		ifelse_pipe(
			nrow(.) > 0,
			~ .x %>% left_join(.data %>% pivot_sample(!!.element), by=quo_names(.element))
		) %>%

		# Add original dataset
		bind_rows(.data %>% anti_join(combo_to_impute, by=c(quo_names(.feature), quo_names(.element)))) %>%
		select(.data %>% colnames) %>%

		# Reattach internals
		reattach_internals(.data)

}


# # Iterative version of Siberg function because fails
# siberg_iterative = function(x) {
# 	if (x %>% unique %>% length %>% st(5))
# 		return(c(NA, NA))
#
#
#
# 	mu = NA
# 	max_i = ceiling(length(x) / 10)
# 	#max_i = 10
# 	i = 0
# 	while (mu %>% is.na | i <= max_i) {
# 		res = SIBERG::SIBER(x, model = 'NB')
#
# 		BI = res[7]
# 		mu = res[1]
# 		x = x[-1]
# 		i = i + 1
#
# 	}
#
#
# 	if (mu %>% is.na & x %>% length %>% st(5))
# 		return(c(NA, NA))
#
# 	return(c(max(res[1], res[2]) / (min(res[1], res[2]) + 1),
# 					 res[7]))
# }


# # Calculate bimodality
# bimodality =
#
# 	counts_non_red %>%
# 	#keep_variable(top = 5000) %>%
# 	tidybulk:::drop_class(c("tidybulk", "tt")) %>%
# 	tidybulk:::drop_internals() %>%
# 	nest(data = -c(`Cell type formatted`, symbol)) %>%
#
# 	#slice(1:10) %>%
# 	mutate(	bimodality_NB =
# 		map(
# 			data,
# 			~ tryCatch(
# 							.x %>% pull(count_scaled) %>% as.integer %>%
# 								siberg_iterative() %>%
# 								`[` (1:2) , error=function(e) c(NA, NA))		%>%
# 							setNames(c("bimodality_NB_diff", "bimodality_NB")) %>%
# 							enframe() %>% spread(name, value)
#
# 		)
# 	) %>%
# 	select(-data) %>%
# 	unnest(bimodality_NB)
#
# bimodality %>% saveRDS("dev/bimodality.rds")
#
# bimodality = readRDS("dev/bimodality.rds")
#
# non_bimodal =
# 	bimodality %>%
# 	add_count(symbol) %>%
# 	filter(n==max(n)) %>%
# 	mutate(bimodal = ((bimodality_NB > 0.8 & bimodality_NB_diff > 20) | bimodality_NB_diff > 100) ) %>%
# 	nest(data = -symbol) %>%
# 	mutate(how_many_bimod = map_int(data, ~ .x %>% pull(bimodal) %>% sum(na.rm=TRUE))) %>%
# 	filter(how_many_bimod == 0)


#' @keywords internal
#' @noRd
#'
#' @importFrom stats p.adjust
entrez_over_to_gsea = function(my_entrez_rank, species, gene_collections  = NULL){

	# From the page
	# https://yulab-smu.github.io/clusterProfiler-book/chapter5.html

	# Check if package is installed, otherwise install
	if (find.package("fastmatch", quiet = TRUE) %>% length %>% equals(0)) {
		message("tidybulk says: Installing fastmatch needed for analyses")
		install.packages("fastmatch", repos = "https://cloud.r-project.org")
	}

	if (find.package("clusterProfiler", quiet = TRUE) %>% length %>% equals(0)) {
		message("clusterProfiler not installed. Installing.")
		BiocManager::install("clusterProfiler", ask = FALSE)
	}




	# Get gene sets signatures
	msigdbr::msigdbr(species = species) %>%

		# Filter specific gene_collections  if specified. This was introduced to speed up examples executionS
		when(
			!is.null(gene_collections ) ~ filter(., gs_cat %in% gene_collections ),
			~ (.)
		) %>%

		# Execute calculation
		nest(data = -gs_cat) %>%
		mutate(test =
					 	map(
					 		data,
					 		~ clusterProfiler::enricher(
					 			my_entrez_rank,
					 			TERM2GENE=.x %>% select(gs_name, entrez_gene),
					 			pvalueCutoff = 1
					 		) %>%
					 			as_tibble
					 	)) %>%
		select(-data) %>%
		unnest(test) %>%

		# Order
		arrange(`p.adjust`) %>%

		# format transcripts
		mutate(entrez = strsplit(geneID, "/")) %>%
		select(-geneID)

}


#'
#' @keywords internal
#' @noRd
#'
#' @importFrom tibble rowid_to_column
#' @importFrom stats p.adjust
#' @importFrom purrr map
#'
entrez_rank_to_gsea = function(my_entrez_rank, species, gene_collections  = NULL){

	# From the page
	# https://yulab-smu.github.io/clusterProfiler-book/chapter5.html

	# Check if package is installed, otherwise install
	if (find.package("fastmatch", quiet = TRUE) %>% length %>% equals(0)) {
		message("tidybulk says: Installing fastmatch needed for analyses")
		install.packages("fastmatch", repos = "https://cloud.r-project.org")
	}

	if (find.package("clusterProfiler", quiet = TRUE) %>% length %>% equals(0)) {
		message("tidybulk says: clusterProfiler not installed. Installing.")
		BiocManager::install("clusterProfiler", ask = FALSE)
	}

	if (find.package("enrichplot", quiet = TRUE) %>% length %>% equals(0)) {
		message("tidybulk says: enrichplot not installed. tidybulk says: Installing.")
		BiocManager::install("enrichplot", ask = FALSE)
	}

	if (find.package("ggplot2", quiet = TRUE) %>% length %>% equals(0)) {
		message("tidybulk says: Installing ggplot2 needed for analyses")
		install.packages("ggplot2", repos = "https://cloud.r-project.org")
	}

  # Get gene sets signatures
  my_gene_collection =
    gene_collections %>%
    when(
      is.null(gene_collections ) ~  msigdbr::msigdbr(species = species) ,
      is(., "character") ~  msigdbr::msigdbr(species = species) %>%  filter( gs_cat %in% gene_collections ),
      is(., "list") ~
        tibble(gs_name=names(.), entrez_gene = . ) %>% unnest(entrez_gene) %>% mutate(gs_cat = "user_defined"),
      ~ stop("tidybulk says: the gene sets should be either a character vector or a named list")
    )
   my_gene_collection %>%


		# Execute calculation
		nest(data = -gs_cat) %>%
		mutate(fit =
					 	map(
					 		data,
					 		~ 	clusterProfiler::GSEA(
					 				my_entrez_rank,
					 				TERM2GENE=.x %>% select(gs_name, entrez_gene),
					 				pvalueCutoff = 1
					 		)

					 	)) %>%
			mutate(test =
					 	map(
					 		fit,
					 		~ .x %>%
					 			ggplot2::fortify(showCategory=Inf) %>%
					 			as_tibble() %>%
					 			rowid_to_column(var = "idx_for_plotting")
					 			#%>%
					 			#	mutate(plot = future_imap(ID, ~ enrichplot::gseaplot2(fit, geneSetID = .y, title = .x)))

					 	)) %>%
		select(-data)


}

# gsea_de = function(.data,
# 								 .formula,
# 								 .sample = NULL,
# 								 .entrez,
# 								 .abundance = NULL,
# 								 .contrasts = NULL,
# 								 species){
#
# 	# Comply with CRAN NOTES
# 	. = NULL
#
# 	# Check packages
# 	# Check if package is installed, otherwise install
# 	if (find.package("EGSEA", quiet = TRUE) %>% length %>% equals(0)) {
# 		message("EGSEA not installed. Please install it with.")
# 		message("BiocManager::install(\"EGSEA\", ask = FALSE)")
# 	}
# 	if (!"EGSEA" %in% (.packages())) {
# 		message("EGSEA package not loaded. Please run library(\"EGSEA\")")
# 	}
#
# 	# Check column type
# 	if (reference %>% sapply(class) %in% c("numeric", "double", "integer") %>% not() %>% any)
# 		stop("tidybulk says: your reference has non-numeric/integer columns.")
#
# 	# Get column names
# 	.sample = enquo(.sample)
# 	.abundance = enquo(.abundance)
# 	col_names = get_sample_counts(.data, .sample, .abundance)
# 	.sample = col_names$.sample
# 	.abundance = col_names$.abundance
#
#
# 	df_for_edgeR %>%
# 		select(!!.sample, !!.entrez,!!.abundance)
#
#
# 	library(msigdbr)
# 	library(clusterProfiler)
#
# 	counts %>% tidybulk(sample, transcript, count) %>% test_differential_abundance(~ condition) %>% arrange(logFC %>% desc) %>% pull(transcript)
#
#
# 	em <- enricher(gene, TERM2GENE=m_df %>% select(gs_name, entrez_gene))
#
# 	emGSEA(geneList, TERM2GENE = m_t2g)
#
#
# 	m_t2g <- msigdbr(species = "Homo sapiens", category = "C6")
#
#
#
# }
