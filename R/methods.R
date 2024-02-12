setOldClass("tidybulk")

#' Creates an annotated `tidybulk` tibble from a `tbl` or 
#' `SummarizedExperiment` object
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description tidybulk() creates an annotated `tidybulk` tibble from a 
#'   `tbl` (with at least three columns for sample, feature and transcript 
#'   abundance) or `SummarizedExperiment` (more convenient if abstracted to 
#'   tibble with library(tidySummarizedExperiment))
#'
#' @importFrom rlang enquo
#' @importFrom rlang quo_is_missing
#' @import readr
#' @import SummarizedExperiment
#' @import methods
#'
#' @name tidybulk
#'
#' @param .data A `tbl` (with at least three columns for sample, feature 
#'   and transcript abundance) or `SummarizedExperiment` (more convenient if 
#'   abstracted to tibble with library(tidySummarizedExperiment))
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param .abundance_scaled The name of the transcript/gene scaled 
#'   abundance column
#'
#' @details This function creates a tidybulk object and is useful if you want
#' to avoid to specify .sample, .transcript and .abundance arguments all the times.
#' The tidybulk object have an attribute called internals where these three
#' arguments are stored as metadata. They can be extracted as 
#' attr(<object>, "internals").
#'
#' @return A `tidybulk` object
#'
#'
#' @examples
#'
#' tidybulk(tidybulk::se_mini)
#'
#'
#' @docType methods
#' @rdname tidybulk-methods
#'
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

	if(
	  quo_is_missing(.sample) |
	  quo_is_missing(.transcript) |
	  quo_is_missing(.abundance)
	) stop("tidybulk says: the arguments .sample, .transcript and .abundance ",
	       "must include column names (not surrounded by quotes)")

	# Validate data frame
	if(do_validate()) validation(.data,
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
setMethod("tidybulk", "spec_tbl_df", .tidybulk)

#' tidybulk
#'
#' @export
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


#' as_SummarizedExperiment
#'
#' @description as_SummarizedExperiment() creates a `SummarizedExperiment` 
#'   object from a `tbl` or `tidybulk` tbl formatted as | <SAMPLE> | 
#'   <TRANSCRIPT> | <COUNT> | <...> |
#'
#'
#' @importFrom utils data
#' @importFrom tidyr pivot_longer
#'
#' @param .data A tibble
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#'
#' @return A `SummarizedExperiment` object
#'
#' @docType methods
#' @rdname as_SummarizedExperiment-methods
#' @export
#'
setGeneric("as_SummarizedExperiment", function(.data,
																							 .sample = NULL,
																							 .transcript = NULL,
																							 .abundance = NULL)
	standardGeneric("as_SummarizedExperiment"))


.as_SummarizedExperiment = function(.data,
																		.sample = NULL,
																		.transcript = NULL,
																		.abundance = NULL) {

  # Fix NOTEs
  . = NULL

	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	col_names = get_sample_transcript_counts(.data, .sample, 
	                                         .transcript, .abundance)
	.sample = col_names$.sample
	.transcript = col_names$.transcript
	.abundance = col_names$.abundance

	# Check if package is installed, otherwise install
	if (find.package("SummarizedExperiment", quiet = TRUE) |> length() |> equals(0)) {
		message("Installing SummarizedExperiment")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("SummarizedExperiment", ask = FALSE)
	}
	if (find.package("S4Vectors", quiet = TRUE) |> length() %>% equals(0)) {
		message("Installing S4Vectors")
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

	colData =
	  .data %>%
	  select(!!.sample, sample_cols) %>%
	  distinct() %>%

	  # Unite if multiple sample columns
	  tidyr::unite(!!sample__$name, !!.sample, remove = FALSE, sep = "___") |>

	  arrange(!!sample__$symbol) %>% {
		S4Vectors::DataFrame(
		  (.) %>% select(-!!sample__$symbol),
		  row.names = (.) %>% pull(!!sample__$symbol)
		)
	}

	rowData =
	  .data %>%
	  select(!!.transcript, feature_cols) %>%
	  distinct() %>%

	  # Unite if multiple sample columns
	  tidyr::unite(!!feature__$name, !!.transcript, remove = FALSE, sep = "___") |>

	  arrange(!!feature__$symbol) %>% {
		S4Vectors::DataFrame(
		  (.) %>% select(-!!feature__$symbol),
		  row.names = (.) %>% pull(!!feature__$symbol)
		)
	}

	my_assays =
		.data %>%

	  # Unite if multiple sample columns
	  tidyr::unite(!!sample__$name, !!.sample, remove = FALSE, sep = "___") |>

	  # Unite if multiple sample columns
	  tidyr::unite(!!feature__$name, !!.transcript, remove = FALSE, sep = "___") |>

	  select(!!sample__$symbol,
	         !!feature__$symbol,
	         !!.abundance,
	         !!.abundance_scaled,
	         counts_cols) %>%
	  distinct() %>%

		pivot_longer( cols=-c(!!feature__$symbol,!!sample__$symbol), 
		              names_to="assay", values_to= ".a") %>%
		nest(`data` = -`assay`) %>%
		mutate(`data` = `data` %>%  map(
			~ .x %>%
			  spread(!!sample__$symbol, .a) %>%

			  # arrange sample
			  select(!!feature__$symbol, rownames(colData)) |>

			  # Arrange symbol
			  arrange(!!feature__$symbol) |>

			  # Convert
			  as_matrix(rownames = feature__$name)
		))

	# Build the object
	SummarizedExperiment::SummarizedExperiment(
		assays = my_assays %>% pull(`data`) %>% setNames(my_assays$assay),
		rowData = rowData,
		colData = colData
	)

}

#' as_SummarizedExperiment
#'
#' @export
#'
#' @inheritParams as_SummarizedExperiment
#'
#' @docType methods
#' @rdname as_SummarizedExperiment-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("as_SummarizedExperiment", "spec_tbl_df", .as_SummarizedExperiment)

#' as_SummarizedExperiment
#'
#' @export
#'
#' @inheritParams as_SummarizedExperiment
#'
#' @docType methods
#' @rdname as_SummarizedExperiment-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("as_SummarizedExperiment", "tbl_df", .as_SummarizedExperiment)

#' as_SummarizedExperiment
#'
#' @export
#'
#' @inheritParams as_SummarizedExperiment
#'
#' @docType methods
#' @rdname as_SummarizedExperiment-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("as_SummarizedExperiment", "tidybulk", .as_SummarizedExperiment)


#' Creates a `tt` object from a list of file names of BAM/SAM
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description tidybulk_SAM_BAM() creates a `tt` object from A `tbl` 
#'   (with at least three columns for sample, feature and transcript abundance) 
#'   or `SummarizedExperiment` (more convenient if abstracted to tibble with 
#'   library(tidySummarizedExperiment))
#'
#' @importFrom rlang enquo
#'
#' @name tidybulk_SAM_BAM
#'
#' @param file_names A character vector
#' @param genome A character string specifying an in-built annotation used for 
#'   read summarization. It has four possible values including "mm10", "mm9", 
#'   "hg38" and "hg19"
#' @param ... Further parameters passed to the function Rsubread::featureCounts
#'
#' @details This function is based on FeatureCounts package 
#'   (DOI: 10.1093/bioinformatics/btt656). This function creates a tidybulk 
#'   object and is useful if you want to avoid to specify .sample, 
#'   .transcript and .abundance arguments all the times. The tidybulk object 
#'   have an attribute called internals where these three arguments are stored 
#'   as metadata. They can be extracted as attr(<object>, "internals").
#'
#' Underlying core function
#' Rsubread::featureCounts(annot.inbuilt = genome,nthreads = n_cores, ...)
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
#'
#' @export
#'
#' @inheritParams tidybulk_SAM_BAM-methods
#'
#' @docType methods
#' @rdname tidybulk_SAM_BAM-methods
#'
#' @return A `tidybulk` object
#'
setMethod("tidybulk_SAM_BAM", c(file_names = "character", genome = "character"),
          function(file_names, genome = "hg38", ...)
	create_tt_from_bam_sam_bulk(file_names = file_names, genome = genome, ...))

#' Scale the counts of transcripts/genes
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description scale_abundance() takes as input A `tbl` (with at least three 
#'   columns for sample, feature and transcript abundance) or 
#'   `SummarizedExperiment` (more convenient if abstracted to tibble with 
#'   library(tidySummarizedExperiment)) and Scales transcript abundance 
#'   compansating for sequencing depth (e.g., with TMM algorithm, 
#'   Robinson and Oshlack doi.org/10.1186/gb-2010-11-3-r25).
#'
#' @importFrom rlang enquo
#' @importFrom stats median
#'
#' @name scale_abundance
#'
#' @param .data A `tbl` (with at least three columns for sample, 
#'   feature and transcript abundance) or `SummarizedExperiment` 
#'   (more convenient if abstracted to tibble with 
#'   library(tidySummarizedExperiment))
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param method A character string. The scaling method passed to the back-end 
#'   function (i.e., edgeR::calcNormFactors; "TMM","TMMwsp","RLE","upperquartile")
#' @param reference_sample A character string. The name of the reference sample.
#'   If NULL the sample with highest total read count will be selected as reference.
#' @param .subset_for_scaling A gene-wise quosure condition. This will be used 
#'  to filter rows (features/genes) of the dataset. For example
#' @param action A character string between "add" (default) and "only". 
#'  "add" joins the new information to the input tbl (default), "only" 
#'  return a non-redundant tbl with the just new information.
#'
#' @param reference_selection_function DEPRECATED. please use reference_sample.
#'
#' @details Scales transcript abundance compensating for sequencing depth
#' (e.g., with TMM algorithm, Robinson and Oshlack doi.org/10.1186/gb-2010-11-3-r25).
#' Lowly transcribed transcripts/genes (defined with minimum_counts 
#' and minimum_proportion parameters) are filtered out from the scaling procedure.
#' The scaling inference is then applied back to all unfiltered data.
#'
#' Underlying method
#' edgeR::calcNormFactors(.data, method = c("TMM","TMMwsp","RLE","upperquartile"))
#'
#'
#' @return A tbl object with additional columns with scaled data 
#' as `<NAME OF COUNT COLUMN>_scaled`
#'
#' @examples
#'  tidybulk::se_mini |>
#'    identify_abundant() |>
#'    scale_abundance()
#'
#' @docType methods
#' @rdname scale_abundance-methods
#' @export

setGeneric("scale_abundance", function(.data,
																			 .sample = NULL,
																			 .transcript = NULL,
																			 .abundance = NULL,
																			 method = "TMM",
																			 reference_sample = NULL,
																			 .subset_for_scaling = NULL,
																			 action = "add",
																			 # DEPRECATED
																			 reference_selection_function = NULL)
	standardGeneric("scale_abundance"))

# Set internal
.scale_abundance = 	function(.data,
														 .sample = NULL,
														 .transcript = NULL,
														 .abundance = NULL,
														 method = "TMM",
														 reference_sample = NULL,
														 .subset_for_scaling = NULL,
														 action = "add",

														 # DEPRECATED
														 reference_selection_function = NULL) {

  # Fix NOTEs
  . = NULL

	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	col_names = get_sample_transcript_counts(.data, .sample, 
	                                         .transcript, .abundance)
	.sample = col_names$.sample
	.transcript = col_names$.transcript
	.abundance = col_names$.abundance

	.subset_for_scaling = enquo(.subset_for_scaling)

	# Set column name for value scaled
	value_scaled = as.symbol(sprintf("%s%s",  quo_name(.abundance), scaled_string))

	# DEPRECATION OF reference function
	if (is_present(reference_selection_function) & !is.null(reference_selection_function)) {
		# Signal the deprecation to the user
		deprecate_warn("1.1.8", "tidybulk::scale_abundance(reference_selection_function = )", 
		  details = "The argument reference_selection_function is now deprecated please use reference_sample. By default the reference selection function is max()")

	}

	# Validate data frame
	if(do_validate()) {
		validation(.data, !!.sample, !!.transcript, !!.abundance)
		warning_if_data_is_not_rectangular(.data, !!.sample, !!.transcript, !!.abundance)
	}

	# Check that reference sample exists
	if(!is.null(reference_sample) && !reference_sample %in% (.data %>% pull(!!.sample)))
		stop("tidybulk says: your reference sample is not among the samples in your data frame")

	.data_norm =
		.data %>%

		# Filter abundant if performed
		when(
			".abundant" %in% colnames(.) ~ filter(., .abundant),
			~ {
				warning("tidybulk says: highly abundant transcripts were not identified",
				" (i.e. identify_abundant()) or filtered (i.e., keep_abundant), ",
				"therefore this operation will be performed on unfiltered data. ",
				"In rare occasions this could be wanted. In standard whole-transcriptome ",
				"workflows is generally unwanted.")
				(.)
			}
		) %>%

	  # filter based on user selection
	  when(
	   !quo_is_null(.subset_for_scaling) ~ filter(., !!.subset_for_scaling),
	   ~ (.)
	 ) %>%

	  # Check I have genes left
	  when(nrow(.) == 0 ~ stop("tidybulk says: there are 0 genes that passes ",
	  "the filters (.abundant and/or .subset_for_scaling). ",
	  "Please check your filtering or your data."), ~ (.)) %>%

		get_scaled_counts_bulk(
			.sample = !!.sample,
			.transcript = !!.transcript,
			.abundance = !!.abundance,
			method = method,
			reference_sample = reference_sample
		) %>%

		# Attach column internals
		add_tt_columns(
			!!.sample,
			!!.transcript,
			!!.abundance,
			!!(function(x, v)	enquo(v))(x,!!value_scaled)
		)


	if (action == "add"){

		.data %>%

			left_join(.data_norm, by=quo_name(.sample)) %>%
			dplyr::mutate(!!value_scaled := !!.abundance * multiplier) %>%

			# Attach attributes
			reattach_internals(.data_norm)

	}
	else if (action == "get"){

		.data %>%

			# Selecting the right columns
			pivot_sample(!!.sample) %>%

			# Join result
			left_join(.data_norm, by=quo_name(.sample)) %>%


			# Attach attributes
			reattach_internals(.data_norm)
	}
	else if (action == "only") .data_norm
	else
		stop(
			"tidybulk says: action must be either \"add\" for adding this ",
			"information to your data frame or \"get\" to just get the information"
		)
}


#' scale_abundance
#'
#' @export
#' @inheritParams scale_abundance
#' @docType methods
#' @rdname scale_abundance-methods
#'
#' @return A tbl object with additional columns with scaled data 
#' as `<NAME OF COUNT COLUMN>_scaled`
#'
setMethod("scale_abundance", "spec_tbl_df", .scale_abundance)

#' scale_abundance
#'
#' @export
#' @inheritParams scale_abundance
#' @docType methods
#' @rdname scale_abundance-methods
#'
#' @return A tbl object with additional columns with scaled data 
#' as `<NAME OF COUNT COLUMN>_scaled`
#'
setMethod("scale_abundance", "tbl_df", .scale_abundance)

#' scale_abundance
#'
#' @export
#'
#' @inheritParams scale_abundance
#'
#' @docType methods
#' @rdname scale_abundance-methods
#'
#' @return A tbl object with additional columns with scaled data 
#' as `<NAME OF COUNT COLUMN>_scaled`
#'
setMethod("scale_abundance", "tidybulk", .scale_abundance)



#' Normalise by quantiles the counts of transcripts/genes
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description quantile_normalise_abundance() takes as input A `tbl`
#'   (with at least three columns for sample, feature and transcript abundance) 
#'   or `SummarizedExperiment` (more convenient if abstracted to tibble with 
#'   library(tidySummarizedExperiment)) and Scales transcript abundance 
#'   compansating for sequencing depth (e.g., with TMM algorithm, Robinson 
#'   and Oshlack doi.org/10.1186/gb-2010-11-3-r25).
#'
#' @importFrom rlang enquo
#' @importFrom stats median
#' @importFrom dplyr join_by
#'
#' @name quantile_normalise_abundance
#'
#' @param .data A `tbl` (with at least three columns for sample, 
#'   feature and transcript abundance) or `SummarizedExperiment` 
#'   (more convenient if abstracted to tibble with 
#'   library(tidySummarizedExperiment))
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param method A character string. Either "limma_normalize_quantiles" 
#'   for limma::normalizeQuantiles or 
#'   "preprocesscore_normalize_quantiles_use_target" for 
#'   preprocessCore::normalize.quantiles.use.target for large-scale dataset, 
#'   where limmma could not be compatible.
#' @param action A character string between "add" (default) and "only". 
#'   "add" joins the new information to the input tbl (default), 
#'   "only" return a non-redundant tbl with the just new information.
#'
#'
#' @details Scales transcript abundance compensating for sequencing depth
#' (e.g., with TMM algorithm, Robinson and Oshlack doi.org/10.1186/gb-2010-11-3-r25).
#' Lowly transcribed transcripts/genes (defined with minimum_counts and 
#' minimum_proportion parameters) are filtered out from the scaling procedure.
#' The scaling inference is then applied back to all unfiltered data.
#'
#' Underlying method
#' edgeR::calcNormFactors(.data, method = c("TMM","TMMwsp","RLE","upperquartile"))
#'
#'
#'
#' @return A tbl object with additional columns with scaled data as 
#'   `<NAME OF COUNT COLUMN>_scaled`
#'
#'
#' @examples
#'
#'
#'  tidybulk::se_mini |>
#'    quantile_normalise_abundance()
#'
#'
#'
#' @docType methods
#' @rdname quantile_normalise_abundance-methods
#' @export

setGeneric("quantile_normalise_abundance", function(
                                          .data,
                                          .sample = NULL,
                                          .transcript = NULL,
                                          .abundance = NULL,
                                          method = "limma_normalize_quantiles",
                                          action = "add")
  standardGeneric("quantile_normalise_abundance"))

# Set internal
.quantile_normalise_abundance <- function(.data,
                                          .sample = NULL,
                                          .transcript = NULL,
                                          .abundance = NULL,
                                          method = "limma_normalize_quantiles",
                                          action = "add") {

  # Fix NOTEs
  . = NULL

  # Get column names
  .sample = enquo(.sample)
  .transcript = enquo(.transcript)
  .abundance = enquo(.abundance)
  col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
  .sample = col_names$.sample
  .transcript = col_names$.transcript
  .abundance = col_names$.abundance

  # Set column name for value scaled
  value_scaled = as.symbol(sprintf("%s%s",  quo_name(.abundance), scaled_string))

  # Reformat input data set
  .data_norm <-
    .data %>%

    # Rename
    dplyr::select(!!.sample,!!.transcript,!!.abundance) %>%

    # Set samples and genes as factors
    dplyr::mutate(!!.sample := factor(!!.sample),!!.transcript := factor(!!.transcript))  |>
    pivot_wider(names_from = !!.sample, values_from = !!.abundance) |>
    as_matrix(rownames=!!.transcript)


  if(tolower(method) == "limma_normalize_quantiles"){

    # Check if package is installed, otherwise install
    if (find.package("limma", quiet = TRUE) %>% length %>% equals(0)) {
      message("tidybulk says: Installing limma needed for analyses")
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install("limma", ask = FALSE)
    }

    .data_norm =
      .data_norm |>
      limma::normalizeQuantiles()
  }
  else if(tolower(method) == "preprocesscore_normalize_quantiles_use_target"){

    # Check if package is installed, otherwise install
    if (find.package("preprocessCore", quiet = TRUE) %>% length %>% equals(0)) {
      message("tidybulk says: Installing preprocessCore needed for analyses")
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install("preprocessCore", ask = FALSE)
    }

    .data_norm_quant =
      .data_norm |>
      preprocessCore::normalize.quantiles.use.target(
        target = preprocessCore::normalize.quantiles.determine.target(.data_norm)
      )

    colnames(.data_norm_quant) = .data_norm |> colnames()
    rownames(.data_norm_quant) = .data_norm |> rownames()

    .data_norm = .data_norm_quant
    rm(.data_norm_quant)

  } else stop("tidybulk says: the methods must be limma_normalize_quantiles or preprocesscore")

  .data_norm =
    .data_norm |>
    as_tibble(rownames = quo_name(.transcript)) |>
    pivot_longer(-!!.transcript, names_to = quo_name(.sample), values_to = quo_names(value_scaled)) |>


    # Attach column internals
    add_tt_columns(
      !!.sample,
      !!.transcript,
      !!.abundance,
      !!(function(x, v)	enquo(v))(x,!!value_scaled)
    )


  if (action %in% c( "add", "get")){

    .data %>%

      left_join(.data_norm, by= join_by(!!.sample, !!.transcript)) %>%

      # Attach attributes
      reattach_internals(.data_norm) |>

      # Add methods
      memorise_methods_used(c("quantile"))

  }
  else if (action == "only") .data_norm
  else
    stop(
      "tidybulk says: action must be either \"add\" for adding this ",
      "information to your data frame or \"get\" to just get the information"
    )
}

#' quantile_normalise_abundance
#'
#' @export
#'
#' @inheritParams quantile_normalise_abundance
#'
#' @docType methods
#' @rdname quantile_normalise_abundance-methods
#'
#' @return A tbl object with additional columns with scaled data 
#' as `<NAME OF COUNT COLUMN>_scaled`
#'
setMethod("quantile_normalise_abundance", "spec_tbl_df", 
          .quantile_normalise_abundance)

#' quantile_normalise_abundance
#'
#' @export
#'
#' @inheritParams quantile_normalise_abundance
#'
#' @docType methods
#' @rdname quantile_normalise_abundance-methods
#'
#' @return A tbl object with additional columns with scaled data 
#'  as `<NAME OF COUNT COLUMN>_scaled`
#'
setMethod("quantile_normalise_abundance", "tbl_df", .quantile_normalise_abundance)

#' quantile_normalise_abundance
#'
#' @export
#'
#' @inheritParams quantile_normalise_abundance
#'
#' @docType methods
#' @rdname quantile_normalise_abundance-methods
#'
#' @return A tbl object with additional columns with scaled data 
#' as `<NAME OF COUNT COLUMN>_scaled`
#'
setMethod("quantile_normalise_abundance", "tidybulk", .quantile_normalise_abundance)


#' Get clusters of elements (e.g., samples or transcripts)
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description cluster_elements() takes as input A `tbl` (with at least 
#'  three columns for sample, feature and transcript abundance) or 
#'  `SummarizedExperiment` (more convenient if abstracted to tibble with 
#'  library(tidySummarizedExperiment)) and identify clusters in the data.
#'
#' @importFrom rlang enquo
#'
#' @name cluster_elements
#'
#' @param .data A `tbl` (with at least three columns for sample, 
#'  feature and transcript abundance) or `SummarizedExperiment` 
#'  (more convenient if abstracted to tibble with 
#'  library(tidySummarizedExperiment))
#' @param .element The name of the element column (normally samples).
#' @param .feature The name of the feature column (normally transcripts/genes)
#' @param .abundance The name of the column including the numerical value the 
#'  clustering is based on (normally transcript abundance)
#'
#' @param method A character string. The cluster algorithm to use, at the 
#'  moment k-means is the only algorithm included.
#' @param of_samples A boolean. In case the input is a tidybulk object, 
#'  it indicates Whether the element column will be sample or transcript column
#' @param transform A function that will tranform the counts, by default it is 
#'  log1p for RNA sequencing data, but for avoinding tranformation you can 
#'  use identity
#' @param action A character string. Whether to join the new information to 
#'  the input tbl (add), or just get the non-redundant tbl with the new information (get).
#' @param ... Further parameters passed to the function kmeans
#'
#' @param log_transform DEPRECATED - A boolean, whether the value should be 
#' log-transformed (e.g., TRUE for RNA sequencing data)
#'
#' @details identifies clusters in the data, normally of samples.
#' This function returns a tibble with additional columns for the cluster annotation.
#' At the moment only k-means (DOI: 10.2307/2346830) and SNN clustering 
#' (DOI:10.1016/j.cell.2019.05.031) is supported, the plan is to introduce more 
#' clustering methods.
#'
#' Underlying method for kmeans
#' do.call(kmeans(.data, iter.max = 1000, ...)
#'
#' Underlying method for SNN
#' .data %>%
#' Seurat::CreateSeuratObject() %>%
#' Seurat::ScaleData(display.progress = TRUE,num.cores = 4, do.par = TRUE) %>%
#' Seurat::FindVariableFeatures(selection.method = "vst") %>%
#' Seurat::RunPCA(npcs = 30) %>%
#' Seurat::FindNeighbors() %>%
#' Seurat::FindClusters(method = "igraph", ...)
#'
#'
#' @return A tbl object with additional columns with cluster labels
#'
#'
#' @examples
#'
#'
#'     cluster_elements(tidybulk::se_mini,	centers = 2, method="kmeans")
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
																				transform = log1p,
																				action = "add",
																				...,
																				# DEPRECATED
																				log_transform = NULL
                                      )
	standardGeneric("cluster_elements"))

# Set internal
.cluster_elements <- function(.data,
															.element = NULL,
															.feature = NULL,
															.abundance = NULL,
															method ,
															of_samples = TRUE,
															transform = log1p,
															action = "add",
															...,
															# DEPRECATED
															log_transform = NULL) {

  # Fix NOTEs
  . = NULL

  # DEPRECATION OF log_transform
  if (is_present(log_transform) & !is.null(log_transform)) {

    # Signal the deprecation to the user
    deprecate_warn("1.7.4", 
      "tidybulk::test_differential_abundance(log_transform = )", 
      details = "The argument log_transform is now deprecated, please use transform.")

    if(log_transform == TRUE) transform = log1p
  }

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
	if(do_validate()) {
	validation(.data, !!.element, !!.feature, !!.abundance)
	error_if_data_is_not_rectangular(.data, !!.element, !!.feature, !!.abundance)
	}


	.data_procesed =

		.data %>%

		# Filter abundant if performed
		when(
			".abundant" %in% colnames(.) ~ filter(., .abundant),
			~ {
				warning("tidybulk says: highly abundant transcripts were not ",
				"identified (i.e. identify_abundant()) or filtered (i.e., keep_abundant)",
				", therefore this operation will be performed on unfiltered data. ",
				"In rare occasions this could be wanted. In standard ",
				"whole-transcriptome workflows is generally unwanted.")
				(.)
			}
		) %>%

		# Choose algorithm
		when(
			method == "kmeans" ~ 	get_clusters_kmeans_bulk(.,
				.abundance = !!.abundance,
				.element = !!.element,
				.feature = !!.feature,
				of_samples = of_samples,
				transform = transform,
				...
			),
			method == "SNN" ~ stop("tidybulk says: Matrix package (v1.3-3) causes ",
			  "an error with Seurat::FindNeighbors used in this method. We are ",
			  "trying to solve this issue. At the moment this option in unaviable."),
			# 	get_clusters_SNN_bulk(.,
			# 	.abundance = !!.abundance,
			# 	.element = !!.element,
			# 	.feature = !!.feature,
			# 	of_samples = of_samples,
			# 	transform = transform,
			# 	...
			# ),
			TRUE ~ stop("tidybulk says: the only supported methods are \"kmeans\" or \"SNN\" ")

		)


	# Actions
		if (action == "add"){

			.data |>
				dplyr::left_join(	.data_procesed,	by=quo_name(.element)	) |>

				# Attach attributes
				reattach_internals(.data)

		}
		else if (action == "get"){

			.data |>

				# Selecting the right columns
				pivot_sample(!!.element) |>

				dplyr::left_join(	.data_procesed,	by=quo_name(.element)	) |>

				# Attach attributes
				reattach_internals(.data)

		}
		else if (action == "only") 	.data_procesed
		else
			stop(
				"tidybulk says: action must be either \"add\" for adding this ",
				"information to your data frame or \"get\" to just get the information"
			)

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


#' Dimension reduction of the transcript abundance data
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description reduce_dimensions() takes as input A `tbl` (with at least 
#'   three columns for sample, feature and transcript abundance) or 
#'   `SummarizedExperiment` (more convenient if abstracted to tibble with 
#'   library(tidySummarizedExperiment)) and calculates the reduced dimensional 
#'   space of the transcript abundance.
#'
#' @importFrom rlang enquo
#'
#' @name reduce_dimensions
#'
#' @param .data A `tbl` (with at least three columns for sample, 
#'  feature and transcript abundance) or `SummarizedExperiment` (more 
#'  convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .element The name of the element column (normally samples).
#' @param .feature The name of the feature column (normally transcripts/genes)
#' @param .abundance The name of the column including the numerical value 
#'  the clustering is based on (normally transcript abundance)
#'
#' @param method A character string. The dimension reduction algorithm to 
#'    use (PCA, MDS, tSNE).
#' @param top An integer. How many top genes to select for dimensionality reduction
#' @param of_samples A boolean. In case the input is a tidybulk object, 
#'    it indicates Whether the element column will be sample or transcript column
#' @param .dims An integer. The number of dimensions your are interested in 
#'    (e.g., 4 for returning the first four principal components).
#' @param transform A function that will tranform the counts, by default 
#'  it is log1p for RNA sequencing data, but for avoinding tranformation you 
#'  can use identity
#' @param scale A boolean for method="PCA", this will be passed to the `prcomp`
#'   function. It is not included in the ... argument because although the 
#'   default for `prcomp` if FALSE, it is advisable to set it as TRUE.
#' @param action A character string. Whether to join the new information to 
#'   the input tbl (add), or just get the non-redundant tbl with the new 
#'   information (get).
#' @param ... Further parameters passed to the function prcomp if you choose 
#'   method="PCA" or Rtsne if you choose method="tSNE"
#'
#' @param log_transform DEPRECATED - A boolean, whether the value should be 
#'   log-transformed (e.g., TRUE for RNA sequencing data)
#'
#' @details This function reduces the dimensions of the transcript abundances.
#' It can use multi-dimensional scaling (MDS; DOI.org/10.1186/gb-2010-11-3-r25),
#' principal component analysis (PCA), or tSNE (Jesse Krijthe et al. 2018)
#'
#' Underlying method for PCA:
#' prcomp(scale = scale, ...)
#'
#' Underlying method for MDS:
#' limma::plotMDS(ndim = .dims, plot = FALSE, top = top)
#'
#' Underlying method for tSNE:
#' Rtsne::Rtsne(data, ...)
#'
#' Underlying method for UMAP:
#'
#'  df_source =
#' .data |>
#'
#'   # Filter NA symbol
#'   filter(!!.feature |> is.na() |> not()) |>
#'
#'   # Prepare data frame
#'   distinct(!!.feature,!!.element,!!.abundance) |>
#'
#'   # Filter most variable genes
#'   keep_variable_transcripts(top) |>
#'   reduce_dimensions(method="PCA",  .dims = calculate_for_pca_dimensions,  action="get" ) |>
#'   as_matrix(rownames = quo_name(.element)) |>
#'   uwot::tumap(...)
#'
#'
#' @return A tbl object with additional columns for the reduced dimensions
#'
#'
#' @examples
#'
#'
#'
#' counts.MDS =
#'  tidybulk::se_mini |>
#'  identify_abundant() |>
#'  reduce_dimensions( method="MDS", .dims = 3)
#'
#'
#' counts.PCA =
#'  tidybulk::se_mini |>
#'  identify_abundant() |>
#'  reduce_dimensions(method="PCA", .dims = 3)
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
																				 transform = log1p,
																				 scale = TRUE,
																				 action = "add",
																				 ...,
																				 # DEPRECATED
																				 log_transform = NULL
                                      )
					 standardGeneric("reduce_dimensions"))

# Set internal
.reduce_dimensions <- function(.data,
															 .element = NULL,
															 .feature = NULL,
															 .abundance = NULL,
															 method,
															 .dims = 2,
															 top = 500,
															 of_samples = TRUE,
															 transform = log1p,
															 scale = TRUE,
															 action = "add",
															 ...,
															 # DEPRECATED
															 log_transform = NULL) {

  # Fix NOTEs
  . = NULL
  
  # DEPRECATION OF log_transform
  if (is_present(log_transform) & !is.null(log_transform)) {

    # Signal the deprecation to the user
    deprecate_warn("1.7.4", 
      "tidybulk::test_differential_abundance(log_transform = )", 
      details = "The argument log_transform is now deprecated, please use transform.")

    if(log_transform == TRUE) transform = log1p
  }

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

	# adjust top for the max number of features I have
	if(top > .data |> distinct(!!.feature) |> nrow()){
	  warning(sprintf(
	    "tidybulk says: the \"top\" argument %s is higher than the number of features %s", 
	    top, 
	    .data |> distinct(!!.feature) |> nrow()
	  ))
	  
	  top = min(top, .data |> distinct(!!.feature) |> nrow())
	}
	
	# Validate data frame
	if(do_validate()) {
	validation(.data, !!.element, !!.feature, !!.abundance)
	warning_if_data_is_not_rectangular(.data, !!.element, !!.feature, !!.abundance)
	if(!check_if_transcript_is_na(.data, !!.feature)) stop("tidybulk says: you have empty transcript names")
	}

	.data_processed =

		.data %>%

		# Filter abundant if performed
		when(
			".abundant" %in% colnames(.) ~ filter(., .abundant),
			~ {
				warning("tidybulk says: highly abundant transcripts were not ",
				"identified (i.e. identify_abundant()) or filtered (i.e., keep_abundant), ",
				"therefore this operation will be performed on unfiltered data. ",
				"In rare occasions this could be wanted. In standard whole-transcriptome ",
				"workflows is generally unwanted.")
				(.) }
		) %>%
		when(
			tolower(method) == tolower("MDS") ~ get_reduced_dimensions_MDS_bulk(.,
				.abundance = !!.abundance,
				.dims = .dims,
				.element = !!.element,
				.feature = !!.feature,
				top = top,
				of_samples = of_samples,
				transform = transform,
				...
			),
			tolower(method) == tolower("PCA") ~ get_reduced_dimensions_PCA_bulk(.,
				.abundance = !!.abundance,
				.dims = .dims,
				.element = !!.element,
				.feature = !!.feature,
				top = top,
				of_samples = of_samples,
				transform = transform,
				scale = scale,
				...
			),
			tolower(method) == tolower("tSNE") ~ 	get_reduced_dimensions_TSNE_bulk(.,
				.abundance = !!.abundance,
				.dims = .dims,
				.element = !!.element,
				.feature = !!.feature,
				top = top,
				of_samples = of_samples,
				transform = transform,
				...
			),
			tolower(method) == tolower("UMAP") ~ 	get_reduced_dimensions_UMAP_bulk(
			  .,
       .abundance = !!.abundance,
       .dims = .dims,
       .element = !!.element,
       .feature = !!.feature,
       top = top,
       of_samples = of_samples,
       transform = transform,
       scale = scale,
       ...
			),
			TRUE ~ 	stop("tidybulk says: method must be either \"MDS\" or \"PCA\" or \"tSNE\"")
		)



	if (action == "add"){

		.data |>	dplyr::left_join(.data_processed,	by = quo_name(.element)) |>

			# Attach attributes
			reattach_internals(.data_processed)

	}
	else if (action == "get"){

		.data |>

			# Selecting the right columns
			pivot_sample(!!.element) |>

			dplyr::left_join(.data_processed,	by = quo_name(.element)) |>

			# Attach attributes
			reattach_internals(.data_processed)

	}

	else if (action == "only") .data_processed
	else
		stop(
			"tidybulk says: action must be either \"add\" for adding this ",
			"information to your data frame or \"get\" to just get the information"
		)


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


#' Rotate two dimensions (e.g., principal components) of an arbitrary angle
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description rotate_dimensions() takes as input a `tbl` formatted as 
#'   | <DIMENSION 1> | <DIMENSION 2> | <...> | and calculates the rotated 
#'   dimensional space of the transcript abundance.
#'
#' @importFrom rlang enquo
#'
#' @name rotate_dimensions
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and 
#'   transcript abundance) or `SummarizedExperiment` (more convenient if 
#'   abstracted to tibble with library(tidySummarizedExperiment))
#' @param .element The name of the element column (normally samples).
#'
#' @param dimension_1_column A character string. The column of the dimension 1
#' @param dimension_2_column  A character string. The column of the dimension 2
#' @param rotation_degrees A real number between 0 and 360
#' @param of_samples A boolean. In case the input is a tidybulk object, 
#'   it indicates Whether the element column will be sample or transcript column
#' @param dimension_1_column_rotated A character string. The column of the 
#'   rotated dimension 1 (optional)
#' @param dimension_2_column_rotated A character string. The column of the 
#'   rotated dimension 2 (optional)
#' @param action A character string. Whether to join the new information to the 
#'   input tbl (add), or just get the non-redundant tbl with the new information (get).
#'
#' @details This function to rotate two dimensions such as the reduced dimensions.
#'
#' Underlying custom method:
#' 	rotation = function(m, d) {
#' 		// r = the angle
#' 		// m data matrix
#'    r = d * pi / 180
#'    ((bind_rows(
#' 	  c(`1` = cos(r), `2` = -sin(r)),
#' 	  c(`1` = sin(r), `2` = cos(r))
#'   ) |> as_matrix()) %*% m)
#'  }
#'
#' @return A tbl object with additional columns for the reduced dimensions. 
#'  additional columns for the rotated dimensions. The rotated dimensions will 
#'  be added to the original data set as `<NAME OF DIMENSION> rotated <ANGLE>` 
#'  by default, or as specified in the input arguments.
#'
#' @examples
#' counts.MDS =
#'  tidybulk::se_mini |>
#'  identify_abundant() |>
#'  reduce_dimensions( method="MDS", .dims = 3)
#'
#' counts.MDS.rotated =  rotate_dimensions(counts.MDS, `Dim1`, `Dim2`, 
#'   rotation_degrees = 45, .element = sample)
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
.rotate_dimensions <- function(.data,
																dimension_1_column,
																dimension_2_column,
																rotation_degrees,
																.element = NULL,
																of_samples = TRUE,
																dimension_1_column_rotated = NULL,
																dimension_2_column_rotated = NULL,
																action =	"add") {

  # Fix NOTEs
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
		.data |>
			dplyr::left_join(	.data_processed,	by = quo_name(.element)	) |>

			# Attach attributes
			reattach_internals(.data)
	}
	else if (action == "get"){

		.data |>

			# Selecting the right columns
			select(
				!!.element,
				get_specific_annotation_columns(.data, !!.element)
			) |>
			distinct() |>

			dplyr::left_join(	.data_processed,	by = quo_name(.element)	) |>

			# Attach attributes
			reattach_internals(.data)

	}
	else if (action == "only") .data_processed
	else
		stop(
			"tidybulk says: action must be either \"add\" for adding this ",
			"information to your data frame or \"get\" to just get the information"
		)
}

#' rotate_dimensions
#' @inheritParams rotate_dimensions
#'
#' @docType methods
#' @rdname rotate_dimensions-methods
#'
#' @return A tbl object with additional columns for the reduced dimensions. 
#'   additional columns for the rotated dimensions. The rotated dimensions will 
#'   be added to the original data set as `<NAME OF DIMENSION> rotated <ANGLE>` 
#'   by default, or as specified in the input arguments.
setMethod("rotate_dimensions", "spec_tbl_df", .rotate_dimensions)

#' rotate_dimensions
#' @inheritParams rotate_dimensions
#'
#' @docType methods
#' @rdname rotate_dimensions-methods
#'
#' @return A tbl object with additional columns for the reduced dimensions. 
#'   additional columns for the rotated dimensions. The rotated dimensions will 
#'   be added to the original data set as `<NAME OF DIMENSION> rotated <ANGLE>` 
#'   by default, or as specified in the input arguments.
setMethod("rotate_dimensions", "tbl_df", .rotate_dimensions)

#' rotate_dimensions
#' @inheritParams rotate_dimensions
#'
#' @docType methods
#' @rdname rotate_dimensions-methods
#'
#' @return A tbl object with additional columns for the reduced dimensions. 
#'   additional columns for the rotated dimensions. The rotated dimensions 
#'   will be added to the original data set as 
#'   `<NAME OF DIMENSION> rotated <ANGLE>` by default, or as 
#'   specified in the input arguments.
setMethod("rotate_dimensions", "tidybulk", .rotate_dimensions)




#' Drop redundant elements (e.g., samples) for which feature (e.g., 
#' transcript/gene) abundances are correlated
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description remove_redundancy() takes as input A `tbl` (with at least 
#'   three columns for sample, feature and transcript abundance) or 
#'   `SummarizedExperiment` (more convenient if abstracted to tibble with 
#'   library(tidySummarizedExperiment)) for correlation method or | 
#'   <DIMENSION 1> | <DIMENSION 2> | <...> | for reduced_dimensions method, 
#'   and returns a consistent object (to the input) with dropped 
#'   elements (e.g., samples).
#'
#' @importFrom rlang enquo
#'
#' @name remove_redundancy
#'
#' @param .data A `tbl` (with at least three columns for sample, feature 
#'   and transcript abundance) or `SummarizedExperiment` (more convenient if 
#'   abstracted to tibble with library(tidySummarizedExperiment))
#' @param .element The name of the element column (normally samples).
#' @param .feature The name of the feature column (normally transcripts/genes)
#' @param .abundance The name of the column including the numerical value the 
#'   clustering is based on (normally transcript abundance)
#'
#' @param method A character string. The method to use, correlation and 
#'   reduced_dimensions are available. The latter eliminates one of the most 
#'   proximar pairs of samples in PCA reduced dimensions.
#' @param of_samples A boolean. In case the input is a tidybulk object, 
#'   it indicates Whether the element column will be sample or transcript column
#' @param transform A function that will tranform the counts, by default it 
#'   is log1p for RNA sequencing data, but for avoinding tranformation you 
#'   can use identity
#' @param correlation_threshold A real number between 0 and 1. For correlation based calculation.
#' @param top An integer. How many top genes to select for correlation based method
#' @param Dim_a_column A character string. For reduced_dimension based 
#'   calculation. The column of one principal component
#' @param Dim_b_column A character string. For reduced_dimension based 
#'   calculation. The column of another principal component
#'
#' @param log_transform DEPRECATED - A boolean, whether the value 
#'   should be log-transformed (e.g., TRUE for RNA sequencing data)
#'
#' @details This function removes redundant elements from the original 
#'   data set (e.g., samples or transcripts).
#' For example, if we want to define cell-type specific signatures with 
#' low sample redundancy.
#' This function returns a tibble with dropped redundant elements (e.g., samples).
#' Two redundancy estimation approaches are supported:
#' (i) removal of highly correlated clusters of elements (keeping a 
#' representative) with method="correlation";
#' (ii) removal of most proximal element pairs in a reduced dimensional space.
#'
#' Underlying method for correlation:
#' widyr::pairwise_cor(sample, transcript,count, sort = TRUE, 
#'                     diag = FALSE, upper = FALSE)
#'
#' Underlying custom method for reduced dimensions:
#' select_closest_pairs = function(df) {
#' 		couples <- df |> head(n = 0)
#'
#' 		while (df |> nrow() > 0) {
#' 			pair <- df |>
#' 			arrange(dist) |>
#' 			head(n = 1)
#' 			couples <- couples |> bind_rows(pair)
#' 			df <- df |>
#' 				filter(
#' 					!`sample 1` %in% (pair |> select(1:2) |> as.character()) &
#' 						!`sample 2` %in% (pair |> select(1:2) |> as.character())
#' 				)
#' 		}
#' 		couples
#' 	}
#'
#'
#' @return A tbl object with with dropped redundant elements (e.g., samples).
#'
#' @examples
#'  tidybulk::se_mini |>
#'  identify_abundant() |>
#'    remove_redundancy(
#' 	   .element = sample,
#' 	   .feature = transcript,
#' 	   	.abundance =  count,
#' 	   	method = "correlation"
#' 	   	)
#'
#' counts.MDS =
#'  tidybulk::se_mini |>
#'  identify_abundant() |>
#'   reduce_dimensions( method="MDS", .dims = 3)
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
																				 transform = identity,
																				 Dim_a_column,
																				 Dim_b_column,

																				 # DEPRECATED
																				 log_transform = NULL
																				 )
					 standardGeneric("remove_redundancy"))

# Set internal
.remove_redundancy <- function(.data,
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

  # DEPRECATION OF log_transform
  if (is_present(log_transform) & !is.null(log_transform)) {

    # Signal the deprecation to the user
    deprecate_warn("1.7.4", 
                   "tidybulk::test_differential_abundance(log_transform = )", 
                   details = "The argument log_transform is now deprecated, please use transform.")

    if(log_transform == TRUE) transform = log1p
  }

	# Make col names
	.abundance = enquo(.abundance)
	.element = enquo(.element)
	.feature = enquo(.feature)

	Dim_a_column = enquo(Dim_a_column)
	Dim_b_column = enquo(Dim_b_column)


	if (method == "correlation") {
		# Validate data frame
		if(do_validate()) {
		validation(.data, !!.element, !!.feature, !!.abundance)
		# warning_if_data_is_not_rectangular(.data, !!.element, !!.feature, !!.abundance)
		}

		remove_redundancy_elements_through_correlation(
			.data,
			.abundance = !!.abundance,
			.element = !!.element,
			.feature = !!.feature,
			correlation_threshold = correlation_threshold,
			top = top,
			of_samples = of_samples,
			transform = transform
		)
	}
	else if (method == "reduced_dimensions") {
		# Validate data frame
		# MISSING because feature not needed. I should build a custom function.

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
			"tidybulk says: method must be either \"correlation\" for dropping ",
			"correlated elements or \"reduced_dimension\" to drop the closest pair ",
			"according to two dimensions (e.g., PCA)"
		)

}

#' remove_redundancy
#' @inheritParams remove_redundancy
#'
#' @docType methods
#' @rdname remove_redundancy-methods
#'
#' @return A tbl object with with dropped redundant elements (e.g., samples).
setMethod("remove_redundancy", "spec_tbl_df", .remove_redundancy)

#' remove_redundancy
#' @inheritParams remove_redundancy
#'
#' @docType methods
#' @rdname remove_redundancy-methods
#'
#' @return A tbl object with with dropped redundant elements (e.g., samples).
setMethod("remove_redundancy", "tbl_df", .remove_redundancy)

#' remove_redundancy
#' @inheritParams remove_redundancy
#'
#' @docType methods
#' @rdname remove_redundancy-methods
#'
#' @return A tbl object with with dropped redundant elements (e.g., samples).
setMethod("remove_redundancy", "tidybulk", .remove_redundancy)



#' Adjust transcript abundance for unwanted variation
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description adjust_abundance() takes as input A `tbl` (with at least three 
#' columns for sample, feature and transcript abundance) or 
#' `SummarizedExperiment` (more convenient if abstracted to tibble with 
#' library(tidySummarizedExperiment)) and returns a consistent object 
#' (to the input) with an additional adjusted abundance column. 
#' This method uses scaled counts if present.
#'
#' @importFrom rlang enquo
#'
#' @name adjust_abundance
#'
#' @param .data A `tbl` (with at least three columns for sample, feature 
#' and transcript abundance) or `SummarizedExperiment` (more convenient if 
#' abstracted to tibble with library(tidySummarizedExperiment))
#' @param .factor_unwanted A tidy select, e.g. column names without double 
#' quotation. c(batch, country) These are the factor that we want to adjust for, 
#' including unwanted batcheffect, and unwanted biological effects.
#' @param .factor_of_interest A tidy select, e.g. column names without double 
#' quotation. c(treatment) These are the factor that we want to preserve.
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param method A character string. Methods include combat_seq (default), 
#' combat and limma_remove_batch_effect.
#'
#' @param action A character string. Whether to join the new information to the 
#' input tbl (add), or just get the non-redundant tbl with the new information (get).
#' @param ... Further parameters passed to the function sva::ComBat
#'
#' @param .formula DEPRECATED - A formula with no response variable, 
#' representing the desired linear model where the first covariate is the 
#' factor of interest and the second covariate is the unwanted variation 
#' (of the kind ~ factor_of_interest + batch)
#' @param transform DEPRECATED - A function that will tranform the counts, 
#' by default it is log1p for RNA sequencing data, but for avoinding 
#' tranformation you can use identity
#' @param inverse_transform DEPRECATED - A function that is the inverse of 
#' transform (e.g. expm1 is inverse of log1p). This is needed to tranform 
#' back the counts after analysis.
#' @param log_transform DEPRECATED - A boolean, whether the value should be 
#' log-transformed (e.g., TRUE for RNA sequencing data)
#'
#' @details This function adjusts the abundance for (known) unwanted variation.
#' At the moment just an unwanted covariate is allowed at a time using 
#' Combat (DOI: 10.1093/bioinformatics/bts034)
#'
#' Underlying method:
#' 	sva::ComBat(data, batch = my_batch,	mod = design,	prior.plots = FALSE, ...)
#'
#' @return A consistent object (to the input) with additional columns for 
#' the adjusted counts as `<COUNT COLUMN>_adjusted`
#'
#' @examples
#'
#' cm = tidybulk::se_mini
#' cm$batch = 0
#' cm$batch[colnames(cm) %in% c("SRR1740035", "SRR1740043")] = 1
#'
#' cm |>
#' identify_abundant() |>
#'	adjust_abundance(	.factor_unwanted = batch, 
#'	.factor_of_interest = condition, method="combat"	)
#'
#'
#' @docType methods
#' @rdname adjust_abundance-methods
#' @export
#'
#'
setGeneric("adjust_abundance", function(.data,
                                        # DEPRECATED
                                        .formula = NULL,
                                        .factor_unwanted =NULL,
                                        .factor_of_interest = NULL,
																				.sample = NULL,
																				.transcript = NULL,
																				.abundance = NULL,
																				method = "combat_seq",

																				action = "add",
																				...,

																				# DEPRECATED
																				log_transform = NULL,
																				transform = NULL,
																				inverse_transform = NULL
																				)
	standardGeneric("adjust_abundance"))

# Set internal
.adjust_abundance <- function(.data,
                              # DEPRECATED
                              .formula = NULL,
                              .factor_unwanted = NULL,
                              .factor_of_interest = NULL,
                              .sample = NULL,
                              .transcript = NULL,
                              .abundance = NULL,
                              method = "combat_seq",
                              action = "add",
                              ...,

                              # DEPRECATED
                              log_transform = NULL,
                              transform = NULL,
                              inverse_transform = NULL
                              ) {

  # Fix NOTEs
  . = NULL

  # Get column names
  .sample = enquo(.sample)
  .transcript = enquo(.transcript)
  col_names = get_sample_transcript(.data, .sample, .transcript)
  .sample = col_names$.sample
  .transcript = col_names$.transcript

  # DEPRECATION OF log_transform
  if (is_present(log_transform) & !is.null(log_transform)) {

    # Signal the deprecation to the user
    deprecate_warn("1.7.4", 
                   "tidybulk::test_differential_abundance(log_transform = )", 
                   details = "The argument log_transform is now deprecated, please use transform.")

    if(log_transform){
      transform = log1p
      inverse_transform = expm1
    }
  }

  # DEPRECATION OF log_transform
  if (
    (is_present(transform) & !is.null(transform)) |
    is_present(inverse_transform) & !is.null(inverse_transform)
    ) {

    # Signal the deprecation to the user
    deprecate_warn("1.11.6", "tidybulk::test_differential_abundance(transform = )", 
                   details = "The argument transform and inverse_transform is now deprecated, please use method argument instead specifying \"combat\", \"combat_seq\" or \"limma_remove_batch_effect\".")

  }

  # DEPRECATION OF .formula
  if (is_present(.formula) & !is.null(.formula)) {

    # Signal the deprecation to the user
    deprecate_warn("1.11.6", "tidybulk::test_differential_abundance(.formula = )", 
                   details = "The argument .formula is now deprecated, please use factor_unwanted and factor_of_interest. Using the formula, the first factor is of interest and the second is unwanted")

    # Check that .formula includes at least two covariates
    if (parse_formula(.formula) %>% length %>% st(2))
      stop(
        "The .formula must contain two covariates, the first being the factor of interest, the second being the factor of unwanted variation"
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



	# Get scaled abundance if present, otherwise get abundance (if present get scaled one)
	.abundance =
		enquo(.abundance) %>%
		when(!quo_is_symbol(.) ~
		       get_abundance_norm_if_exists(.data, .)$.abundance, ~ (.))

	# Validate data frame
	if(do_validate()) {
	validation(.data, !!.sample, !!.transcript, !!.abundance)
	warning_if_data_is_not_rectangular(.data, !!.sample, !!.transcript, !!.abundance)
	}

	.data_processed =

		.data %>%

		# Filter abundant if performed
		when(
			".abundant" %in% colnames(.) ~ filter(., .abundant),
			~ {
				warning("tidybulk says: highly abundant transcripts were not ",
				"identified (i.e. identify_abundant()) or filtered (i.e., keep_abundant), ",
				"therefore this operation will be performed on unfiltered data. ",
				"In rare occasions this could be wanted. In standard whole-transcriptome ",
				"workflows is generally unwanted.")
				(.) }
		) |>

		get_adjusted_counts_for_unwanted_variation_bulk(
		  .factor_unwanted = !!.factor_unwanted,
		  .factor_of_interest = !!.factor_of_interest,
			.sample = !!.sample,
			.transcript = !!.transcript,
			.abundance = !!.abundance,
			method = method,
			...
		)

	if (action == "add"){

		.data |>

			# Add adjusted column
			dplyr::left_join(.data_processed,	by = c(quo_name(.transcript), 
			                                         quo_name(.sample))) |>

			# Attach attributes
			reattach_internals(.data)

	}
	else if (action == "get"){

		.data |>

			# Selecting the right columns
			pivot_sample(!!.sample) |>
			#
			# select(
			# 	!!.sample,
			# 	get_x_y_annotation_columns(.data, !!.sample,!!.transcript, !!.abundance, NULL)$horizontal_cols
			# ) |>
			# distinct() |>

			# Add adjusted column
			dplyr::left_join(.data_processed,	by = quo_name(.sample)) |>

			# Attach attributes
			reattach_internals(.data)

	}
	else if (action == "only") .data_processed
	else
		stop(
			"tidybulk says: action must be either \"add\" for adding this information ",
			"to your data frame or \"get\" to just get the information"
		)
}

#' adjust_abundance
#' @inheritParams adjust_abundance
#'
#' @docType methods
#' @rdname adjust_abundance-methods
#'
#' @return A consistent object (to the input) with additional columns for the 
#' adjusted counts as `<COUNT COLUMN>_adjusted`
setMethod("adjust_abundance", "spec_tbl_df", .adjust_abundance)

#' adjust_abundance
#' @inheritParams adjust_abundance
#'
#' @docType methods
#' @rdname adjust_abundance-methods
#'
#' @return A consistent object (to the input) with additional columns for the 
#' adjusted counts as `<COUNT COLUMN>_adjusted`
setMethod("adjust_abundance", "tbl_df", .adjust_abundance)

#' adjust_abundance
#' @inheritParams adjust_abundance
#'
#' @docType methods
#' @rdname adjust_abundance-methods
#'
#' @return A consistent object (to the input) with additional columns for the 
#' adjusted counts as `<COUNT COLUMN>_adjusted`
setMethod("adjust_abundance", "tidybulk", .adjust_abundance)



#' Aggregates multiple counts from the same samples (e.g., from isoforms), 
#' concatenates other character columns, and averages other numeric columns
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description aggregate_duplicates() takes as input A `tbl` (with at least 
#' three columns for sample, feature and transcript abundance) or 
#' `SummarizedExperiment` (more convenient if abstracted to tibble with 
#' library(tidySummarizedExperiment)) and returns a consistent object 
#' (to the input) with aggregated transcripts that were duplicated.
#'
#' @importFrom rlang enquo
#'
#' @name aggregate_duplicates
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and 
#' transcript abundance) or `SummarizedExperiment` (more convenient if 
#' abstracted to tibble with library(tidySummarizedExperiment))
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#'
#' @param aggregation_function A function for counts aggregation (e.g., sum,  
#' median, or mean)
#' @param keep_integer A boolean. Whether to force the aggregated counts to integer
#'
#' @details This function aggregates duplicated transcripts (e.g., isoforms, ensembl).
#' For example, we often have to convert ensembl symbols to gene/transcript symbol,
#'  but in doing so we have to deal with duplicates. `aggregate_duplicates` takes a tibble
#'  and column names (as symbols; for `sample`, `transcript` and `count`) as arguments and
#'  returns a tibble with aggregate transcript with the same name. All the rest of the column
#'  are appended, and factors and boolean are appended as characters.
#'
#'  Underlying custom method:
#'  data |>
#' 		filter(n_aggr > 1) |>
#' 		group_by(!!.sample,!!.transcript) |>
#' 		dplyr::mutate(!!.abundance := !!.abundance |> aggregation_function())
#'
#' @return A consistent object (to the input) with aggregated transcript 
#' abundance and annotation
#'
#'
#'
#'
#' @examples
#'
#' # Create a aggregation column
#' se_mini = tidybulk::se_mini
#' SummarizedExperiment::rowData(se_mini )$gene_name = rownames(se_mini )
#'
#'    aggregate_duplicates(
#'      se_mini,
#'    .transcript = gene_name
#'    )
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
.aggregate_duplicates <- function(.data,
																	.sample = NULL,
																	.transcript = NULL,
																	.abundance = NULL,
																	aggregation_function = sum,
																	keep_integer = TRUE)  {

  # Fix NOTEs
  . = NULL

	# Make col names
  # Get column names
  .sample = enquo(.sample)
  .transcript = enquo(.transcript)
  .abundance = enquo(.abundance)
  col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
  .sample = col_names$.sample
  .transcript = col_names$.transcript
  .abundance = col_names$.abundance

	# Validate data frame
	if(do_validate()) validation(.data,
						 !!.sample,
						 !!.transcript,
						 !!.abundance,
						 skip_dupli_check = TRUE)

	# If I have a small normal data set
	if(.data |> pull(!!.sample) |> unique() |> length() |> st(100)){
  	aggregate_duplicated_transcripts_bulk(
  		.data,

  		.sample = !!.sample,
  		.transcript = !!.transcript,
  		.abundance = !!.abundance,
  		aggregation_function = aggregation_function,
  		keep_integer = TRUE
  	)
}
	# If I have a big data set
  else {

    message("tidybulk says: for big data sets (>100 samples) this efficient ",
    "implementation aggregates count columns and keeps the first instance for sample and transcript annotations")

    aggregate_duplicated_transcripts_DT(
      .data,
      .sample = !!.sample,
      .transcript = !!.transcript,
      .abundance = !!.abundance,
      aggregation_function = aggregation_function,
      keep_integer = TRUE
    )
  }

}

#' aggregate_duplicates
#' @inheritParams aggregate_duplicates
#'
#' @docType methods
#' @rdname aggregate_duplicates-methods
#'
#' @return A consistent object (to the input) with aggregated transcript 
#' abundance and annotation
setMethod("aggregate_duplicates", "spec_tbl_df", .aggregate_duplicates)

#' aggregate_duplicates
#' @inheritParams aggregate_duplicates
#'
#' @docType methods
#' @rdname aggregate_duplicates-methods
#'
#' @return A consistent object (to the input) with aggregated transcript 
#' abundance and annotation
setMethod("aggregate_duplicates", "tbl_df", .aggregate_duplicates)

#' aggregate_duplicates
#' @inheritParams aggregate_duplicates
#'
#' @docType methods
#' @rdname aggregate_duplicates-methods
#'
#' @return A consistent object (to the input) with aggregated transcript 
#' abundance and annotation
setMethod("aggregate_duplicates", "tidybulk", .aggregate_duplicates)



#' Get cell type proportions from samples
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description deconvolve_cellularity() takes as input A `tbl` (with at least 
#' three columns for sample, feature and transcript abundance) or 
#' `SummarizedExperiment` (more convenient if abstracted to tibble with 
#' library(tidySummarizedExperiment)) and returns a consistent object (to the 
#' input) with the estimated cell type abundance for each sample
#'
#' @importFrom rlang enquo
#'
#' @name deconvolve_cellularity
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and 
#' transcript abundance) or `SummarizedExperiment` (more convenient if 
#' abstracted to tibble with library(tidySummarizedExperiment))
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param reference A data frame. The methods cibersort and llsr can accept a 
#' custom rectangular dataframe with genes as rows names, cell types as column 
#' names and gene-transcript abundance as values. For exampler tidybulk::X_cibersort. 
#' The transcript/cell_type data frame of integer transcript abundance. If NULL, 
#' the default reference for each algorithm will be used. For llsr will be LM22.
#' @param method A character string. The method to be used. At the moment 
#' Cibersort (default, can accept custom reference), epic (can accept custom 
#' reference) and llsr (linear least squares regression, can accept custom 
#' reference), mcp_counter, quantiseq, xcell are available.
#' @param prefix A character string. The prefix you would like to add to the 
#' result columns. It is useful if you want to reshape data.
#' @param action A character string. Whether to join the new information to the 
#' input tbl (add), or just get the non-redundant tbl with the new information (get).
#' @param ... Further parameters passed to the function Cibersort
#'
#' @details This function infers the cell type composition of our samples
#' (with the algorithm Cibersort; Newman et al., 10.1038/nmeth.3337).
#'
#' Underlying method:
#' CIBERSORT(Y = data, X = reference, ...)
#'
#' @return A consistent object (to the input) including additional columns 
#' for each cell type estimated
#'
#' @examples
#'
#' # Subsetting for time efficiency
#' tidybulk::se_mini |> deconvolve_cellularity(cores = 1)
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
																							reference = NULL,
																							method = "cibersort",
																							prefix = "",
																							action = "add",
																							...)
	standardGeneric("deconvolve_cellularity"))

# Set internal
.deconvolve_cellularity <- function(.data,
																		 .sample = NULL,
																		 .transcript = NULL,
																		 .abundance = NULL,
																		 reference = NULL,
																		 method = "cibersort",
																		 prefix = "",
																		 action = "add",
																		 ...)  {

  # Fix NOTEs
  . = NULL

	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
	.sample = col_names$.sample
	.transcript = col_names$.transcript
	.abundance = col_names$.abundance

	# # Check that reference is matrix
	# if(!is.null(reference) & reference %>% class %>% equals("data.frame") %>% not())
	# 	stop("tidybulk says: reference must be NULL or a data.frame")

	# Validate data frame
	if(do_validate()) {
	validation(.data, !!.sample, !!.transcript, !!.abundance)
	warning_if_data_is_not_rectangular(.data, !!.sample, !!.transcript, !!.abundance)
	}

	.data_processed =
		get_cell_type_proportions(
		.data,
		.sample = !!.sample,
		.transcript = !!.transcript,
		.abundance = !!.abundance,
		reference = reference,
		method = method,
		prefix = prefix,
		...
	)

	if (action == "add"){
		.data |>
			# Add new annotation
			dplyr::left_join(.data_processed, by = quo_name(.sample)) |>

			# Attach attributes
			reattach_internals(.data_processed)
	}

	else if (action == "get"){
		.data |>
			# Selecting the right columns
			pivot_sample(!!.sample) |>
			#
			# select(
			# 	!!.sample,
			# 	get_x_y_annotation_columns(.data, !!.sample,!!.transcript, !!.abundance, NULL)$horizontal_cols
			# ) |>
			# distinct() |>

			# Add new annotation
			dplyr::left_join(.data_processed, by = quo_name(.sample)) |>

			# Attach attributes
			reattach_internals(.data_processed)
	}

	else if (action == "only") .data_processed
	else
		stop(
			"tidybulk says: action must be either \"add\" for adding this ",
			"information to your data frame or \"get\" to just get the information"
		)
}

#' deconvolve_cellularity
#' @inheritParams deconvolve_cellularity
#'
#' @docType methods
#' @rdname deconvolve_cellularity-methods
#'
#' @return A consistent object (to the input) including additional columns 
#' for each cell type estimated
setMethod("deconvolve_cellularity",
					"spec_tbl_df",
					.deconvolve_cellularity)

#' deconvolve_cellularity
#' @inheritParams deconvolve_cellularity
#'
#' @docType methods
#' @rdname deconvolve_cellularity-methods
#'
#' @return A consistent object (to the input) including additional columns 
#' for each cell type estimated
setMethod("deconvolve_cellularity", "tbl_df", .deconvolve_cellularity)

#' deconvolve_cellularity
#' @inheritParams deconvolve_cellularity
#'
#' @docType methods
#' @rdname deconvolve_cellularity-methods
#'
#' @return A consistent object (to the input) including additional columns 
#' for each cell type estimated
setMethod("deconvolve_cellularity",
					"tidybulk",
					.deconvolve_cellularity)


#' Get ENTREZ id from gene SYMBOL
#'
#' @param .data A tt or tbl object.
#' @param .transcript A character. The name of the gene symbol column.
#' @param .sample The name of the sample column
#'
#' @return A tbl
#'
#' @examples
#'
#' # This function was designed for data.frame
#' # Convert from SummarizedExperiment for this example. It is NOT reccomended.
#'
#' tidybulk::se_mini |> tidybulk() |> as_tibble() |> 
#' symbol_to_entrez(.transcript = .feature, .sample = .sample)
#'
#' @export
#'
symbol_to_entrez = function(.data,
														.transcript = NULL,
														.sample = NULL) {

  # Fix NOTEs
  . = NULL

	# Get column names
	.transcript = enquo(.transcript)
	.sample = enquo(.sample)
	col_names = get_sample_transcript(.data, .sample, .transcript)
	.transcript = col_names$.transcript

	# Check if package is installed, otherwise install
	if (find.package("org.Hs.eg.db", quiet = TRUE) |> length() |> equals(0)) {
		message("Installing org.Hs.eg.db needed for annotation")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("org.Hs.eg.db", ask = FALSE)
	}

	.data |>

		# Solve the lower case
		mutate(transcript_upper := (!!.transcript) |> toupper()) %>%

		# Join
		dplyr::left_join(
			# Get entrez mapping 1:1
			AnnotationDbi::mapIds(
				org.Hs.eg.db::org.Hs.eg.db,
				(.) |> pull(transcript_upper) |> as.character() |> unique(),
				'ENTREZID',
				'SYMBOL'
			) |>
				enframe(name = "transcript_upper", value = "entrez") |>
				filter(entrez |> is.na() |> not()) |>
				group_by(transcript_upper) |>
				slice(1) |>
				ungroup(),
			by = "transcript_upper"
		) |>

		# Eliminate the upper case
		select(-transcript_upper)

}


#' Get DESCRIPTION from gene SYMBOL for Human and Mouse
#'
#' @param .data A tt or tbl object.
#' @param .transcript A character. The name of the gene symbol column.
#'
#' @return A tbl
#'
#' @examples
#'
#' describe_transcript(tidybulk::se_mini)
#'
#' @docType methods
#' @rdname describe_transcript-methods
#' @export
#'
#'
setGeneric("describe_transcript", function(.data,
																					 .transcript = NULL)
	standardGeneric("describe_transcript"))

#'
.describe_transcript = function(.data,
														.transcript = NULL) {

  # Fix NOTEs
  . = NULL

	# Get column names
	.transcript = enquo(.transcript)
	col_names = get_transcript(.data, .transcript)
	.transcript = col_names$.transcript


	# Check if package is installed, otherwise install
	if (find.package("org.Hs.eg.db", quiet = TRUE) |> length() |> equals(0)) {
		message("Installing org.Hs.eg.db needed for differential transcript abundance analyses")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("org.Hs.eg.db", ask = FALSE)
	}

	# Check if package is installed, otherwise install
	if (find.package("org.Mm.eg.db", quiet = TRUE) |> length() |> equals(0)) {
		message("Installing org.Mm.eg.db needed for differential transcript abundance analyses")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("org.Mm.eg.db", ask = FALSE)
	}

	# Check if package is installed, otherwise install
	if (find.package("AnnotationDbi", quiet = TRUE) |> length() |> equals(0)) {
		message("Installing AnnotationDbi needed for differential transcript abundance analyses")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("AnnotationDbi", ask = FALSE)
	}

	description_df =


		# Human
		tryCatch(suppressMessages(AnnotationDbi::mapIds(
			org.Hs.eg.db::org.Hs.eg.db,
			keys = pull(.data, !!.transcript) |> unique() |> as.character(),  #ensembl_symbol_mapping$transcript %>% unique,
			column = "GENENAME",
			keytype = "SYMBOL",
			multiVals = "first"
		))  %>%
			.[!is.na(.)], error = function(x){}) %>%

		# Mouse
		c(
			tryCatch(suppressMessages(AnnotationDbi::mapIds(
				org.Mm.eg.db::org.Mm.eg.db,
				keys = pull(.data, !!.transcript) |> unique() |> as.character(),  #ensembl_symbol_mapping$transcript %>% unique,
				column = "GENENAME",
				keytype = "SYMBOL",
				multiVals = "first"
			)) %>% .[!is.na(.)], error = function(x){})

		) |>

		# Parse
		unlist() |>
		#unique() |>
		enframe(name = quo_name(.transcript), value = "description") |>

		# Select just one per transcript
		distinct() |>
		group_by(!!.transcript) |>
		slice(1) |>
		ungroup()

	.data |>
		left_join(description_df, by = quo_name(.transcript))
}


#' describe_transcript
#' @inheritParams describe_transcript
#'
#' @docType methods
#' @rdname describe_transcript-methods
#'
#' @return A consistent object (to the input) including additional columns for transcript symbol
setMethod("describe_transcript", "spec_tbl_df", .describe_transcript)

#' describe_transcript
#' @inheritParams describe_transcript
#'
#' @docType methods
#' @rdname describe_transcript-methods
#'
#' @return A consistent object (to the input) including additional columns for transcript symbol
setMethod("describe_transcript", "tbl_df", .describe_transcript)

#' describe_transcript
#' @inheritParams describe_transcript
#'
#' @docType methods
#' @rdname describe_transcript-methods
#'
#' @return A consistent object (to the input) including additional columns for transcript symbol
setMethod("describe_transcript", "tidybulk", .describe_transcript)


#' Add transcript symbol column from ensembl id for human and mouse data
#'
#' \lifecycle{questioning}
#'
#' @description ensembl_to_symbol() takes as input a `tbl` (with at least 
#' three columns for sample, feature and transcript abundance) or 
#' `SummarizedExperiment` (more convenient if abstracted to tibble with 
#' library(tidySummarizedExperiment)) and returns a consistent object 
#' (to the input) with the additional transcript symbol column
#'
#' @importFrom rlang enquo
#'
#' @name ensembl_to_symbol
#'
#' @param .data a `tbl` (with at least three columns for sample, feature 
#' and transcript abundance) or `SummarizedExperiment` (more convenient if 
#' abstracted to tibble with library(tidySummarizedExperiment))
#' @param .ensembl A character string. The column that is represents 
#' ensembl gene id
#'
#' @param action A character string. Whether to join the new information 
#' to the input tbl (add), or just get the non-redundant tbl with the new 
#' information (get).
#'
#' @details This is useful since different resources use ensembl IDs while 
#' others use gene symbol IDs. At the moment this work for human (genes and 
#' transcripts) and mouse (genes) data.
#'
#' @return A consistent object (to the input) including additional columns 
#' for transcript symbol
#'
#'
#' @examples
#'
#' # This function was designed for data.frame
#' # Convert from SummarizedExperiment for this example. It is NOT reccomended.
#'
#' tidybulk::se_mini |> tidybulk() |> as_tibble() |> ensembl_to_symbol(.feature)
#'
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
.ensembl_to_symbol <- function(.data,
															.ensembl,
															action = "add") {

  # Fix NOTEs
  . = NULL

	# Make col names
	.ensembl = enquo(.ensembl)

	.data_processed = get_symbol_from_ensembl(.data,!!.ensembl)

	if (action == "add"){

		# Add new symbols column
		.data |>
			dplyr::left_join(.data_processed, by=quo_name(.ensembl)) |>

			# Attach attributes
			reattach_internals(.data)

	}
	# else if (action == "get"){
	#
	# 	# Add new symbols column
	# 	.data |>
	#
	#
	# 		dplyr::left_join(.data_processed) |>
	#
	# 		# Attach attributes
	# 		reattach_internals(.data)
	#
	# }

	else if (action == "only") .data_processed

	else
		stop(
			"tidybulk says: action must be either \"add\" for adding this ",
			"information to your data frame or \"get\" to just get the information"
		)

}

#' ensembl_to_symbol
#' @inheritParams ensembl_to_symbol
#'
#' @docType methods
#' @rdname ensembl_to_symbol-methods
#'
#' @return A consistent object (to the input) including additional columns 
#' for transcript symbol
setMethod("ensembl_to_symbol", "spec_tbl_df", .ensembl_to_symbol)

#' ensembl_to_symbol
#' @inheritParams ensembl_to_symbol
#'
#' @docType methods
#' @rdname ensembl_to_symbol-methods
#'
#' @return A consistent object (to the input) including additional columns 
#' for transcript symbol
setMethod("ensembl_to_symbol", "tbl_df", .ensembl_to_symbol)

#' ensembl_to_symbol
#' @inheritParams ensembl_to_symbol
#'
#' @docType methods
#' @rdname ensembl_to_symbol-methods
#'
#' @return A consistent object (to the input) including additional columns 
#' for transcript symbol
setMethod("ensembl_to_symbol", "tidybulk", .ensembl_to_symbol)


#' Perform differential transcription testing using edgeR quasi-likelihood 
#' (QLT), edgeR likelihood-ratio (LR), limma-voom, 
#' limma-voom-with-quality-weights or DESeq2
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description test_differential_abundance() takes as input A `tbl` 
#' (with at least three columns for sample, feature and transcript abundance) 
#' or `SummarizedExperiment` (more convenient if abstracted to tibble with 
#' library(tidySummarizedExperiment)) and returns a consistent object (to the 
#' input) with additional columns for the statistics from the hypothesis test.
#'
#' @importFrom rlang enquo
#'
#' @name test_differential_abundance
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and 
#' transcript abundance) or `SummarizedExperiment` (more convenient if 
#' abstracted to tibble with library(tidySummarizedExperiment))
#' @param .formula A formula representing the desired linear model. If there 
#' is more than one factor, they should be in the order factor of interest + 
#' additional factors.
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#'
#' @param contrasts This parameter takes the format of the contrast parameter 
#' of the method of choice. For edgeR and limma-voom is a character vector. 
#' For DESeq2 is a list including a character vector of length three. The first 
#' covariate is the one the model is tested against (e.g., ~ factor_of_interest)
#' @param method A string character. Either "edgeR_quasi_likelihood" (i.e., QLF), 
#' "edgeR_likelihood_ratio" (i.e., LRT), "edger_robust_likelihood_ratio", 
#' "DESeq2", "limma_voom", "limma_voom_sample_weights"
#' @param test_above_log2_fold_change A positive real value. This works for edgeR 
#' and limma_voom methods. It uses the `treat` function, which tests that the 
#' difference in abundance is bigger than this threshold rather than zero 
#' \url{https://pubmed.ncbi.nlm.nih.gov/19176553}.
#' @param scaling_method A character string. The scaling method passed to the 
#' back-end functions: edgeR and limma-voom (i.e., edgeR::calcNormFactors; 
#' "TMM","TMMwsp","RLE","upperquartile"). Setting the parameter to \"none\" 
#' will skip the compensation for sequencing-depth for the method edgeR or limma-voom.
#' @param omit_contrast_in_colnames If just one contrast is specified you can 
#' choose to omit the contrast label in the colnames.
#' @param prefix A character string. The prefix you would like to add to the 
#' result columns. It is useful if you want to compare several methods.
#' @param action A character string. Whether to join the new information to the 
#' input tbl (add), or just get the non-redundant tbl with the new information (get).
#' @param significance_threshold DEPRECATED - A real between 0 and 1 (usually 0.05).
#' @param fill_missing_values DEPRECATED - A boolean. Whether to fill missing 
#' sample/transcript values with the median of the transcript. This is rarely needed.
#' @param .contrasts DEPRECATED - This parameter takes the format of the 
#' contrast parameter of the method of choice. For edgeR and limma-voom is a 
#' character vector. For DESeq2 is a list including a character vector of length 
#' three. The first covariate is the one the model is tested against (e.g., ~ factor_of_interest)
#' @param ... Further arguments passed to some of the internal functions. 
#' Currently, it is needed just for internal debug.
#'
#'
#' @details This function provides the option to use edgeR 
#' \url{https://doi.org/10.1093/bioinformatics/btp616}, limma-voom 
#' \url{https://doi.org/10.1186/gb-2014-15-2-r29}, limma_voom_sample_weights 
#' \url{https://doi.org/10.1093/nar/gkv412} or  DESeq2 
#' \url{https://doi.org/10.1186/s13059-014-0550-8} to perform the testing.
#' All methods use raw counts, irrespective of if scale_abundance or 
#' adjust_abundance have been calculated, therefore it is essential to add 
#' covariates such as batch effects (if applicable) in the formula.
#'
#' Underlying method for edgeR framework:
#' 
#' 	.data |>
#'
#' 	# Filter
#'	keep_abundant(
#'			factor_of_interest = !!(as.symbol(parse_formula(.formula)[1])),
#'			minimum_counts = minimum_counts,
#'			minimum_proportion = minimum_proportion
#'		) |>
#'
#'			# Format
#'			select(!!.transcript,!!.sample,!!.abundance) |>
#'			spread(!!.sample,!!.abundance) |>
#'			as_matrix(rownames = !!.transcript) %>%
#'
#'			# edgeR
#'			edgeR::DGEList(counts = .) |>
#'			edgeR::calcNormFactors(method = scaling_method) |>
#'			edgeR::estimateDisp(design) |>
#'
#'			# Fit
#'			edgeR::glmQLFit(design) |> // or glmFit according to choice
#'			edgeR::glmQLFTest(coef = 2, contrast = my_contrasts) // or glmLRT according to choice
#'
#'
#'
#'	Underlying method for DESeq2 framework:
#'	
#'	keep_abundant(
#'			factor_of_interest = !!as.symbol(parse_formula(.formula)[[1]]),
#'			minimum_counts = minimum_counts,
#'			minimum_proportion = minimum_proportion
#'	) |>
#'
#'	# DESeq2
#'	DESeq2::DESeqDataSet(design = .formula) |>
#'	DESeq2::DESeq() |>
#'	DESeq2::results()
#'
#'
#'
#' Underlying method for glmmSeq framework:
#'
#' counts =
#' .data %>%
#'   assay(my_assay)
#' 
#' # Create design matrix for dispersion, removing random effects
#' design =
#'   model.matrix(
#'     object = .formula |> eliminate_random_effects(),
#'     data = metadata
#'   )
#' 
#' dispersion = counts |> edgeR::estimateDisp(design = design) %$% tagwise.dispersion |> setNames(rownames(counts))
#' 
#'   glmmSeq( .formula,
#'            countdata = counts ,
#'            metadata =   metadata |> as.data.frame(),
#'            dispersion = dispersion,
#'            progress = TRUE,
#'            method = method |> str_remove("(?i)^glmmSeq_" ),
#'   )
#'   
#' 
#' @return A consistent object (to the input) with additional columns for the 
#' statistics from the test (e.g.,  log fold change, p-value and false discovery rate).
#'
#'
#' @examples
#'
#'  # edgeR
#'
#'  tidybulk::se_mini |>
#'  identify_abundant() |>
#' 	test_differential_abundance( ~ condition )
#'
#' 	# The function `test_differential_abundance` operates with contrasts too
#'
#'  tidybulk::se_mini |>
#'  identify_abundant(factor_of_interest = condition) |>
#'  test_differential_abundance(
#' 	    ~ 0 + condition,
#' 	    contrasts = c( "conditionTRUE - conditionFALSE")
#'  )
#'
#'  # DESeq2 - equivalent for limma-voom
#'
#' my_se_mini = tidybulk::se_mini
#' my_se_mini$condition  = factor(my_se_mini$condition)
#'
#' # demontrating with `fitType` that you can access any arguments to DESeq()
#' my_se_mini  |>
#'    identify_abundant(factor_of_interest = condition) |>
#'        test_differential_abundance( ~ condition, method="deseq2", fitType="local")
#'
#' # testing above a log2 threshold, passes along value to lfcThreshold of results()
#' res <- my_se_mini  |>
#'    identify_abundant(factor_of_interest = condition) |>
#'         test_differential_abundance( ~ condition, method="deseq2",
#'             fitType="local",
#'             test_above_log2_fold_change=4 )
#'
#' # Use random intercept and random effect models
#'
#'  se_mini[1:50,] |>
#'   identify_abundant(factor_of_interest = condition) |>
#'   test_differential_abundance(
#'     ~ condition + (1 + condition | time),
#'     method = "glmmseq_lme4", cores = 1
#'   )
#'
#' # confirm that lfcThreshold was used
#' \dontrun{
#'     res |>
#'         mcols() |>
#'         DESeq2::DESeqResults() |>
#'         DESeq2::plotMA()
#' }
#'
#' # The function `test_differential_abundance` operates with contrasts too
#'
#'  my_se_mini |>
#'  identify_abundant() |>
#'  test_differential_abundance(
#' 	    ~ 0 + condition,
#' 	    contrasts = list(c("condition", "TRUE", "FALSE")),
#' 	    method="deseq2",
#'          fitType="local"
#'  )
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
																									 contrasts = NULL,
																									 method = "edgeR_quasi_likelihood",
																									 test_above_log2_fold_change = NULL,
																									 scaling_method = "TMM",
																									 omit_contrast_in_colnames = FALSE,
																									 prefix = "",
																									 action = "add",
																									 ...,

																									 # DEPRECATED
																									 significance_threshold = NULL,
																									 fill_missing_values = NULL,
																									 .contrasts = NULL
																									)
					 standardGeneric("test_differential_abundance"))

# Set internal
#' @importFrom rlang inform
.test_differential_abundance <- function(.data,
																					.formula,
																					.sample = NULL,
																					.transcript = NULL,
																					.abundance = NULL,
																					contrasts = NULL,
																					method = "edgeR_quasi_likelihood",
																					test_above_log2_fold_change = NULL,
																					scaling_method = "TMM",
																					omit_contrast_in_colnames = FALSE,
																					prefix = "",
																					action = "add",
																					...,
																					# DEPRECATED
																					significance_threshold = NULL,
																					fill_missing_values = NULL,
																					.contrasts = NULL
																				) {

  # Fix NOTEs
  . = NULL

	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
	.sample = col_names$.sample
	.transcript = col_names$.transcript
	.abundance = col_names$.abundance

	# DEPRECATION OF significance_threshold
	if (is_present(significance_threshold) & !is.null(significance_threshold)) {

		# Signal the deprecation to the user
		deprecate_warn("1.1.7", 
		               "tidybulk::test_differential_abundance(significance_threshold = )", 
		               details = "The argument significance_threshold is now deprecated, tigether with the column significance.")

	}

	# DEPRECATION OF fill_missing_values
	if (is_present(fill_missing_values) & !is.null(significance_threshold)) {

		# Signal the deprecation to the user
		deprecate_warn("1.1.7", 
     "tidybulk::test_differential_abundance(fill_missing_values = )", 
     details = "The argument fill_missing_values is now deprecated, you will receive a warning/error instead. Please use externally the methods fill_missing_abundance or impute_missing_abundance instead.")
	}

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
=====================================", .frequency_id = "All testing methods use raw counts",  .frequency = "once")

	# Validate data frame
	if(do_validate()) {
	validation(.data, !!.sample, !!.transcript, !!.abundance)
	warning_if_data_is_not_rectangular(.data, !!.sample, !!.transcript, !!.abundance)
	}

	# Test test_above_log2_fold_change
	if(!is.null(test_above_log2_fold_change) && test_above_log2_fold_change < 0)
		stop("tidybulk says: test_above_log2_fold_change should be a positive real or NULL")

	.data_processed =
		.data %>%

		# Filter abundant if performed
		when(
			".abundant" %in% colnames(.) ~ filter(., .abundant),
			~ {
				warning("tidybulk says: highly abundant transcripts were not identified ",
				"(i.e. identify_abundant()) or filtered (i.e., keep_abundant), therefore ",
				"this operation will be performed on unfiltered data. In rare occasions ",
				"this could be wanted. In standard whole-transcriptome workflows is generally unwanted.")
				(.)
			}
		) %>%

		# Choose method
		when(

			# edgeR
			tolower(method) %in% c("edger_quasi_likelihood", "edger_likelihood_ratio", 
			                       "edger_robust_likelihood_ratio") ~
			get_differential_transcript_abundance_bulk(
				.,
				.formula,
				.sample = !!.sample,
				.transcript = !!.transcript,
				.abundance = !!.abundance,
				.contrasts = contrasts,
				method = method,
				test_above_log2_fold_change = test_above_log2_fold_change,
				scaling_method = scaling_method,
				omit_contrast_in_colnames = omit_contrast_in_colnames,
				prefix = prefix,
				...
			),

			# Voom
			grepl("limma_voom", method) ~ get_differential_transcript_abundance_bulk_voom(
					.,
					.formula,
					.sample = !!.sample,
					.transcript = !!.transcript,
					.abundance = !!.abundance,
					.contrasts = contrasts,
					method = method,
					test_above_log2_fold_change = test_above_log2_fold_change,
					scaling_method = scaling_method,
					omit_contrast_in_colnames = omit_contrast_in_colnames,
					prefix = prefix,
					...
				),

			# DESeq2
			tolower(method)=="deseq2" ~ get_differential_transcript_abundance_deseq2(
				.,
				.formula,
				.sample = !!.sample,
				.transcript = !!.transcript,
				.abundance = !!.abundance,
				.contrasts = contrasts,
				method = method,
				test_above_log2_fold_change = test_above_log2_fold_change,
				scaling_method = scaling_method,
				omit_contrast_in_colnames = omit_contrast_in_colnames,
				prefix = prefix,
				...
			),

			# glmmseq
			tolower(method) %in% c("glmmseq_lme4", "glmmseq_glmmtmb") ~ get_differential_transcript_abundance_glmmSeq(
			  .,
			  .formula,
			  .sample = !!.sample,
			  .transcript = !!.transcript,
			  .abundance = !!.abundance,
			  .contrasts = contrasts,
			  method = method,
			  test_above_log2_fold_change = test_above_log2_fold_change,
			  scaling_method = scaling_method,
			  omit_contrast_in_colnames = omit_contrast_in_colnames,
			  prefix = prefix,
			  ...
			),

			# Else error
			TRUE ~  stop('tidybulk says: the only methods supported at the moment are',
			             '"edgeR_quasi_likelihood" (i.e., QLF), "edgeR_likelihood_ratio"',
			             '(i.e., LRT), "limma_voom", "limma_voom_sample_weights", ',
			             '"DESeq2", "glmmseq_lme4", "glmmseq_glmmTMB"')
		)


	if (action == "add"){

		.data |>
			dplyr::left_join(.data_processed, by = quo_name(.transcript)) |>

			# Attach attributes
			reattach_internals(.data_processed)

	}
	else if (action == "get"){

		.data |>

			# Selecting the right columns
			pivot_transcript(!!.transcript) |>

			dplyr::left_join(.data_processed, by = quo_name(.transcript)) |>

			# Attach attributes
			reattach_internals(.data_processed)

	}
	else if (action == "only") .data_processed
	else
		stop(
			"tidybulk says: action must be either \"add\" for adding this ",
			"information to your data frame or \"get\" to just get the information"
		)
}

#' test_differential_abundance
#' @inheritParams test_differential_abundance
#'
#' @docType methods
#' @rdname test_differential_abundance-methods
#'
#' @return A consistent object (to the input) with additional columns for the 
#' statistics from the test (e.g.,  log fold change, p-value and false discovery rate).
setMethod("test_differential_abundance",
					"spec_tbl_df",
					.test_differential_abundance)

#' test_differential_abundance
#' @inheritParams test_differential_abundance
#'
#' @docType methods
#' @rdname test_differential_abundance-methods
#'
#' @return A consistent object (to the input) with additional columns for the 
#' statistics from the hypothesis test (e.g.,  log fold change, p-value and 
#' false discovery rate).
setMethod("test_differential_abundance",
					"tbl_df",
					.test_differential_abundance)

#' test_differential_abundance
#' @inheritParams test_differential_abundance
#'
#' @docType methods
#' @rdname test_differential_abundance-methods
#'
#' @return A consistent object (to the input) with additional columns for the 
#' statistics from the hypothesis test (e.g.,  log fold change, p-value and 
#' false discovery rate).
setMethod("test_differential_abundance",
					"tidybulk",
					.test_differential_abundance)





#' Keep variable transcripts
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description keep_variable() takes as input A `tbl` (with at least three 
#' columns for sample, feature and transcript abundance) or 
#' `SummarizedExperiment` (more convenient if abstracted to tibble with 
#' library(tidySummarizedExperiment)) and returns a consistent object 
#' (to the input) with additional columns for the statistics from the 
#' hypothesis test.
#'
#' @importFrom rlang enquo
#'
#' @name keep_variable
#'
#' @param .data A `tbl` (with at least three columns for sample, feature 
#' and transcript abundance) or `SummarizedExperiment` (more convenient if 
#' abstracted to tibble with library(tidySummarizedExperiment))
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param top Integer. Number of top transcript to consider
#' @param transform A function that will tranform the counts, by default it is 
#' log1p for RNA sequencing data, but for avoinding tranformation you can use 
#' identity
#'
#' @param log_transform DEPRECATED - A boolean, whether the value should be 
#' log-transformed (e.g., TRUE for RNA sequencing data)
#'
#' @details At the moment this function uses edgeR 
#' \url{https://doi.org/10.1093/bioinformatics/btp616}
#'
#' @return A consistent object (to the input) with additional columns for the 
#' statistics from the hypothesis test (e.g.,  log fold change, p-value and 
#' false discovery rate).
#'
#' Underlying method:
#' 	s <- rowMeans((x - rowMeans(x)) ^ 2)
#'	o <- order(s, decreasing = TRUE)
#'	x <- x[o[1L:top], , drop = FALSE]
#'	variable_trancripts = rownames(x)
#'
#' @examples
#'
#' 	keep_variable(tidybulk::se_mini, top = 500)
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
																			 transform = log1p,

																			 # DEPRECATED
																			 log_transform = TRUE
																			 )
	standardGeneric("keep_variable"))

# Set internal
.keep_variable <- function(.data,
  													.sample = NULL,
  													.transcript = NULL,
  													.abundance = NULL,
  													top = 500,
  													transform = log1p,
  
  													# DEPRECATED
  													log_transform = NULL) {

  # Fix NOTEs
  . = NULL

  # DEPRECATION OF log_transform
  if (is_present(log_transform) & !is.null(log_transform)) {

    # Signal the deprecation to the user
    deprecate_warn("1.7.4", 
                   "tidybulk::test_differential_abundance(log_transform = )", 
                   details = "The argument log_transform is now deprecated, please use transform.")

    if(log_transform == TRUE) transform = log1p
  }

	# Make col names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	# Validate data frame
	if(do_validate()) {
	validation(.data, !!.sample, !!.transcript, !!.abundance)
	warning_if_data_is_not_rectangular(.data, !!.sample, !!.transcript, !!.abundance)
	}

	keep_variable_transcripts(
		.data,
		.sample = !!.sample,
		.transcript = !!.transcript,
		.abundance = !!.abundance,
		top = top,
		transform = transform
	)
}

#' keep_variable
#' @inheritParams keep_variable
#'
#' @docType methods
#' @rdname keep_variable-methods
#'
#' @return A consistent object (to the input) with additional columns for 
#' the statistics from the hypothesis test (e.g.,  log fold change, p-value 
#' and false discovery rate).
setMethod("keep_variable", "spec_tbl_df", .keep_variable)

#' keep_variable
#' @inheritParams keep_variable
#'
#' @docType methods
#' @rdname keep_variable-methods
#'
#' @return A consistent object (to the input) with additional columns for the 
#' statistics from the hypothesis test (e.g.,  log fold change, p-value and 
#' false discovery rate).
setMethod("keep_variable", "tbl_df", .keep_variable)

#' keep_variable
#' @inheritParams keep_variable
#'
#' @docType methods
#' @rdname keep_variable-methods
#'
#' @return A consistent object (to the input) with additional columns for the 
#' statistics from the hypothesis test (e.g.,  log fold change, p-value and 
#' false discovery rate).
setMethod("keep_variable", "tidybulk", .keep_variable)

#' find abundant transcripts
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description identify_abundant() takes as input A `tbl` (with at least three 
#' columns for sample, feature and transcript abundance) or 
#' `SummarizedExperiment` (more convenient if abstracted to tibble with 
#' library(tidySummarizedExperiment)) and returns a consistent object 
#' (to the input) with additional columns for the statistics from the hypothesis test.
#'
#' @importFrom rlang enquo
#' @importFrom dplyr filter
#' @importFrom tidyr drop_na
#'
#' @name identify_abundant
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and 
#' transcript abundance) or `SummarizedExperiment` (more convenient if 
#' abstracted to tibble with library(tidySummarizedExperiment))
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param factor_of_interest The name of the column of the factor of interest. 
#' This is used for defining sample groups for the filtering process. It uses 
#' the filterByExpr function from edgeR.
#' @param minimum_counts A real positive number. It is the threshold of count 
#' per million that is used to filter transcripts/genes out from the scaling procedure.
#' @param minimum_proportion A real positive number between 0 and 1. It is the 
#' threshold of proportion of samples for each transcripts/genes that have to 
#' be characterised by a cmp bigger than the threshold to be included for 
#' scaling procedure.
#'
#' @details At the moment this function uses edgeR (DOI: 10.1093/bioinformatics/btp616)
#'
#'  Underlying method:
#'  edgeR::filterByExpr(
#'    data,
#'		min.count = minimum_counts,
#'		group = string_factor_of_interest,
#'		min.prop = minimum_proportion
#'	)
#'
#' @return A consistent object (to the input) with additional columns for the 
#' statistics from the hypothesis test (e.g.,  log fold change, p-value and 
#' false discovery rate).
#'
#' @examples
#'
#' 	identify_abundant(
#' 	tidybulk::se_mini
#' 	)
#'
#'
#' @docType methods
#' @rdname identify_abundant-methods
#' @export
#'
setGeneric("identify_abundant", function(.data,
																		 .sample = NULL,
																		 .transcript = NULL,
																		 .abundance = NULL,
																		 factor_of_interest = NULL,
																		 minimum_counts = 10,
																		 minimum_proportion = 0.7)
	standardGeneric("identify_abundant"))

# Set internal
.identify_abundant <- function(.data,
														.sample = NULL,
														.transcript = NULL,
														.abundance = NULL,
														factor_of_interest = NULL,
														minimum_counts = 10,
														minimum_proportion = 0.7) {

  # Fix NOTEs
  . = NULL

  factor_of_interest = enquo(factor_of_interest)

	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
	.sample = col_names$.sample
	.transcript = col_names$.transcript
	.abundance = col_names$.abundance

	if (minimum_counts < 0)
	  stop("The parameter minimum_counts must be > 0")
	if (minimum_proportion < 0 |	minimum_proportion > 1)
	  stop("The parameter minimum_proportion must be between 0 and 1")


	# Validate data frame
	if(do_validate()) {
	validation(.data, !!.sample, !!.transcript, !!.abundance)
	warning_if_data_is_not_rectangular(.data, !!.sample, !!.transcript, !!.abundance)
	}


	if(".abundant" %in% colnames(.data)) return(.data |> reattach_internals(.data))


  # Check if package is installed, otherwise install
  if (find.package("edgeR", quiet = TRUE) %>% length %>% equals(0)) {
    message("Installing edgeR needed for differential transcript abundance analyses")
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install("edgeR", ask = FALSE)
  }

	# If character fail
	if(
	  !is.null(factor_of_interest) &&
	  !factor_of_interest |> quo_is_null() &&
	  !factor_of_interest |> quo_is_symbolic()
	) stop("tidybulk says: factor_of_interest must be symbolic (i.e. ",
	  "column name/s not surrounded by single or double quotes) and not a character.")


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
	      .data |>
	      select(factor_of_interest) |>
	      lapply(class) |>
	      unlist() |>
	      as.character()%in% c("numeric", "integer", "double") |>
	      any()
	    )
	      stop("tidybulk says: The factor(s) of interest must not include continuous variables (e.g., integer,numeric, double).")

	    string_factor_of_interest =
	      .data %>%
	      select(!!.sample, factor_of_interest) |>
	      distinct() |>
	      arrange(!!.sample) |>
	      select(factor_of_interest) |>
	      pull(1)


	} else {
	  string_factor_of_interest = NULL
	}

	gene_to_exclude =
	.data %>%
	  select(!!.sample,!!.transcript, !!.abundance) |>
	  spread(!!.sample, !!.abundance) |>

	  # Drop if transcript have missing value
	  drop_na() %>%

	  # If I don't have any transcript with all samples give meaningful error
	  when(
	    nrow(.) == 0 ~ stop("tidybulk says: you don't have any transcript that is in all samples. Please consider using impute_missing_abundance."),
	    ~ (.)
	  ) %>%

	  # Call edgeR
	  as_matrix(rownames = !!.transcript) |>
	  edgeR::filterByExpr(
	    min.count = minimum_counts,
	    group = string_factor_of_interest,
	    min.prop = minimum_proportion
	  ) %>%
	  not() |>
	  which() |>
	  names()

	.data |>
	  dplyr::mutate(.abundant := (!!.transcript %in% gene_to_exclude) |> not()) |>

	  # Attach attributes
	  reattach_internals(.data)

}

#' keep_abundant
#' @inheritParams identify_abundant
#'
#' @docType methods
#' @rdname identify_abundant-methods
#'
#' @return A consistent object (to the input) with additional columns for the 
#' statistics from the hypothesis test (e.g.,  log fold change, p-value and 
#' false discovery rate).
setMethod("identify_abundant", "spec_tbl_df", .identify_abundant)

#' identify_abundant
#' @inheritParams identify_abundant
#'
#' @docType methods
#' @rdname identify_abundant-methods
#'
#' @return A consistent object (to the input) with additional columns for the 
#' statistics from the hypothesis test (e.g.,  log fold change, p-value and 
#' false discovery rate).
setMethod("identify_abundant", "tbl_df", .identify_abundant)

#' identify_abundant
#' @inheritParams identify_abundant
#'
#' @docType methods
#' @rdname identify_abundant-methods
#'
#' @return A consistent object (to the input) with additional columns for the 
#' statistics from the hypothesis test (e.g.,  log fold change, p-value and 
#' false discovery rate).
setMethod("identify_abundant", "tidybulk", .identify_abundant)


#' Keep abundant transcripts
#'
#' \lifecycle{questioning}
#'
#' @description keep_abundant() takes as input A `tbl` (with at least three 
#' columns for sample, feature and transcript abundance) or 
#' `SummarizedExperiment` (more convenient if abstracted to tibble with 
#' library(tidySummarizedExperiment)) and returns a consistent object 
#' (to the input) with additional columns for the statistics from the 
#' hypothesis test.
#'
#' @importFrom rlang enquo
#' @importFrom dplyr filter
#'
#' @name keep_abundant
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and 
#' transcript abundance) or `SummarizedExperiment` (more convenient if 
#' abstracted to tibble with library(tidySummarizedExperiment))
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param factor_of_interest The name of the column of the factor of interest. 
#' This is used for defining sample groups for the filtering process. It uses 
#' the filterByExpr function from edgeR.
#' @param minimum_counts A real positive number. It is the threshold of count 
#' per million that is used to filter transcripts/genes out from the scaling 
#' procedure.
#' @param minimum_proportion A real positive number between 0 and 1. It is the 
#' threshold of proportion of samples for each transcripts/genes that have to 
#' be characterised by a cmp bigger than the threshold to be included for 
#' scaling procedure.
#'
#' @details At the moment this function uses edgeR (DOI: 10.1093/bioinformatics/btp616)
#'
#'  Underlying method:
#'  edgeR::filterByExpr(
#'    data,
#'		min.count = minimum_counts,
#'		group = string_factor_of_interest,
#'		min.prop = minimum_proportion
#'	)
#'
#' @return A consistent object (to the input) with additional columns for the 
#' statistics from the hypothesis test (e.g.,  log fold change, p-value and 
#' false discovery rate).
#'
#'
#' @examples
#' 	keep_abundant(
#' 	tidybulk::se_mini
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
.keep_abundant <- function(.data,
														.sample = NULL,
														.transcript = NULL,
														.abundance = NULL,
														factor_of_interest = NULL,
														minimum_counts = 10,
														minimum_proportion = 0.7) {
  # Fix NOTEs
  . = NULL

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
	if(do_validate()) {
	validation(.data, !!.sample, !!.transcript, !!.abundance)
	warning_if_data_is_not_rectangular(.data, !!.sample, !!.transcript, !!.abundance)
	}

	.data |>

		# Filter
		identify_abundant(
			.sample = !!.sample,
			.transcript = !!.transcript,
			.abundance = !!.abundance,
			factor_of_interest = !!factor_of_interest,
			minimum_counts = minimum_counts,
			minimum_proportion = minimum_proportion
		) |>
		dplyr::filter(.abundant) |>

		# Attach attributes
		reattach_internals(.data)
}

#' keep_abundant
#' @inheritParams keep_abundant
#'
#' @docType methods
#' @rdname keep_abundant-methods
#'
#' @return A consistent object (to the input) with additional columns for 
#' the statistics from the hypothesis test (e.g.,  log fold change, p-value 
#' and false discovery rate).
setMethod("keep_abundant", "spec_tbl_df", .keep_abundant)

#' keep_abundant
#' @inheritParams keep_abundant
#'
#' @docType methods
#' @rdname keep_abundant-methods
#'
#' @return A consistent object (to the input) with additional columns for the 
#' statistics from the hypothesis test (e.g.,  log fold change, p-value and 
#' false discovery rate).
setMethod("keep_abundant", "tbl_df", .keep_abundant)

#' keep_abundant
#' @inheritParams keep_abundant
#'
#' @docType methods
#' @rdname keep_abundant-methods
#'
#' @return A consistent object (to the input) with additional columns for the 
#' statistics from the hypothesis test (e.g.,  log fold change, p-value and 
#' false discovery rate).
setMethod("keep_abundant", "tidybulk", .keep_abundant)



#' analyse gene enrichment with EGSEA
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description test_gene_enrichment() takes as input a `tbl` (with at least 
#' three columns for sample, feature and transcript abundance) or 
#' `SummarizedExperiment` (more convenient if abstracted to tibble with 
#' library(tidySummarizedExperiment)) and returns a `tbl` of gene set information
#'
#' @importFrom rlang enquo
#'
#' @name test_gene_enrichment
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and 
#' transcript abundance) or `SummarizedExperiment` (more convenient if 
#' abstracted to tibble with library(tidySummarizedExperiment))
#' @param .formula A formula with no response variable, representing the 
#' desired linear model
#' @param .sample The name of the sample column
#' @param .entrez The ENTREZ ID of the transcripts/genes
#' @param .abundance The name of the transcript/gene abundance column
#' @param contrasts This parameter takes the format of the contrast parameter 
#' of the method of choice. For edgeR and limma-voom is a character vector. 
#' For DESeq2 is a list including a character vector of length three. The first 
#' covariate is the one the model is tested against (e.g., ~ factor_of_interest)
#' @param methods A character vector. One or 3 or more methods to use in the 
#' testing (currently EGSEA errors if 2 are used). Type EGSEA::egsea.base() to 
#' see the supported GSE methods.
#' @param gene_sets A character vector or a list. It can take one or more of 
#' the following built-in collections as a character vector: c("h", "c1", "c2", 
#' "c3", "c4", "c5", "c6", "c7", "kegg_disease", "kegg_metabolism", 
#' "kegg_signaling"), to be used with EGSEA buildIdx. c1 is human specific. 
#' Alternatively, a list of user-supplied gene sets can be provided, to be used 
#' with EGSEA buildCustomIdx. In that case, each gene set is a character vector
#' of Entrez IDs and the names of the list are the gene set names.
#' @param species A character. It can be human, mouse or rat.
#' @param cores An integer. The number of cores available
#'
#' @param method DEPRECATED. Please use methods.
#' @param .contrasts DEPRECATED - This parameter takes the format of the 
#' contrast parameter of the method of choice. For edgeR and limma-voom is a 
#' character vector. For DESeq2 is a list including a character vector of length 
#' three. The first covariate is the one the model is tested against (e.g., ~ factor_of_interest)
#'
#' @details This wrapper executes ensemble gene enrichment analyses of the 
#' dataset using EGSEA (DOI:0.12688/f1000research.12544.1)
#'
#' dge =
#' 	data |>
#' 	keep_abundant(
#' 		factor_of_interest = !!as.symbol(parse_formula(.formula)[[1]]),
#' 		!!.sample, !!.entrez, !!.abundance
#' 	) %>%
#'
#' 	# Make sure transcript names are adjacent
#' 	[...] %>%
#' 	as_matrix(rownames = !!.entrez) %>%
#' 	edgeR::DGEList(counts = .)
#'
#' idx =  buildIdx(entrezIDs = rownames(dge), species = species, msigdb.gsets = msigdb.gsets,
#'	               kegg.exclude = kegg.exclude)
#'
#' dge |>
#'
#' 	# Calculate weights
#' 	limma::voom(design, plot = FALSE) |>
#'
#' 	# Execute EGSEA
#' 	egsea(
#' 		contrasts = my_contrasts,
#' 		baseGSEAs = methods,
#' 		gs.annots = idx,
#' 		sort.by = "med.rank",
#' 		num.threads = cores,
#' 		report = FALSE
#' 	)
#'
#' @return A consistent object (to the input)
#'
#'
#'
#'
#' @examples
#' \dontrun{
#'
#' library(SummarizedExperiment)
#' se = tidybulk::se_mini
#' rowData( se)$entrez = rownames(se )
#' df_entrez = aggregate_duplicates(se,.transcript = entrez )
#'
#' library("EGSEA")
#'
#' 	test_gene_enrichment(
#'			df_entrez,
#'			~ condition,
#'			.sample = sample,
#'			.entrez = entrez,
#'			.abundance = count,
#'          methods = c("roast" , "safe", "gage"  ,  "padog" , "globaltest", "ora" ),
#'          gene_sets = c("h", "c1", "c2", "c3", "c4", "c5", "c6", 
#'           "c7", "kegg_disease", "kegg_metabolism", "kegg_signaling"),
#'			species="human",
#'			cores = 2
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
																					 contrasts = NULL,
																					 methods = c("camera" , "roast" , "safe", "gage", "padog" , "globaltest", "ora"),
																					 gene_sets = c("h", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "kegg_disease", "kegg_metabolism", "kegg_signaling"),
																					 species,
																					 cores = 10,
																					 # DEPRECATED
																						method = NULL,
																						.contrasts = NULL
																						)
	standardGeneric("test_gene_enrichment"))

# Set internal
#' @importFrom lifecycle deprecate_warn
.test_gene_enrichment <- function(
          .data,
					.formula,
					.sample = NULL,
					.entrez,
					.abundance = NULL,
					contrasts = NULL,
			    methods = c("camera" , "roast" , "safe", "gage", "padog", "globaltest", "ora" ),
					gene_sets = c("h", "c1", "c2", "c3", "c4", "c5", "c6", "c7", 
					              "kegg_disease", "kegg_metabolism", "kegg_signaling"),
					species,
					cores = 10,
					# DEPRECATED
					 method = NULL,
					 .contrasts = NULL
					 ) {

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

	# Make col names
	.sample = enquo(.sample)
	.abundance = enquo(.abundance)
	col_names = get_sample_counts(.data, .sample, .abundance)
	.sample = col_names$.sample
	.abundance = col_names$.abundance

	.entrez = enquo(.entrez)

	# Check that there are no entrez missing
	.data =
		.data %>%
		when(
			filter(., !!.entrez |> is.na()) |> nrow() |> gt(0) ~ {
				warning("tidybulk says: There are NA entrez IDs. Those genes will be filtered")
				filter(., !!.entrez |> is.na() |> not())
			},
			~ (.)
		)

	# Validate data frame
	if(do_validate()) {
	validation(.data, !!.sample, !!.entrez, !!.abundance)
	warning_if_data_is_not_rectangular(.data, !!.sample, !!.entrez, !!.abundance)
	}

	.data %>%

		# Filter abundant if performed
		when(
			".abundant" %in% colnames(.) ~ filter(., .abundant),
			~ {
				warning("tidybulk says: highly abundant transcripts were not identified ",
				"(i.e. identify_abundant()) or filtered (i.e., keep_abundant), therefore ",
				"this operation will be performed on unfiltered data. In rare occasions this ",
				"could be wanted. In standard whole-transcriptome workflows is generally unwanted.")
				(.)
			}
		) |>

		test_gene_enrichment_bulk_EGSEA(
			.formula,
			.sample = !!.sample,
			.entrez = !!.entrez,
			.abundance = !!.abundance,
			.contrasts = contrasts,
			methods = methods,
			gene_sets = gene_sets,
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
#' @return A consistent object (to the input)
setMethod("test_gene_enrichment",
					"spec_tbl_df",
					.test_gene_enrichment)

#' test_gene_enrichment
#' @inheritParams test_gene_enrichment
#'
#' @docType methods
#' @rdname test_gene_enrichment-methods
#'
#' @return A consistent object (to the input)
setMethod("test_gene_enrichment",
					"tbl_df",
					.test_gene_enrichment)

#' test_gene_enrichment
#' @inheritParams test_gene_enrichment
#'
#' @docType methods
#' @rdname test_gene_enrichment-methods
#'
#' @return A consistent object (to the input)
setMethod("test_gene_enrichment",
					"tidybulk",
					.test_gene_enrichment)

#' analyse gene over-representation with GSEA
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description test_gene_overrepresentation() takes as input a `tbl` 
#' (with at least three columns for sample, feature and transcript abundance) 
#' or `SummarizedExperiment` (more convenient if abstracted to tibble with 
#' library(tidySummarizedExperiment)) and returns a `tbl` with the GSEA statistics
#'
#' @importFrom rlang enquo
#' @importFrom rlang quo_is_missing
#'
#' @name test_gene_overrepresentation
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and 
#' transcript abundance) or `SummarizedExperiment` (more convenient if 
#' abstracted to tibble with library(tidySummarizedExperiment))
#' @param .sample The name of the sample column
#' @param .entrez The ENTREZ ID of the transcripts/genes
#' @param .do_test A boolean column name symbol. It indicates the transcript to check
#' @param species A character. For example, human or mouse. MSigDB uses the 
#' latin species names (e.g., \"Mus musculus\", \"Homo sapiens\")
#' @param gene_sets  A character vector. The subset of MSigDB datasets you want 
#' to test against (e.g. \"C2\"). If NULL all gene sets are used (suggested). 
#' This argument was added to avoid time overflow of the examples.
#'
#' @param gene_set DEPRECATED. Use gene_sets instead.
#'
#' @details This wrapper execute gene enrichment analyses of the dataset using 
#' a list of transcripts and GSEA.
#' This wrapper uses clusterProfiler (DOI: doi.org/10.1089/omi.2011.0118) on the back-end.
#'
#' Undelying method:
#'  msigdbr::msigdbr(species = species) |>
#' 	nest(data = -gs_cat) |>
#' 	mutate(test =
#' 			map(
#' 				data,
#' 				~ clusterProfiler::enricher(
#' 					my_entrez_rank,
#' 				 	TERM2GENE=.x |> select(gs_name, entrez_gene),
#' 					pvalueCutoff = 1
#' 					) |>	as_tibble()
#' 			))
#'
#' @return A consistent object (to the input)
#'
#'
#'
#'
#' @examples
#'
#' print("Not run for build time.")
#'
#' #se_mini = aggregate_duplicates(tidybulk::se_mini, .transcript = entrez)
#' #df_entrez = mutate(df_entrez, do_test = feature %in% c("TNFRSF4", "PLCH2", "PADI4", "PAX7"))
#'
#' \dontrun{
#' 	test_gene_overrepresentation(
#' 		df_entrez,
#' 		.sample = sample,
#' 		.entrez = entrez,
#' 		.do_test = do_test,
#' 		species="Homo sapiens",
#'    gene_sets =c("C2")
#' 	)
#' }
#'
#' @docType methods
#' @rdname test_gene_overrepresentation-methods
#' @export
#'
#'
setGeneric("test_gene_overrepresentation", function(.data,
																										.entrez,
																										.do_test,
																										species,
																										.sample = NULL,
																										gene_sets  = NULL,
																										gene_set = NULL # DEPRECATED
																										)
	standardGeneric("test_gene_overrepresentation"))

# Set internal
.test_gene_overrepresentation <- function(.data,
																					 .entrez,
																					 .do_test,
																					 species,
																					 .sample = NULL,
																					 gene_sets  = NULL,
																					 gene_set = NULL  # DEPRECATED
																					 )	{
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
	if (.data %>% mutate(my_do_test = !!.do_test) %>% pull(my_do_test) |> is("logical") |> not() )
		stop("tidybulk says: .do_test column must be logical (i.e., TRUE or FALSE)")

	# Check packages msigdbr
	# Check if package is installed, otherwise install
	if (find.package("msigdbr", quiet = TRUE) |> length() |> equals(0)) {
		message("msigdbr not installed. Installing.")
		BiocManager::install("msigdbr", ask = FALSE)
	}

	# Check is correct species name
	if(species %in% msigdbr::msigdbr_species()$species_name |> not())
		stop(sprintf("tidybulk says: wrong species name. MSigDB uses the latin species names (e.g., %s)", 
		             paste(msigdbr::msigdbr_species()$species_name, collapse=", ")))

	.data |>
		#filter(!!.entrez %in% unique(m_df$entrez_gene)) |>
		filter(!!.do_test) |>
		distinct(!!.entrez) |>
		pull(!!.entrez) |>
		entrez_over_to_gsea(species, gene_collections  = gene_sets ) |>

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
#' @return A `spec_tbl_df` object
setMethod("test_gene_overrepresentation",
					"spec_tbl_df",
					.test_gene_overrepresentation)

#' test_gene_overrepresentation
#' @inheritParams test_gene_overrepresentation
#'
#' @docType methods
#' @rdname test_gene_overrepresentation-methods
#'
#' @return A `tbl_df` object
setMethod("test_gene_overrepresentation",
					"tbl_df",
					.test_gene_overrepresentation)

#' test_gene_overrepresentation
#' @inheritParams test_gene_overrepresentation
#'
#' @docType methods
#' @rdname test_gene_overrepresentation-methods
#'
#' @return A `tidybulk` object
setMethod("test_gene_overrepresentation",
					"tidybulk",
					.test_gene_overrepresentation)

#' analyse gene rank with GSEA
#'
#' \lifecycle{maturing}
#'
#' @description test_gene_rank() takes as input a `tbl` (with at least three 
#' columns for sample, feature and transcript abundance) or 
#' `SummarizedExperiment` (more convenient if abstracted to tibble with 
#' library(tidySummarizedExperiment)) and returns a `tbl` with the GSEA statistics
#'
#' @importFrom rlang enquo
#' @importFrom rlang quo_is_missing
#'
#' @name test_gene_rank
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and 
#' transcript abundance) or `SummarizedExperiment` (more convenient if 
#' abstracted to tibble with library(tidySummarizedExperiment))
#' @param .sample The name of the sample column
#' @param .entrez The ENTREZ ID of the transcripts/genes
#' @param .arrange_desc A column name of the column to arrange in decreasing order
#' @param species A character. For example, human or mouse. MSigDB uses the 
#' latin species names (e.g., \"Mus musculus\", \"Homo sapiens\")
#' @param gene_sets A character vector or a list. It can take one or more of 
#' the following built-in collections as a character vector: c("h", "c1", "c2", 
#' "c3", "c4", "c5", "c6", "c7", "kegg_disease", "kegg_metabolism", 
#' "kegg_signaling"), to be used with EGSEA buildIdx. c1 is human specific. 
#' Alternatively, a list of user-supplied gene sets can be provided, to be used 
#' with EGSEA buildCustomIdx. In that case, each gene set is a character vector 
#' of Entrez IDs and the names of the list are the gene set names.
#'
#' @param gene_set DEPRECATED. Use gene_sets instead.
#'
#' @details This wrapper execute gene enrichment analyses of the dataset using 
#' a list of transcripts and GSEA.
#' This wrapper uses clusterProfiler (DOI: doi.org/10.1089/omi.2011.0118) on the back-end.
#'
#' Undelying method:
#'# Get gene sets signatures
#'msigdbr::msigdbr(species = species) %>%
#'
#'	# Filter specific gene_sets  if specified. This was introduced to speed up 
#'	examples executionS
#'	when(
#'		!is.null(gene_sets ) ~ filter(., gs_cat %in% gene_sets ),
#'		~ (.)
#'	) |>
#'
#'	# Execute calculation
#'	nest(data = -gs_cat) |>
#'	mutate(fit =
#'				 	map(
#'				 		data,
#'				 		~ 	clusterProfiler::GSEA(
#'				 			my_entrez_rank,
#'				 			TERM2GENE=.x |> select(gs_name, entrez_gene),
#'				 			pvalueCutoff = 1
#'				 		)
#'
#'				 	))
#'
#' @return A consistent object (to the input)
#'
#'
#' @examples
#'
#' print("Not run for build time.")
#'
#' \dontrun{
#'
#' df_entrez = tidybulk::se_mini
#' df_entrez = mutate(df_entrez, do_test = .feature %in% c("TNFRSF4", "PLCH2", "PADI4", "PAX7"))
#' df_entrez  = df_entrez |> test_differential_abundance(~ condition)
#'
#'
#'	test_gene_rank(
#'		df_entrez,
#' 		.sample = .sample,
#'		.entrez = entrez,
#' 		species="Homo sapiens",
#'    gene_sets =c("C2"),
#'  .arrange_desc = logFC
#' 	)
#' }
#'
#' @docType methods
#' @rdname test_gene_rank-methods
#' @export
#'
#'
setGeneric("test_gene_rank", function(.data,
																			.entrez,
																			.arrange_desc,
																			species,
																			.sample = NULL,
																			gene_sets  = NULL,
																			gene_set = NULL  # DEPRECATED
																			)
	standardGeneric("test_gene_rank"))

# Set internal
.test_gene_rank <- function(.data,
														 .entrez,
														 .arrange_desc,
														 species,
														 .sample = NULL,
														 gene_sets = c("h", "c1", "c2", "c3", "c4", "c5", "c6", "c7"),
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

	# DEPRECATION OF reference function
	if (is_present(.sample) & !is.null(.sample)) {
	  
	  # Signal the deprecation to the user
	  deprecate_warn("1.13.2", "tidybulk::test_gene_rank(.sample = )", 
	                 details = "The argument .sample is now deprecated and not needed anymore.")

	}
	
	# Get column names
	.arrange_desc = enquo(.arrange_desc)
	.entrez = enquo(.entrez)

	# Check if ranking is set
	if(quo_is_missing(.arrange_desc))
		stop("tidybulk says: the .arrange_desc parameter appears to no be set")

	# Check if entrez is set
	if(quo_is_missing(.entrez))
		stop("tidybulk says: the .entrez parameter appears to no be set")

	# Check packages msigdbr
	# Check if package is installed, otherwise install
	if (find.package("msigdbr", quiet = TRUE) |> length() |> equals(0)) {
		message("msigdbr not installed. Installing.")
		BiocManager::install("msigdbr", ask = FALSE)
	}

	# Check is correct species name
	if(species %in% msigdbr::msigdbr_species()$species_name |> not())
		stop(sprintf("tidybulk says: wrong species name. MSigDB uses the latin species names (e.g., %s)", paste(msigdbr::msigdbr_species()$species_name, collapse=", ")))

	# Check if missing entrez
	if(.data |> filter(is.na(!!.entrez)) |> nrow() > 0 ){
		warning("tidybulk says: there are .entrez that are NA. Those will be removed")
		.data = .data |>	filter(is.na(!!.entrez) |> not())
	}

	# Check if missing .arrange_desc
	if(.data |> filter(is.na(!!.arrange_desc)) |> nrow() > 0 ){
		warning("tidybulk says: there are .arrange_desc that are NA. Those will be removed")
		.data = .data |>	filter(is.na(!!.arrange_desc ) |> not())
	}

	.data |>
		select(!!.entrez, !!.arrange_desc) |>
	  distinct() |> 
	  
	  # Select one entrez - NEEDED?
	  with_groups(c(!!.entrez,!!.arrange_desc ), slice, 1) |> 

	  # arrange 
	  arrange(desc(!!.arrange_desc)) |>
	  
	  # Format
		deframe() |>
		entrez_rank_to_gsea(species, gene_collections  = gene_sets ) |>

	  # Add methods used. It is here and not in fucntions because I need the original .data
	  memorise_methods_used(c("clusterProfiler", "enrichplot"), object_containing_methods = .data) %>%
	  when(
	    gene_sets |> is("character") ~ (.) |> memorise_methods_used("msigdbr"),
	    ~ (.)
	  )


}

#' test_gene_rank
#' @inheritParams test_gene_rank
#'
#' @docType methods
#' @rdname test_gene_rank-methods
#'
#' @return A `spec_tbl_df` object
setMethod("test_gene_rank",
					"spec_tbl_df",
					.test_gene_rank)

#' test_gene_rank
#' @inheritParams test_gene_rank
#'
#' @docType methods
#' @rdname test_gene_rank-methods
#'
#' @return A `tbl_df` object
setMethod("test_gene_rank",
					"tbl_df",
					.test_gene_rank)

#' test_gene_rank
#' @inheritParams test_gene_rank
#'
#' @docType methods
#' @rdname test_gene_rank-methods
#'
#' @return A `tidybulk` object
setMethod("test_gene_rank",
					"tidybulk",
					.test_gene_rank)


#' Extract sample-wise information
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description pivot_sample() takes as input a `tbl` 
#' (with at least three columns for sample, feature and transcript abundance) 
#' or `SummarizedExperiment` (more convenient if abstracted to tibble with 
#' library(tidySummarizedExperiment)) and returns a `tbl` with only sample-related columns
#'
#'
#' @name pivot_sample
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and 
#' transcript abundance) or `SummarizedExperiment` (more convenient if 
#' abstracted to tibble with library(tidySummarizedExperiment))
#' @param .sample The name of the sample column
#'
#'
#' @details This functon extracts only sample-related information for 
#' downstream analysis (e.g., visualisation). It is disruptive in the sense 
#' that it cannot be passed anymore to tidybulk function.
#'
#' @return A `tbl` with transcript-related information
#'
#'
#'
#'
#' @examples
#'
#'
#' 	pivot_sample(tidybulk::se_mini )
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

	.data |>

		# Selecting the right columns
		select(
			!!.sample,
			get_specific_annotation_columns(.data, !!.sample)
		) |>
		distinct() |>

		drop_class(c("tidybulk", "tt")) |>
		drop_internals()


}

#' pivot_sample
#' @inheritParams pivot_sample
#'
#' @docType methods
#' @rdname pivot_sample-methods
#'
setMethod("pivot_sample",
					"spec_tbl_df",
					.pivot_sample)

#' pivot_sample
#' @inheritParams pivot_sample
#'
#' @docType methods
#' @rdname pivot_sample-methods
#'
setMethod("pivot_sample",
					"tbl_df",
					.pivot_sample)

#' pivot_sample
#' @inheritParams pivot_sample
#'
#' @docType methods
#' @rdname pivot_sample-methods
#'
setMethod("pivot_sample",
					"tidybulk",
					.pivot_sample)

#' Extract transcript-wise information
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description pivot_transcript() takes as input a `tbl` (with at least three 
#' columns for sample, feature and transcript abundance) or 
#' `SummarizedExperiment` (more convenient if abstracted to tibble with 
#' library(tidySummarizedExperiment)) and returns a `tbl` with only 
#' transcript-related columns
#'
#'
#' @name pivot_transcript
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and 
#' transcript abundance) or `SummarizedExperiment` (more convenient if 
#' abstracted to tibble with library(tidySummarizedExperiment))
#' @param .transcript The name of the transcript column
#'
#'
#' @details This functon extracts only transcript-related information for 
#' downstream analysis (e.g., visualisation). It is disruptive in the sense 
#' that it cannot be passed anymore to tidybulk function.
#'
#' @return A `tbl` with transcript-related information
#'
#'
#'
#'
#' @examples
#'
#'
#' 	pivot_transcript(tidybulk::se_mini 	)
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

	.data |>

		# Selecting the right columns
		select(
			!!.transcript,
			get_specific_annotation_columns(.data, !!.transcript)
		) |>
		distinct() |>

		drop_class(c("tidybulk", "tt")) |>
		drop_internals()


}

#' pivot_transcript
#' @inheritParams pivot_transcript
#'
#' @docType methods
#' @rdname pivot_transcript-methods
#'
setMethod("pivot_transcript",
					"spec_tbl_df",
					.pivot_transcript)

#' pivot_transcript
#' @inheritParams pivot_transcript
#'
#' @docType methods
#' @rdname pivot_transcript-methods
#'
setMethod("pivot_transcript",
					"tbl_df",
					.pivot_transcript)

#' pivot_transcript
#' @inheritParams pivot_transcript
#'
#' @docType methods
#' @rdname pivot_transcript-methods
#'
setMethod("pivot_transcript",
					"tidybulk",
					.pivot_transcript)


#' Fill transcript abundance if missing from sample-transcript pairs
#'
#' \lifecycle{questioning}
#'
#' @description fill_missing_abundance() takes as input A `tbl` (with at least 
#' three columns for sample, feature and transcript abundance) or 
#' `SummarizedExperiment` (more convenient if abstracted to tibble with 
#' library(tidySummarizedExperiment)) and returns a consistent object 
#' (to the input) with new observations
#'
#' @importFrom rlang enquo
#'
#' @name fill_missing_abundance
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT>  | <...> |
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript column
#' @param .abundance The name of the transcript abundance column
#' @param fill_with A numerical abundance with which fill the missing data points
#'
#' @details This function fills the abundance of missing sample-transcript 
#' pair using the median of the sample group defined by the formula
#'
#' @return A consistent object (to the input) non-sparse abundance
#'
#'
#'
#'
#' @examples
#'
#' print("Not run for build time.")
#'
#' # tidybulk::se_mini |>  fill_missing_abundance( fill_with = 0)
#'
#'
#' @docType methods
#' @rdname fill_missing_abundance-methods
#'
#' @export
#'
#'
setGeneric("fill_missing_abundance", function(.data,
																		.sample= NULL,
																		.transcript= NULL,
																		.abundance= NULL,
																		fill_with)
	standardGeneric("fill_missing_abundance"))

# Set internal
.fill_missing_abundance = 	function(.data,
													.sample = NULL,
													.transcript= NULL,
													.abundance= NULL,
													fill_with)
{



	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
	.sample = col_names$.sample
	.transcript = col_names$.transcript
	.abundance = col_names$.abundance


	# Check the abundance is set
	if(length(fill_with)==0) stop("nanny says: the argument fill_with must not be empty.")

	# Validate data frame
	if(do_validate()) validation(.data, !!.sample, !!.transcript, !!.abundance)

	fill_NA_using_value(
		.data,
		.sample = !!.sample,
		.transcript = !!.transcript,
		.abundance = !!.abundance,
		fill_with = fill_with)
}

#' fill_missing_abundance
#' @inheritParams fill_missing_abundance
#'
#' @docType methods
#' @rdname fill_missing_abundance-methods
#'
#' @return A consistent object (to the input) with filled abundance
setMethod("fill_missing_abundance", "spec_tbl_df", .fill_missing_abundance)

#' fill_missing_abundance
#' @inheritParams fill_missing_abundance
#'
#' @docType methods
#' @rdname fill_missing_abundance-methods
#'
#' @return A consistent object (to the input) with filled abundance
setMethod("fill_missing_abundance", "tbl_df", .fill_missing_abundance)

#' fill_missing_abundance
#' @inheritParams fill_missing_abundance
#'
#' @docType methods
#' @rdname fill_missing_abundance-methods
#'
#' @return A consistent object (to the input) with filled abundance
setMethod("fill_missing_abundance", "tidybulk", .fill_missing_abundance)



#' impute transcript abundance if missing from sample-transcript pairs
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description impute_missing_abundance() takes as input A `tbl` (with at least 
#' three columns for sample, feature and transcript abundance) or 
#' `SummarizedExperiment` (more convenient if abstracted to tibble with 
#' library(tidySummarizedExperiment)) and returns a consistent object 
#' (to the input) with additional sample-transcript pairs with imputed 
#' transcript abundance.
#'
#' @importFrom rlang enquo
#'
#' @name impute_missing_abundance
#'
#' @param .data A `tbl` (with at least three columns for sample, feature 
#' and transcript abundance) or `SummarizedExperiment` (more convenient if 
#' abstracted to tibble with library(tidySummarizedExperiment))
#' @param .formula A formula with no response variable, representing the 
#' desired linear model where the first covariate is the factor of interest 
#' and the second covariate is the unwanted variation (of the kind ~ factor_of_interest + batch)
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param suffix A character string. This is added to the imputed count column 
#' names. If empty the count column are overwritten
#' @param force_scaling A boolean. In case a abundance-containing column is not 
#' scaled (columns with _scale suffix), setting force_scaling = TRUE will result 
#' in a scaling by library size, to compensating for a possible difference in 
#' sequencing depth.
#'
#' @details This function imputes the abundance of missing sample-transcript 
#' pair using the median of the sample group defined by the formula
#'
#' @return A consistent object (to the input) non-sparse abundance
#'
#'
#'
#'
#' @examples
#'
#'
#' res =
#' 	impute_missing_abundance(
#' 		tidybulk::se_mini,
#' 	~ condition
#' )
#'
#'
#' @docType methods
#' @rdname impute_missing_abundance-methods
#'
#' @export
#'
#'
setGeneric("impute_missing_abundance", function(.data,
																				.formula,
																				.sample = NULL,
																				.transcript = NULL,
																				.abundance = NULL,
																				suffix = "",
																				force_scaling = FALSE)
	standardGeneric("impute_missing_abundance"))

# Set internal
.impute_missing_abundance <- function(.data,
															.formula,
															.sample = NULL,
															.transcript = NULL,
															.abundance = NULL,
															suffix = "",
															force_scaling = FALSE)
{

  # Fix NOTEs
  . = NULL

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
		.data |> get_tt_columns() |> is.null() |> not() &&
		".abundance_scaled" %in% (.data |> get_tt_columns() |> names()) &&
		quo_name(.data |> get_tt_columns() %$% .abundance_scaled) %in% (.data |> colnames()) &&
		quo_name(.data |> get_tt_columns() %$% .abundance_scaled) != quo_name(.abundance)
	)
		.abundance_scaled = get_tt_columns(.data)$.abundance_scaled

	# Validate data frame
	if(do_validate())  validation(.data, !!.sample, !!.transcript, !!.abundance)

	fill_NA_using_formula(
			.data,
			.formula,
			.sample = !!.sample,
			.transcript = !!.transcript,
			.abundance = !!.abundance,
			.abundance_scaled = !!.abundance_scaled,
			suffix = suffix,
			force_scaling = force_scaling) |>

		# Reattach internals
		reattach_internals(.data)

}

#' impute_missing_abundance
#' @inheritParams impute_missing_abundance
#'
#' @docType methods
#' @rdname impute_missing_abundance-methods
#'
#' @return A consistent object (to the input) with imputed abundance
setMethod("impute_missing_abundance", "spec_tbl_df", .impute_missing_abundance)

#' impute_missing_abundance
#' @inheritParams impute_missing_abundance
#'
#' @docType methods
#' @rdname impute_missing_abundance-methods
#'
#' @return A consistent object (to the input) with imputed abundance
setMethod("impute_missing_abundance", "tbl_df", .impute_missing_abundance)

#' impute_missing_abundance
#' @inheritParams impute_missing_abundance
#'
#' @docType methods
#' @rdname impute_missing_abundance-methods
#'
#' @return A consistent object (to the input) with imputed abundance
setMethod("impute_missing_abundance", "tidybulk", .impute_missing_abundance)




#' Add differential tissue composition information to a tbl
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description test_differential_cellularity() takes as input A `tbl` (with at 
#' least three columns for sample, feature and transcript abundance) or 
#' `SummarizedExperiment` (more convenient if abstracted to tibble with 
#' library(tidySummarizedExperiment)) and returns a consistent object 
#' (to the input) with additional columns for the statistics from the hypothesis test.
#'
#' @importFrom rlang enquo
#' @importFrom stringr str_detect
#'
#' @name test_differential_cellularity
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and 
#' transcript abundance) or `SummarizedExperiment` (more convenient if 
#' abstracted to tibble with library(tidySummarizedExperiment))
#' @param .formula A formula representing the desired linear model. The formula 
#' can be of two forms: multivariable (recommended) or univariable 
#' Respectively: \"factor_of_interest ~ .\" or \". ~ factor_of_interest\". 
#' The dot represents cell-type proportions, and it is mandatory. If censored 
#' regression is desired (coxph) the formula should be of the form \"survival::Surv\(y, dead\) ~ .\"
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param method A string character. Either \"cibersort\", \"epic\" or \"llsr\". 
#' The regression method will be chosen based on being multivariable: lm or 
#' cox-regression (both on logit-transformed proportions); or univariable: 
#' beta or cox-regression (on logit-transformed proportions). See .formula 
#' for multi- or univariable choice.
#' @param reference A data frame. The transcript/cell_type data frame of 
#' integer transcript abundance
#' @param significance_threshold A real between 0 and 1 (usually 0.05).
#' @param ... Further parameters passed to the method deconvolve_cellularity
#'
#' @details This routine applies a deconvolution method (e.g., Cibersort; DOI: 10.1038/nmeth.3337)
#' and passes the proportions inferred into a generalised linear model 
#' (DOI:dx.doi.org/10.1007/s11749-010-0189-z)
#' or a cox regression model (ISBN: 978-1-4757-3294-8)
#'
#' Underlying method for the generalised linear model:
#' data |>
#' deconvolve_cellularity(
#' 	!!.sample, !!.transcript, !!.abundance,
#' 	method=method,
#' 	reference = reference,
#' 	action="get",
#' 	...
#' )  %>%
#' 	[..] %>%
#' 	betareg::betareg(.my_formula, .)
#'
#' Underlying method for the cox regression:
#' data |>
#' deconvolve_cellularity(
#' 	!!.sample, !!.transcript, !!.abundance,
#' 	method=method,
#' 	reference = reference,
#' 	action="get",
#' 	...
#' )  %>%
#' 	[..] %>%
#' 	mutate(.proportion_0_corrected = .proportion_0_corrected  |> boot::logit()) %>%
#' 	survival::coxph(.my_formula, .)
#'
#' @return A consistent object (to the input) with additional columns for 
#' the statistics from the hypothesis test (e.g.,  log fold change, p-value 
#' and false discovery rate).
#'
#'
#'
#'
#' @examples
#'
#'  # Regular regression
#' 	test_differential_cellularity(
#' 	 tidybulk::se_mini ,
#' 	    . ~ condition,
#' 	    cores = 1
#' 	)
#'
#' 	# Cox regression - multiple
#'
#'	tidybulk::se_mini |>
#'
#'		# Test
#'		test_differential_cellularity(
#'		    survival::Surv(days, dead) ~ .,
#'		    cores = 1
#'		)
#'
#'
#'
#' @docType methods
#' @rdname test_differential_cellularity-methods
#' @export
#'
setGeneric("test_differential_cellularity", function(.data,
																										 .formula,
																										 .sample = NULL,
																										 .transcript = NULL,
																										 .abundance = NULL,
																										 method = "cibersort",
																										 reference = X_cibersort,
																										 significance_threshold = 0.05,
																										 ...)
					 standardGeneric("test_differential_cellularity"))

# Set internal
.test_differential_cellularity <- function(.data,
																						.formula,
																						.sample = NULL,
																						.transcript = NULL,
																						.abundance = NULL,
																						method = "cibersort",
																						reference = X_cibersort,
																						significance_threshold = 0.05,
																						...) {
  # Fix NOTEs
  . = NULL

	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
	.sample = col_names$.sample
	.transcript = col_names$.transcript
	.abundance = col_names$.abundance

	# Validate data frame
	if(do_validate()) validation(.data, !!.sample, !!.transcript, !!.abundance)

	# Validate formula
	if(.formula |> format() |> str_detect(" \\.|\\. ", negate = TRUE))
		stop("tidybulk says: in the formula a dot must be present in either these ",
		     "forms \". ~\" or \"~ .\" with a white-space after or before respectively")

	test_differential_cellularity_(
		.data,
		.formula = .formula,
		.sample = !!.sample,
		.transcript = !!.transcript,
		.abundance = !!.abundance,
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
setMethod("test_differential_cellularity",
					"spec_tbl_df",
					.test_differential_cellularity)

#' test_differential_cellularity
#' @inheritParams test_differential_cellularity
#'
#' @docType methods
#' @rdname test_differential_cellularity-methods
#'
setMethod("test_differential_cellularity",
					"tbl_df",
					.test_differential_cellularity)

#' test_differential_cellularity
#' @inheritParams test_differential_cellularity
#'
#' @docType methods
#' @rdname test_differential_cellularity-methods
#'
setMethod("test_differential_cellularity",
					"tidybulk",
					.test_differential_cellularity)

#' Test of stratification of biological replicates based on tissue composition, 
#' one cell-type at the time, using Kaplan-meier curves.
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description test_stratification_cellularity() takes as input A `tbl` 
#' (with at least three columns for sample, feature and transcript abundance) 
#' or `SummarizedExperiment` (more convenient if abstracted to tibble with 
#' library(tidySummarizedExperiment)) and returns a consistent object 
#' (to the input) with additional columns for the statistics from the hypothesis test.
#'
#' @importFrom rlang enquo
#' @importFrom stringr str_detect
#'
#' @name test_stratification_cellularity
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and 
#' transcript abundance) or `SummarizedExperiment` (more convenient if 
#' abstracted to tibble with library(tidySummarizedExperiment))
#' @param .formula A formula representing the desired linear model. The formula 
#' can be of two forms: multivariable (recommended) or univariable 
#' Respectively: \"factor_of_interest ~ .\" or \". ~ factor_of_interest\". 
#' The dot represents cell-type proportions, and it is mandatory. If censored 
#' regression is desired (coxph) the formula should be of the form \"survival::Surv\(y, dead\) ~ .\"
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param method A string character. Either \"cibersort\", \"epic\" or \"llsr\". 
#' The regression method will be chosen based on being multivariable: lm or 
#' cox-regression (both on logit-transformed proportions); or univariable: beta 
#' or cox-regression (on logit-transformed proportions). See .formula for 
#' multi- or univariable choice.
#' @param reference A data frame. The transcript/cell_type data frame of integer 
#' transcript abundance
#' @param ... Further parameters passed to the method deconvolve_cellularity
#'
#' @details This routine applies a deconvolution method (e.g., Cibersort; DOI: 10.1038/nmeth.3337)
#' and passes the proportions inferred into a generalised linear model 
#' (DOI:dx.doi.org/10.1007/s11749-010-0189-z)
#' or a cox regression model (ISBN: 978-1-4757-3294-8)
#'
#'
#' Underlying method for the test:
#' data |>
#' deconvolve_cellularity(
#' 	!!.sample, !!.transcript, !!.abundance,
#' 	method=method,
#' 	reference = reference,
#' 	action="get",
#' 	...
#' )  %>%
#' 	[..] |>
#' 	mutate(.high_cellularity = .proportion > median(.proportion)) |>
#' 	survival::survdiff(data = data, .my_formula)
#'
#' @return A consistent object (to the input) with additional columns for the 
#' statistics from the hypothesis test (e.g.,  log fold change, p-value and 
#' false discovery rate).
#'

#' @examples
#'	tidybulk::se_mini |>
#'	test_stratification_cellularity(
#'		survival::Surv(days, dead) ~ .,
#'		cores = 1
#'	)
#'
#'
#'
#' @docType methods
#' @rdname test_stratification_cellularity-methods
#' @export
#'
setGeneric("test_stratification_cellularity", function(.data,
																											 .formula,
																											 .sample = NULL,
																											 .transcript = NULL,
																											 .abundance = NULL,
																											 method = "cibersort",
																											 reference = X_cibersort,
																											 ...)
	standardGeneric("test_stratification_cellularity"))

# Set internal
.test_stratification_cellularity <- function(.data,
																							.formula,
																							.sample = NULL,
																							.transcript = NULL,
																							.abundance = NULL,
																							method = "cibersort",
																							reference = X_cibersort,
																							...) {
  # Fix NOTEs
  . = NULL

	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
	.sample = col_names$.sample
	.transcript = col_names$.transcript
	.abundance = col_names$.abundance

	# Validate data frame
	if(do_validate()) validation(.data, !!.sample, !!.transcript, !!.abundance)

	# Validate formula
	if(.formula |> format() %>% str_detect(" \\.|\\. ", negate = TRUE))
		stop("tidybulk says: in the formula a dot must be present in either ",
		     "these forms \". ~\" or \"~ .\" with a white-space after or before respectively")

	test_stratification_cellularity_(
		.data,
		.formula = .formula,
		.sample = !!.sample,
		.transcript = !!.transcript,
		.abundance = !!.abundance,
		method = method,
		reference = reference,
		...
	)

}

#' test_stratification_cellularity
#' @inheritParams test_stratification_cellularity
#'
#' @docType methods
#' @rdname test_stratification_cellularity-methods
#'
setMethod("test_stratification_cellularity",
					"spec_tbl_df",
					.test_stratification_cellularity)

#' test_stratification_cellularity
#' @inheritParams test_stratification_cellularity
#'
#' @docType methods
#' @rdname test_stratification_cellularity-methods
#'
setMethod("test_stratification_cellularity",
					"tbl_df",
					.test_stratification_cellularity)

#' test_stratification_cellularity
#' @inheritParams test_stratification_cellularity
#'
#' @docType methods
#' @rdname test_stratification_cellularity-methods
#'
setMethod("test_stratification_cellularity",
					"tidybulk",
					.test_stratification_cellularity)



#' Produces the bibliography list of your workflow
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description get_bibliography() takes as input a `tidybulk`
#'
#' @importFrom rlang enquo
#'
#' @name get_bibliography
#'
#' @param .data A `tbl` (with at least three columns for sample, feature 
#' and transcript abundance) or `SummarizedExperiment` (more convenient if 
#' abstracted to tibble with library(tidySummarizedExperiment))
#'
#' @details This methods returns the bibliography list of your workflow from 
#' the internals of a tidybulk object (attr(., "internals"))
#'
#'
#' @examples
#' get_bibliography(tidybulk::se_mini)
#'
#'
#' @docType methods
#' @rdname get_bibliography-methods
#'
#' @return NULL. It prints a list of bibliography references for the software 
#' used through the workflow.
#' @export
#'
setGeneric("get_bibliography", function(.data)
	standardGeneric("get_bibliography"))

# Set internal
.get_bibliography = 		function(.data)
{

  # Fix NOTEs
  . = NULL

	default_methods = c("tidybulk", "tidyverse")

	# If there is not attributes parameter
	my_methods =
			.data %>%
			when(
				!(
					!"internals" %in% (attributes(.) |> names()) &&
						!"methods_used" %in% (attr(., "internals") |> names())
				) ~ 	attr(., "internals") %>% .[["methods_used"]],
				~ ""
			)


	my_bibliography() %>%
		.[c(default_methods, my_methods)] |>
		unlist() |>
		writeLines()

}

#' get_bibliography
#' @inheritParams get_bibliography
#'
#' @docType methods
#' @rdname get_bibliography-methods
#'
setMethod("get_bibliography",
					"tbl",
					.get_bibliography)

#' get_bibliography
#' @inheritParams get_bibliography
#'
#' @docType methods
#' @rdname get_bibliography-methods
#'
setMethod("get_bibliography",
					"tbl_df",
					.get_bibliography)

#' get_bibliography
#' @inheritParams get_bibliography
#'
#' @docType methods
#' @rdname get_bibliography-methods
#'
setMethod("get_bibliography",
					"spec_tbl_df",
					.get_bibliography)

#' get_bibliography
#' @inheritParams get_bibliography
#'
#' @docType methods
#' @rdname get_bibliography-methods
#'
setMethod("get_bibliography",
					"tidybulk",
					.get_bibliography)

#' Get matrix from tibble
#'
#'
#'
#'
#' @importFrom magrittr set_rownames
#' @importFrom rlang quo_is_null
#'
#' @param tbl A tibble
#' @param rownames The column name of the input tibble that will become 
#' the rownames of the output matrix
#' @param do_check A boolean
#'
#' @return A matrix
#'
#' @examples
#'
#'
#' tibble(.feature = "CD3G", count=1) |> as_matrix(rownames=.feature)
#'
#' @export
as_matrix <- function(tbl,
                      rownames = NULL,
                      do_check = TRUE) {

  # Fix NOTEs
  . = NULL

  rownames = enquo(rownames)
  tbl %>%
    # Through warning if data frame is not numerical beside the rownames column (if present)
    ifelse_pipe(
      do_check &&
        tbl %>%
        # If rownames defined eliminate it from the data frame
        ifelse_pipe(!quo_is_null(rownames), ~ .x[,-1], ~ .x) |>
        dplyr::summarise_all(class) |>
        tidyr::gather(variable, class) |>
        pull(class) |>
        unique() %>%
        `%in%`(c("numeric", "integer")) |> not() |> any(),
      ~ {
        warning("tidybulk says: there are NON-numerical columns, ",
                "the matrix will NOT be numerical")
        .x
      }
    ) |>
    as.data.frame() |>

    # Deal with rownames column if present
    ifelse_pipe(
      !quo_is_null(rownames),
      ~ .x |>
        magrittr::set_rownames(tbl |> pull(!!rownames)) |>
        select(-1)
    ) |>

    # Convert to matrix
    as.matrix()
}

