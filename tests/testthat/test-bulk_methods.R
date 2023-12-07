context('Bulk methods')

library(dplyr)

data("se_mini")
data("breast_tcga_mini_SE")
input_df = se_mini |> tidybulk() |> as_tibble() |> setNames(c("b","a",  "c", "Cell type", "time" , "condition", "days",  "dead", "entrez"))
input_df_breast =   breast_tcga_mini_SE |> tidybulk() |> as_tibble() |> setNames(c( "b","a", "c", "c norm", "call"))

test_that("Creating tt object from tibble, number of parameters, methods",{

	expect_equal(

		length(
			attr(
				tidybulk(
					input_df,
					.sample = a,
					.transcript = b,
					.abundance = c
				) ,
				"internals"
			)$tt_columns
		),
		3
	)

})

test_that("Test class identity of tt object",{

	expect_equal(
		class(
			tidybulk(
				input_df,
				.sample = a,
				.transcript = b,
				.abundance = c
			)
		)[1],
		"tidybulk"
	)

})

test_that("Only scaled counts - no object",{

	res =
		scale_abundance(
			input_df |> identify_abundant(a, b, c),
			.sample = a,
			.transcript = b,
			.abundance = c,
			action = "only"
		)

	expect_equal(
		unique(res$multiplier),
		c(1.2994008, 1.1781297, 2.6996428, 0.9702628, 1.8290148),
		tolerance=1e-3
	)

	expect_equal(
		ncol(res),
		3
	)

	internals = attr(scale_abundance(tidybulk(input_df, a, b, c)), "internals")

	expect_equal(length(internals$tt_columns), 4 )

	expect_equal(rlang::quo_name(internals$tt_columns[[4]]), "c_scaled" )

	# With factor of interest
	res =
		scale_abundance(
			input_df |> identify_abundant(a, b, c, factor_of_interest = condition),
			.sample = a,
			.transcript = b,
			.abundance = c,
			action = "only"
		)

	expect_equal(
		unique(res$multiplier),
		c(1.3078113, 1.1929933, 1.9014731, 0.9678922, 1.4771970),
		tolerance=1e-3
	)

	expect_equal(
		ncol(res),
		3
	)

	# Warnings on continuous
	sam = distinct(input_df, a)
	sam = mutate(sam, condition_cont = c(-0.4943428,  0.2428346,  0.7500223, -1.2440371,  1.4582024))

	expect_error(
		scale_abundance(
			left_join(input_df, sam) |> identify_abundant(a, b, c, factor_of_interest = condition_cont),
			.sample = a,
			.transcript = b,
			.abundance = c
		))

})


test_that("Adding scaled counts - no object",{

	res =
		scale_abundance(
			input_df |> identify_abundant(a, b, c),
			.sample = a,
			.transcript = b,
			.abundance = c,
			action = "add"
		)

	expect_equal(
		unique(res$multiplier),
		c(1.2994008, 1.1781297, 2.6996428, 0.9702628, 1.8290148),
		tolerance=1e-3
	)

	expect_equal(
		ncol(res),
		13
	)

})

test_that("Scaled counts - subset",{

  res =
    input_df |>
    identify_abundant(a, b, c) |>
    scale_abundance(
      .sample = a,
      .transcript = b,
      .abundance = c,
      action = "get",
      .subset_for_scaling = .abundant & grepl("^A", b)
    )

  expect_equal(
    unique(res$multiplier),
    c(1.253886, 1.099169, 1.469270, 1.418187, 1.616244),
    tolerance=1e-3
  )

})

test_that("quantile normalisation",{

  res =
    input_df |>
    quantile_normalise_abundance(
      .sample = a,
      .transcript = b,
      .abundance = c,
      action = "get"
    )

  res |>
    pull(c_scaled) |>
    head() |>
  expect_equal(
    c(1052.8 ,  63.8 ,7229.0 ,   2.9, 2143.6, 9272.8),
    tolerance=1e-3
  )

})


test_that("filter variable - no object",{

	res =
		keep_variable(
			input_df,
			.sample = a,
			.transcript = b,
			.abundance = c,
			top = 5
		)

	expect_equal(
		nrow(res),
		25,
		tolerance=1e-3
	)

	expect_equal(
		as.character(sort(unique(res$b))),
		c("FCN1",  "IGHD",  "IGHM",  "IGKC",  "TCL1A")
	)

})

test_that("Only differential trancript abundance - no object",{

	res =
		test_differential_abundance(
			input_df |> identify_abundant(a, b, c, factor_of_interest = condition),
			~ condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			method = "edgeR_likelihood_ratio",
			action="only"
		)

	expect_equal(
		unique(res$logFC)[1:4],
		c(-12.19303, -11.57989, -12.57969, -11.88829),
		tolerance=1e-3
	)

	expect_equal(
		ncol(res),
		6
	)

	expect_equal(	class(attr(res, "internals")$edgeR)[1], 	"DGEGLM"  )

	# Robust version
	res =
		test_differential_abundance(
			input_df |> identify_abundant(a, b, c, factor_of_interest = condition),
			~ condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			method = "edger_robust_likelihood_ratio",
			action="only"
		)

	expect_equal(
		unique(res$logFC)[1:4],
		c(-12.58107, -12.19281, -11.58286, -11.19910),
		tolerance=0.5
	)

	expect_equal(
		ncol(res),
		6
	)

	expect_equal(	class(attr(res, "internals")$edgeR)[1], 	"DGEGLM"  )


	# Continuous covariate
	sam = distinct(input_df, a)
	sam = mutate(sam, condition_cont = c(-0.4943428,  0.2428346,  0.7500223, -1.2440371,  1.4582024))

	res =
		test_differential_abundance(
			left_join(input_df , sam) |> identify_abundant(a, b, c),
			~ condition_cont,
			.sample = a,
			.transcript = b,
			.abundance = c,
			method = "edgeR_likelihood_ratio",
			action="only"
		)

	expect_equal(
		unique(res$logFC)[1:4],
		c(-3.673399, -3.251067, -3.042633,  2.833111),
		tolerance=1e-3
	)

	expect_equal(
		ncol(res),
		6
	)

	expect_equal(	class(attr(res, "internals")$edgeR)[1], 	"DGEGLM"  )

	# Continuous and discrete
	res =
		test_differential_abundance(
			left_join(input_df , sam) |> identify_abundant(a, b, c),
			~ condition_cont + condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			method = "edgeR_likelihood_ratio",
			action="only"
		)

	expect_equal(
		unique(res$logFC)[1:4],
		c(-2.406553, -2.988076, -4.990209, -4.286571),
		tolerance=1e-3
	)

	expect_equal(
		ncol(res),
		6
	)

	expect_equal(	class(attr(res, "internals")$edgeR)[1], 	"DGEGLM"  )

	# Just one covariate error
	expect_error(
		test_differential_abundance(
			filter(input_df |> identify_abundant(a, b, c, factor_of_interest = condition), condition),
			~ condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			method = "edgeR_likelihood_ratio",
			action="only"
		),
		"Design matrix not of full rank"
	)

	# Change scaling method
	res =
		test_differential_abundance(
			input_df |> identify_abundant(a, b, c, factor_of_interest = condition),
			~ condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			scaling_method = "TMM",
			method = "edgeR_likelihood_ratio",
			action="only"
		)

	# Treat
	input_df |>
		identify_abundant(a, b, c, factor_of_interest = condition) |>
		test_differential_abundance(
			~ condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			scaling_method = "TMM",
			method = "edgeR_likelihood_ratio",
			test_above_log2_fold_change = 1,
			action="get"
		) |>
		filter(FDR<0.05) |>
		nrow() |>
		expect_equal(169)

})

test_that("Only differential trancript abundance - no object - with contrasts",{

	res =
		test_differential_abundance(
			input_df |> identify_abundant(a, b, c, factor_of_interest = condition),
			~ 0 + condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			contrasts = c( "conditionTRUE - conditionFALSE",  "conditionFALSE - conditionTRUE"),
			method = "edgeR_likelihood_ratio",
			action="only"
		)

	expect_equal(
		unique(res$`logFC___conditionTRUE - conditionFALSE`)[1:4],
		c(-12.19303, -11.57989, -12.57969, -11.88829),
		tolerance=1e-3
	)

	expect_equal(
		ncol(res),
		11
	)

	expect_equal(	class(attr(res, "internals")$edgeR)[1], 	"DGEGLM"  )

})

test_that("Add differential trancript abundance - no object",{

	res =
		test_differential_abundance(
			input_df |> identify_abundant(a, b, c, factor_of_interest = condition),
			~ condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			method = "edgeR_likelihood_ratio",
			action="add"
		)

	expect_equal(
		dplyr::pull(dplyr::slice(distinct(res, b, logFC), 1:4) , "logFC"),
		c(3.597633, 2.473975, 2.470380,       NA),
		tolerance=1e-3
	)

	expect_equal(
		ncol(res),
		15
	)

	expect_equal(	class(attr(res, "internals")$edgeR)[1], 	"DGEGLM"  )

})

test_that("Only differential trancript abundance voom - no object",{

	res =
		test_differential_abundance(
			input_df |> identify_abundant(a, b, c, factor_of_interest = condition),
			~ condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			method = "limma_voom",
			action="only"
		)

	expect_equal(
		unique(res$logFC)[1:4],
		c(-12.25012, -11.48490, -10.29393, -11.69070),
		tolerance=1e-3
	)

	expect_equal(
		ncol(res),
		7
	)

	expect_equal(	class(attr(res, "internals")$voom)[1], 	"MArrayLM"  )

	# Continuous covariate
	sam = distinct(input_df, a)
	sam = mutate(sam, condition_cont = c(-0.4943428,  0.2428346,  0.7500223, -1.2440371,  1.4582024))

	res =
		test_differential_abundance(
			left_join(input_df , sam) |> identify_abundant(a, b, c),
			~ condition_cont,
			.sample = a,
			.transcript = b,
			.abundance = c,
			method = "limma_voom",
			action="only"
		)

	expect_equal(
		unique(res$logFC)[1:4],
		c(-3.659127, -3.233295, -3.750860 ,-3.033059),
		tolerance=1e-3
	)

	expect_equal(
		ncol(res),
		7
	)

	expect_equal(	class(attr(res, "internals")$voom)[1], 	"MArrayLM"  )

	# Continuous and discrete
	res =
		test_differential_abundance(
			left_join(input_df , sam) |> identify_abundant(a, b, c),
			~ condition_cont + condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			method = "limma_voom",
			action="only"
		)

	expect_equal(
		unique(res$logFC)[1:4],
		c(-2.981563, -4.883692, -1.702294,  2.423231),
		tolerance=1e-3
	)

	expect_equal(
		ncol(res),
		7
	)

	expect_equal(	class(attr(res, "internals")$voom)[1], 	"MArrayLM"  )

	# Just one covariate error
	expect_error(
		test_differential_abundance(
			filter(input_df |> identify_abundant(a, b, c, factor_of_interest = condition), condition),
			~ condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			method = "limma_voom",
			action="only"
		),
		"Subsetting to non-estimable coefficients is not allowed"
	)

	# Change scaling method
	res =
		test_differential_abundance(
			input_df |> identify_abundant(a, b, c, factor_of_interest = condition),
			~ condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			scaling_method = "TMM",
			method = "limma_voom",
			action="only"
		)
})

test_that("Only differential trancript abundance - no object - with contrasts",{

	res =
		test_differential_abundance(
			input_df |> identify_abundant(a, b, c, factor_of_interest = condition),
			~ 0 + condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			contrasts = c( "conditionTRUE - conditionFALSE",  "conditionFALSE - conditionTRUE"),
			method = "limma_voom",
			action="only"
		)

	expect_equal(
		unique(res$`logFC___conditionTRUE - conditionFALSE`)[1:4],
		c(-12.25012, -11.48490 ,-10.29393, -11.69070),
		tolerance=1e-3
	)

	expect_equal(
		ncol(res),
		13
	)

	expect_equal(	class(attr(res, "internals")$voom)[1], 	"MArrayLM"  )

})

test_that("Voom with sample weights method",{

	res =
		test_differential_abundance(
			input_df |> identify_abundant(a, b, c, factor_of_interest = condition),
			~ condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			method = "limma_voom_sample_weights",
			action="only"
		)

	expect_equal(
		unique(res$logFC)[1:4],
		c(-10.357682, -11.624146, -12.121186, -9.287714),
		tolerance=1e-3
	)

	expect_equal(
		ncol(res),
		7
	)
})

test_that("Voom with treat method",{

	input_df |>
        identify_abundant(a, b, c, factor_of_interest = condition) |>
		test_differential_abundance(
			~ condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			method = "limma_voom",
			test_above_log2_fold_change = 1,
			action="only"
		) |>
	    filter(adj.P.Val<0.05) |>
		nrow() |>
		expect_equal(97)

        # with multiple contrasts
    	res <-
    	input_df |>
    	rename(cell_type = `Cell type`) |>
        identify_abundant(a, b, c, factor_of_interest = cell_type) |>
		test_differential_abundance(
			~ 0 + cell_type,
			.sample = a,
			.transcript = b,
			.abundance = c,
			contrasts = c("cell_typeb_cell-cell_typemonocyte", "cell_typeb_cell-cell_typet_cell"),
			method = "limma_voom",
			test_above_log2_fold_change = 1,
			action="only"
		)

    	res |>
    	filter(`adj.P.Val___cell_typeb_cell-cell_typemonocyte` < 0.05) |>
		nrow() |>
		expect_equal(293)

    	res |>
    	filter(`adj.P.Val___cell_typeb_cell-cell_typet_cell`<0.05) |>
		nrow() |>
		expect_equal(246)

})

test_that("New method choice",{

	res =
		test_differential_abundance(
			input_df |> identify_abundant(a, b, c, factor_of_interest = condition),
			~ condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			method = "edgeR_quasi_likelihood",
			action="only"
		)

	expect_equal(
		unique(res$logFC)[1:4],
		c(-11.583849, -12.192713,  -8.927257,  -7.779931),
		tolerance=1e-3
	)

	expect_equal(
		ncol(res),
		6
	)

	expect_equal(	class(attr(res, "internals")$edgeR)[1], 	"DGEGLM"  )

	# Wrong method
	expect_error(
		test_differential_abundance(
			filter(input_df |> identify_abundant(a, b, c, factor_of_interest = condition), condition),
			~ condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			method = "WRONG_METHOD",
			action="only"
		),
		"the only methods supported"
	)
})

test_that("DESeq2 differential trancript abundance - no object",{

  if (find.package("DESeq2", quiet = TRUE) |> length() |> equals(0)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install("DESeq2", ask = FALSE)
  }

  test_deseq2_df = DESeq2::DESeqDataSet(se_mini,design=~condition)
  colData(test_deseq2_df)$condition = factor(colData(test_deseq2_df)$condition)

  res_deseq2 =
	  test_deseq2_df |>
		DESeq2::DESeq() |>
		DESeq2::results()

	res_tidybulk =
		test_deseq2_df |>
		tidybulk() |>
		test_differential_abundance(~condition, method="DESeq2", action="get")


	expect_equal(
		res_tidybulk |> dplyr::slice(c(1, 3,4, 6)) |> dplyr::pull(log2FoldChange),
		res_deseq2[c(1, 3,4, 6),2],
		tolerance =0.005
	)

	res =
		test_differential_abundance(
			input_df |> identify_abundant(a, b, c, factor_of_interest = condition),
			~ condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			method = "deseq2",
			action="only"
		)

	expect_equal(
		unique(res$log2FoldChange)[1:4],
		c(3.449740, 2.459516, 2.433466, 1.951263),
		tolerance=1e-3
	)

	expect_equal(
		ncol(res),
		7
	)

	expect_equal(	class(attr(res, "internals")$DESeq2)[1], 	"DESeqDataSet"  )

	# Continuous covariate
	sam = distinct(input_df, a)
	sam = mutate(sam, condition_cont = c(-0.4943428,  0.2428346,  0.7500223, -1.2440371,  1.4582024))

	res =
	  left_join(input_df , sam) |> identify_abundant(a, b, c) |>
		test_differential_abundance(
			~ condition_cont,
			.sample = a,
			.transcript = b,
			.abundance = c,
			method = "deseq2",
			action="only"
		)

	expect_equal(
		unique(res$log2FoldChange)[1:4],
		c(-1.1906895, -0.4422231,  0.9656122, -0.3328515),
		tolerance=1e-3
	)

	expect_equal(
		ncol(res),
		7
	)

	expect_equal(	class(attr(res, "internals")$DESeq2)[1], 	"DESeqDataSet"  )

	# Continuous and discrete
	res =
		test_differential_abundance(
			left_join(input_df , sam) |> identify_abundant(a, b, c),
			~ condition_cont + condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			method = "deseq2",
			action="only"
		)

	expect_equal(
		unique(res$log2FoldChange)[1:4],
		c(1.1110210, 3.0077779, 0.9024256, 4.7259792),
		tolerance=1e-3
	)

	expect_equal(
		ncol(res),
		7
	)

	expect_equal(	class(attr(res, "internals")$DESeq2)[1], 	"DESeqDataSet"  )

	# Just one covariate error
	expect_error(
		test_differential_abundance(
			filter(input_df |> identify_abundant(a, b, c, factor_of_interest = condition), condition),
			~ condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			method = "deseq2",
			action="only"
		),
		"design has a single variable"
	)

	# # Contrasts
	# input_df |>
	# 	identify_abundant(a, b, c, factor_of_interest = condition) |>
	# 	test_differential_abundance(
	# 		~ 0 + condition,
	# 		.sample = a,
	# 		.transcript = b,
	# 		.abundance = c,
	# 		method = "deseq2",
	# 		contrasts = "this_is - wrong",
	# 		action="only"
	# 	) |>
	# expect_error("for the moment, the contrasts argument")

	deseq2_contrasts =
		input_df |>
		identify_abundant(a, b, c, factor_of_interest = condition) |>
		test_differential_abundance(
			~ 0 + condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			method = "deseq2",
			contrasts = list(c("condition", "TRUE", "FALSE")),
			action="only"
		)

	edger_contrasts =
		input_df |>
		identify_abundant(a, b, c, factor_of_interest = condition) |>
		test_differential_abundance(
			~ 0 + condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			contrasts = "conditionTRUE - conditionFALSE",
			action="only"
		)


	expect_gt(
		(deseq2_contrasts |> filter(b=="ABCB4") |> pull(3)) *
			(edger_contrasts |> filter(b=="ABCB4") |> pull(2)),
		0
	)

        # DESeq2 with lfcThreshold using test_above_log2_fold_change
        res_lfc <- test_deseq2_df |>
          keep_abundant(factor_of_interest = condition) |>
          test_differential_abundance(~condition, method="DESeq2", test_above_log2_fold_change=1) |>
          mcols()

        # DESeq2::plotMA(DESeq2::DESeqResults(res_lfc))

        # significant are outside of LFC 1 / -1
        idx <- which(res_lfc$padj < .1)
        expect_true(all(abs(res_lfc$log2FoldChange[idx]) > 1))


})

test_that("differential trancript abundance - random effects",{

  my_input =
    input_df |>
    identify_abundant(a, b, c, factor_of_interest = condition) |>
    mutate(time = time |> stringr::str_replace_all(" ", "_")) |>

    filter(b %in% c("ABCB4" , "ABCB9" , "ACAP1",  "ACHE",   "ACP5" ,  "ADAM28"))

  my_input |>
    test_differential_abundance(
      ~ condition + (1 + condition | time),
      .sample = a,
      .transcript = b,
      .abundance = c,
      method = "glmmseq_lme4",
      action="only",
      cores = 1
    ) |>
    pull(P_condition_adjusted) |>
    head(4) |>
    expect_equal(
      c(0.03643414, 0.02938584, 0.02938584, 0.03643414),
      tolerance=1e-3
    )

  # Custom dispersion
  my_input =
    my_input |>
    left_join(
      my_input |> pivot_transcript(b) |> mutate(disp_ = 2 ),
      by = join_by(b, entrez, .abundant)
    )


    my_input |>
    test_differential_abundance(
      ~ condition + (1 + condition | time),
      .sample = a,
      .transcript = b,
      .abundance = c,
      method = "glmmseq_lme4",
      action="only",
      cores = 1,
      .dispersion = disp_
    ) |>
    pull(P_condition_adjusted) |>
    head(4) |>
    expect_equal(
      c(0.1081176, 0.1303558, 0.1303558, 0.1693276),
      tolerance=1e-2
    )

})


test_that("test prefix",{

  library(stringr)

	df = input_df |> tidybulk(a, b, c, ) |> identify_abundant(factor_of_interest = condition)

	res_DeSEQ2 =
		df |>
		test_differential_abundance(~condition, method="DESeq2", action="only", prefix = "prefix_")

	res_voom =
		df |>
		test_differential_abundance(~condition, method="limma_voom", action="only", prefix = "prefix_")

	res_voom_sample_weights =
	    df |>
	    test_differential_abundance(~condition, method="limma_voom_sample_weights", action="only", prefix = "prefix_")

	res_edger =
		df |>
		test_differential_abundance(~condition, method="edgeR_likelihood_ratio", action="only", prefix = "prefix_")

	expect_gt(colnames(res_DeSEQ2) |> str_which("prefix_") |> length(), 0)
	expect_gt(colnames(res_voom) |> str_which("prefix_") |> length(), 0)
	expect_gt(colnames(res_voom_sample_weights) |> str_which("prefix_") |> length(), 0)
	expect_gt(colnames(res_edger) |> str_which("prefix_") |> length(), 0)
})


test_that("Get entrez from symbol - no object",{

	res =
	  input_df |>
	  select(-entrez) |>
	  symbol_to_entrez(.transcript = b, .sample = a)

	expect_equal(
		res$entrez[1:4],
		c( "5244",  "23457", "9744",  "43" )
	)

})

# test_that("Get gene enrichment - no object",{
#
# 	if (find.package("EGSEA", quiet = TRUE) |> length() |> equals(0)) {
# 		message("Installing EGSEA needed for differential transcript abundance analyses")
# 		if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
# 		BiocManager::install("EGSEA")
# 	}
#
# 	library(EGSEA)
#
# 	res =
# 		test_gene_enrichment(
# 			aggregate_duplicates(
# 				dplyr::rename(symbol_to_entrez(
# 					#dplyr::filter(input_df, grepl("^B", b)),
# 					input_df,
# 					.transcript = b, .sample = a), d = entrez
# 				),
# 				.transcript = d,
# 				.sample = a,
# 				.abundance = c
# 			) |> identify_abundant(a, b, c, factor_of_interest = condition),
# 			~ condition,
# 			.sample = a,
# 			.entrez = d,
# 			.abundance = c,
# 			species="human"
# 		)
#
# 	expect_equal(
# 		res$pathway[1:4],
# 		c("GNF2_HCK"    ,  "GSE10325_LUPUS_BCELL_VS_LUPUS_MYELOID_DN"   ,"Amino sugar and nucleotide sugar metabolism", "Phagosome"  )
# 	)
#
# 	expect_equal(
# 		ncol(res),
# 		20
# 	)
#
# })
#

test_that("Only adjusted counts - no object",{

	cm = input_df
	cm$batch = 0
	cm$batch[cm$a %in% c("SRR1740035", "SRR1740043")] = 1

	res =
	  cm |>
	  identify_abundant(a, b, c) |>
		adjust_abundance(
			~ condition + batch,
			.sample = a,
			.transcript = b,
			.abundance = c,
			method = "combat",
			action="only"
		)

	expect_equal(
		unique(res$`c_adjusted`)[c(1, 2, 3, 5)],
		c( 7948 ,2193 , 262, 8152),
		tolerance=1e-3
	)

	expect_equal( ncol(res), 3 )

})

test_that("Get adjusted counts - no object",{

	cm = input_df
	cm$batch = 0
	cm$batch[cm$a %in% c("SRR1740035", "SRR1740043")] = 1

	res =
		adjust_abundance(
			cm |> identify_abundant(a, b, c),
			 ~ condition + batch,
			.sample = a,
			.transcript = b,
			.abundance = c,
			method = "combat",
			action="get"
		)

	expect_equal(
		unique(res$`c_adjusted`)[c(1, 2, 3, 5)],
		c( 7948 ,2193 , 262, 8152),
		tolerance=1e-3
	)

	expect_equal(
		ncol(res),
		9
	)

})

test_that("Add adjusted counts - no object",{

	cm = input_df
	cm$batch = 0
	cm$batch[cm$a %in% c("SRR1740035", "SRR1740043")] = 1

	res =
		adjust_abundance(
			cm |> identify_abundant(a, b, c),
			~ condition + batch,
			.sample = a,
			.transcript = b,
			.abundance = c,
			method = "combat",
			action="add"
		)

	expect_equal(
		unique(res$`c_adjusted`)[c(1, 2, 3, 5)],
		c(NA, 7948, 2193, 2407),
		tolerance=1e-3
	)

	expect_equal(
		ncol(res),
		12
	)

})


test_that("Get adjusted counts multiple factors - SummarizedExperiment",{

  cm = input_df
  cm$batch = 0
  cm$batch[cm$a %in% c("SRR1740035", "SRR1740043")] = 1

  cm =
    cm |>
    identify_abundant(a, b, c)
  cm$c  = as.integer(cm$c )

  res =
    cm |>
    identify_abundant(a, b, c) |>
    adjust_abundance(.factor_unwanted = c(time, batch),
                     .factor_of_interest = dead,
                     .sample = a,
                     .transcript = b,
                     .abundance = c,
                     method = "combat_seq",
                     shrink.disp = TRUE,
                     shrink = TRUE,
                     gene.subset.n = 100
    )

  # Usa MCMC so it is stokastic
  # expect_equal(
  #   unique(res$`c_adjusted`)[c(1, 2, 3, 5)],
  #   c(NA, 5059, 1942, 2702),
  #   tolerance=1e-3
  # )



})

test_that("Only cluster lables based on Kmeans - no object",{

	res =
		cluster_elements(
			input_df,
			.abundance = c,
			.element = a,
			.feature = b,
			method="kmeans",
			centers = 2,
			action="only"
		)

	expect_equal(
		typeof(res$`cluster kmeans`),
		"integer"
	)

	expect_equal(
		ncol(res),
		2
	)

})

test_that("Get cluster lables based on Kmeans - no object",{

	res =
		cluster_elements(
			input_df,
			.abundance = c,
			.element = a,
			.feature = b,
			method="kmeans",
			centers = 2,
			action="get"
		)

	expect_equal(
		typeof(res$`cluster kmeans`),
		"integer"
	)

	expect_equal(
		ncol(res),
		7
	)
	expect_equal(
		nrow(res),
		5
	)

})

test_that("Add cluster lables based on Kmeans - no object",{

	res =
		cluster_elements(
			input_df,
			.abundance = c,
			.element = a,
			.feature = b,
			method="kmeans",
			centers = 2,
			action="add"
		)

	expect_equal(
		typeof(res$`cluster kmeans`),
		"integer"
	)

	expect_equal(
		ncol(res),
		10
	)

})

# test_that("Only cluster lables based on SNN - no object",{
#
# 	res =
# 		cluster_elements(
# 			input_df_breast,
# 			.element = a,
# 			.feature = b,
# 			.abundance = `c norm`,
# 			method="SNN",
# 			resolution=0.8,
# 			action="only"
# 		)
#
# 	expect_equal(
# 		typeof(res$`cluster SNN`),
# 		"integer"
# 	)
#
# 	expect_equal(
# 		ncol(res),
# 		2
# 	)
#
# })

# test_that("Get cluster lables based on SNN - no object",{
#
# 	res =
# 		cluster_elements(
# 			input_df_breast,
# 			.element = a,
# 			.feature = b,
# 			.abundance = `c norm`,
# 			method="SNN",
# 			resolution=0.8,
# 			action="get"
# 		)
#
# 	expect_equal(
# 		typeof(res$`cluster SNN`),
# 		"integer"
# 	)
#
# 	expect_equal(
# 		ncol(res),
# 		3
# 	)
# 	expect_equal(
# 		nrow(res),
# 		251
# 	)
#
# })

# test_that("Add cluster lables based on SNN - no object",{
#
# 	res =
# 		cluster_elements(
# 			input_df_breast,
# 			.element = a,
# 			.feature = b,
# 			.abundance = `c norm`,
# 			method="SNN",
# 			resolution=0.8,
# 			action="add"
# 		)
#
# 	expect_equal(
# 		typeof(res$`cluster SNN`),
# 		"integer"
# 	)
#
# 	expect_equal(
# 		ncol(res),
# 		6
# 	)
#
# })

test_that("Only reduced dimensions MDS - no object",{

	res =
		reduce_dimensions(
			input_df |> identify_abundant(a, b, c),
			method = "MDS",
			.abundance = c,
			.element = a,
			.feature = b,
			action="only"
		)

	expect_equal(
		res$`Dim1`,
		c(1.4048441,  1.3933490, -2.0138120 , 0.8832354, -1.6676164),
		tolerance=10
	)

	expect_equal(
		ncol(res),
		3
	)

	expect_equal(	class(attr(res, "internals")$MDS[[1]])[1], 	"MDS"  )

	inp = input_df |> identify_abundant(a, b, c)

	# Duplicate genes/samples
	expect_error(
		reduce_dimensions(
		  inp |> bind_rows( inp|> dplyr::slice(1) |> mutate(c = c+1) ),
			method = "MDS",
			.abundance = c,
			.element = a,
			.feature = b, action="only"
		),
		"Your dataset include duplicated "
	)
})

test_that("Get reduced dimensions MDS - no object",{

	res =
		reduce_dimensions(
			input_df |> identify_abundant(a, b, c),
			method = "MDS",
			.abundance = c,
			.element = a,
			.feature = b,
			action="get"
		)

	expect_equal(
		(res$`Dim1`)[1:4],
		c( -0.8794274, -0.8976436 , 1.4564831 ,-1.0074328),
		tolerance=10
	)

	expect_equal(
		ncol(res),
		8
	)
	expect_equal(
		nrow(res),
		5
	)
	expect_equal(	class(attr(res, "internals")$MDS[[1]])[1], 	"MDS"  )
})

test_that("Add reduced dimensions MDS - no object",{

	res =
		reduce_dimensions(
			input_df |> identify_abundant(a, b, c),
			method = "MDS",
			.abundance = c,
			.element = a,
			.feature = b,
			action="add"
		)

	expect_equal(
		(res$`Dim1`)[1:4],
		c( 1.404844, 1.404844, 1.404844, 1.404844),
		tolerance=10
	)

	expect_equal(
		ncol(res),
		12
	)

	expect_equal(	class(attr(res, "internals")$MDS[[1]])[1], 	"MDS"  )
})

test_that("Only reduced dimensions PCA - no object",{

	res =
		reduce_dimensions(
			input_df |> identify_abundant(a, b, c),
			method="PCA",
			.abundance = c,
			.element = a,
			.feature = b,
			action="only"
		)

	expect_equal(
		res$PC1,
		c( -7.214337, -7.563078, 11.638986, -8.069080 ,11.207509),
		tolerance=1e-1
	)

	expect_equal(
		ncol(res),
		3
	)

	expect_equal(	class(attr(res, "internals")$PCA), 	"prcomp"  )

	# Duplicate genes/samples
	inp = input_df |> identify_abundant(a, b, c)

	expect_error(
		reduce_dimensions(
		  inp |> bind_rows( inp |> dplyr::slice(1) |> mutate(c = c+1) ),
			method = "PCA",
			.abundance = c,
			.element = a,
			.feature = b, action="only"
		),
		"include duplicated sample/gene pairs"
	)
})

test_that("Get reduced dimensions PCA - no object",{

	res =
		reduce_dimensions(
			input_df |> identify_abundant(a, b, c),
			method="PCA",
			.abundance = c,
			.element = a,
			.feature = b,
			action="get"
		)

	expect_equal(
		typeof(res$`PC1`),
		"double"
	)

	expect_equal(
		ncol(res),
		8
	)

	expect_equal(	class(attr(res, "internals")$PCA), 	"prcomp"  )

})


test_that("Add reduced dimensions PCA - no object",{

	res =
		reduce_dimensions(
			input_df |> identify_abundant(a, b, c),
			method="PCA",
			.abundance = c,
			.element = a,
			.feature = b,
			action="add"
		)

	res |>
	  pull(PC1) |>
	  magrittr::extract2(1) |>
	  expect_equal(-7.214337, tolerance = 0.01)

	expect_equal(
		typeof(res$`PC1`),
		"double"
	)

	expect_equal(
		ncol(res),
		12
	)

	expect_equal(	class(attr(res, "internals")$PCA), 	"prcomp"  )

})

test_that("Get reduced dimensions tSNE - no object",{

	set.seed(132)

	res =
		input_df_breast |>
		identify_abundant(a, b, c) |>

		reduce_dimensions(
			method="tSNE",
			.abundance = c,
			.element = a,
			.feature = b,
			action="get",
			verbose=FALSE
		) |>
	  suppressMessages()

	# DOES NOT REPRODUCE MAYBE BECAUSE OF DIFFERENT TSNE VERSIONS
	# res |>
	#   pull(tSNE1) |>
	#   magrittr::extract2(1) |>
	#   expect_equal(2.432608, tolerance = 0.01)

	expect_equal(
		typeof(res$`tSNE1`),
		"double",
		tolerance=1e-1
	)

	expect_equal(
		ncol(res),
		4
	)
	expect_equal(
		nrow(res),
		251
	)

	# Duplicate genes/samples
	inp = input_df |> identify_abundant(a, b, c)

	expect_error(
		reduce_dimensions(
		  inp |> bind_rows( inp |> dplyr::slice(1) |> mutate(c = c+1) ),
			method = "tSNE",
			.abundance = c,
			.element = a,
			.feature = b, action="get",
			verbose=FALSE
		),
		"include duplicated sample/gene pairs"
	)

})

test_that("Add reduced dimensions UMAP - no object",{

  set.seed(132)

  res =
    input_df_breast |>
    identify_abundant(a, b, c) |>

    reduce_dimensions(
      method="UMAP",
      .abundance = c,
      .element = a,
      .feature = b,
      action="add"
    )

  res |>
    pull(UMAP1) |>
    magrittr::extract2(1) |>
    expect_equal(-2.12, tolerance = 0.01)

  expect_equal(ncol(res), 8)

})

test_that("Only rotated dimensions - no object",{

	res.pca =
		reduce_dimensions(
			input_df |> identify_abundant(a, b, c),
			method="PCA",
			.abundance = c,
			.element = a,
			.feature = b,
			action="add"
		)

	res =
		rotate_dimensions(
			res.pca,
			dimension_1_column = PC1,
			dimension_2_column = PC2,
			rotation_degrees = 45,
			.element = a,
			action="only"
		)

	expect_equal(
		res$`PC1_rotated_45`,
		c(-9.450807, -9.739338,  8.853659,  2.741059,  7.595427),
		tolerance=1e-1
	)

	expect_equal(
		ncol(res),
		3
	)

})

test_that("Get rotated dimensions - no object",{

	res.pca =
		reduce_dimensions(
			input_df |> identify_abundant(a, b, c),
			method="PCA",
			.abundance = c,
			.element = a,
			.feature = b,
			action="get"
		)

	res =
		rotate_dimensions(
			res.pca,
			dimension_1_column = PC1,
			dimension_2_column = PC2,
			rotation_degrees = 45,
			.element = a,
			action="get"
		)

	expect_equal(
		res$`PC1_rotated_45`[1:4],
		c(  -9.450807, -9.739338 , 8.853659,  2.741059),
		tolerance=1e-1
	)

	expect_equal(
		ncol(res),
		10
	)
	expect_equal(
		nrow(res),
		5
	)
})


test_that("Add rotated dimensions - no object",{

	res.pca =
		reduce_dimensions(
			input_df |> identify_abundant(a, b, c),
			method="PCA",
			.abundance = c,
			.element = a,
			.feature = b,
			action="add"
		)

	res =
		rotate_dimensions(
			res.pca,
			dimension_1_column = PC1,
			dimension_2_column = PC2,
			rotation_degrees = 45,
			.element = a,
			action="add"
		)

	expect_equal(
		res$`PC1_rotated_45`[1:4],
		c( -9.450807 ,-9.450807, -9.450807 ,-9.450807),
		tolerance=1e-1
	)

	expect_equal(
		ncol(res),
		14
	)

})

test_that("Aggregate duplicated transcript - no object",{

	res =
		aggregate_duplicates(
			input_df,
			.sample = a,
			.transcript = b,
			.abundance = c
		)

	expect_equal(
		res$b[1:4],
		c("ABCB4", "ABCB9", "ACAP1", "ACHE"  )
	)

	expect_equal(
		ncol(res),
		10
	)

})

test_that("Drop redundant correlated - no object",{

	res =
		remove_redundancy(
			input_df,
			method = "correlation",
			.abundance = c,
			.element = a,
			.feature = b
		)

	expect_equal(
		nrow(res),
		2108
	)

	expect_equal(
		ncol(res),
		9
	)

})


test_that("Only symbol from ensambl - no object",{

	# Human
	res =
		ensembl_to_symbol(
			tidybulk::counts_ensembl,
			.ensembl = ens,
			action="only"
		)

	expect_equal(
		as.character(res$transcript),
		"TSPAN6"
	)

	expect_equal(
		ncol(res),
		3
	)

	# Mouse
	# Human
	res =
		ensembl_to_symbol(
			tibble(ens = c("ENSMUSG00000000001",
												 "ENSMUSG00000000003",
												 "ENSMUSG00000000028",
												 "ENSMUSG00000000031",
												 "ENSMUSG00000000037",
												 "ENSMUSG00000000049"
			)),
			.ensembl = ens,
			action="only"
		)

	expect_equal(
		as.character(res$transcript)[1],
		"Gnai3"
	)

	expect_equal(
		ncol(res),
		3
	)
})

test_that("Add description to symbol",{

	# Human
	res =
		describe_transcript(
			input_df,
			.transcript = b
		)


	expect_equal(
		ncol(res),
		10
	)

	res =
		describe_transcript(
			input_df |> tidybulk(a, b, c)
		)


	expect_equal(
		ncol(res),
		10
	)

})

test_that("Add symbol from ensambl - no object",{

	res =
		ensembl_to_symbol(
			tidybulk::counts_ensembl,
			.ensembl = ens,
			action="add"
		)

	expect_equal(
		res$`read count`[1:4],
		c(144,   72,    0 ,1099)
	)

	expect_equal(
		ncol(res),
		8
	)

})

test_that("Only cell type proportions - no object",{

	# Cibersort
	res =
		deconvolve_cellularity(
			input_df,
			.sample = a,
			.transcript = b,
			.abundance = c,
			action="only", cores=1
		)

	expect_equal(
		as.numeric(res[1,2:5]),
		c(0.6196622, 0.2525598, 0.0000000, 0.0000000),
		tolerance=1e-3
	)

	expect_equal(
		ncol(res),
		23
	)

	# LLSR
	res =
		deconvolve_cellularity(
			input_df,
			.sample = a,
			.transcript = b,
			.abundance = c,
			method = "llsr",
			action="only"
		)

	expect_equal(
		as.numeric(res[1,2:5]),
		c(0.6702025807, 0.0000000000, 0.0000000000, 0.0005272016),
		tolerance=1e-3
	)

	expect_equal(
		ncol(res),
		23
	)


	# # EPIC
	# library(EPIC)
	# res =
	# 	deconvolve_cellularity(
	# 		input_df,
	# 		.sample = a,
	# 		.transcript = b,
	# 		.abundance = c,
	# 		method = "epic",
	# 		action="only", cores=1
	# 	)
	#
	# expect_equal(
	# 	as.numeric(res[1,2:5]),
	# 	c(7.750225e-01, 2.199743e-01, 5.808702e-05, 1.944729e-06),
	# 	tolerance=1e-3
	# )
	#
	# expect_equal(
	# 	ncol(res),
	# 	24
	# )

})

test_that("differential composition",{

	# Cibersort
	test_differential_cellularity(
			input_df,
			. ~ condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			cores = 1
		) |>
		pull(`estimate_(Intercept)`) |>
    magrittr::extract2(1) |>
		as.integer() |>
		expect_equal(	-2, 	tollerance =1e-3)

	# llsr
	test_differential_cellularity(
		input_df,
		. ~ condition,
		.sample = a,
		.transcript = b,
		.abundance = c,
		method="llsr"
	) |>
		pull(`estimate_(Intercept)`) |>
		magrittr::extract2(1) |>
		as.integer() |>
		expect_equal(	-2, 	tollerance =1e-3)

	# Survival analyses
	input_df |>
	select(a, b, c) |>
	nest(data = -a) |>
	mutate(
		days = c(1, 10, 500, 1000, 2000),
		dead = c(1, 1, 1, 0, 1)
	) |>
	unnest(data) |>
	test_differential_cellularity(
		survival::Surv(days, dead) ~ .,
		.sample = a,
		.transcript = b,
		.abundance = c,
		cores = 1
	) |>
	pull(estimate) |>
	magrittr::extract2(1) |>
	expect_equal(26.2662279, tolerance = 30)
	# round() %in% c(
	# 	26,  # 97 is the github action MacOS that has different value
	# 	26, # 112 is the github action UBUNTU that has different value
	# 	26 # 93 is the github action Windows that has different value
	# ) |>
	# expect_true()

})

test_that("test_stratification_cellularity",{

	# Cibersort
	input_df |>
		select(a, b, c) |>
		nest(data = -a) |>
		mutate(
			days = c(1, 10, 500, 1000, 2000),
			dead = c(1, 1, 1, 0, 1)
		) |>
		unnest(data) |>
		test_stratification_cellularity(
			survival::Surv(days, dead) ~ .,
			.sample = a,
			.transcript = b,
			.abundance = c,
			cores = 1
		) |>
		pull(.low_cellularity_expected) |>
		magrittr::extract2(1) |>
		expect_equal(3.35, tolerance  =1e-1)

	# llsr
	input_df |>
		select(a, b, c) |>
		nest(data = -a) |>
		mutate(
			days = c(1, 10, 500, 1000, 2000),
			dead = c(1, 1, 1, 0, 1)
		) |>
		unnest(data) |>
		test_stratification_cellularity(
			survival::Surv(days, dead) ~ .,
			.sample = a,
			.transcript = b,
			.abundance = c,
			method = "llsr"
		) |>
		pull(.low_cellularity_expected) |>
		magrittr::extract2(1) |>
		expect_equal(3.35, tolerance  =1e-1)
})

test_that("filter abundant - no object",{

	res1 =
		identify_abundant(
			input_df,
			.sample = a,
			.transcript = b,
			.abundance = c
		)  |>
	  filter(.abundant) |>
	  nrow()

	res2 =
		identify_abundant(
			input_df,
			.sample = a,
			.transcript = b,
			.abundance = c,
			minimum_proportion = 0.5,
			minimum_counts = 30
		) |>
	  filter(.abundant) |>
	  nrow()

  expect_equal(res1, 910)
  expect_equal(res2, 625)

	expect_gt(res1 ,res2 	)


		keep_abundant(
			input_df,
			.sample = a,
			.transcript = b,
			.abundance = c
		) |>
	  nrow() |>

	expect_equal(910	)

})

test_that("filter abundant with design - no object",{

		identify_abundant(
			input_df,
			.sample = a,
			.transcript = b,
			.abundance = c,
			factor_of_interest = condition
		) |>
    filter(.abundant) |>
    nrow() |>
	expect_equal(1965)



})

test_that("nest - no object",{

	expect_equal(	class(nest(tidybulk(input_df, a, b, c), data = a))[1],	"nested_tidybulk"	)

})

test_that("pivot",{

	expect_equal(	ncol(pivot_sample(tidybulk(input_df, a, b, c))),	6	)

	expect_equal(	ncol(pivot_sample(input_df, a)),	6	)

	expect_equal(	ncol(pivot_transcript(tidybulk(input_df, a, b, c))),	2	)

	expect_equal(	ncol(pivot_transcript(input_df, b)),	2	)

})

# test_that("impute missing - no object",{
#
# 	res =
# 		impute_missing_abundance(
# 			dplyr::slice(input_df, -1),
# 			~ condition,
# 			.sample = a,
# 			.transcript = b,
# 			.abundance = c
# 		)
#
# 	expect_equal(	dplyr::pull(filter(res, b=="TNFRSF4" & a == "SRR1740034"), c),	6	)
#
# 	expect_equal(	ncol(res),	ncol(input_df)	)
#
# 	expect_equal(	nrow(res),	nrow(input_df)	)
#
# })

test_that("gene over representation",{

	df_entrez =  se_mini |> tidybulk() |> as_tibble()
	df_entrez = aggregate_duplicates(df_entrez, aggregation_function = sum, .sample = .sample, .transcript = entrez, .abundance = count)
	df_entrez = mutate(df_entrez, do_test = .feature %in% c("TNFRSF4", "PLCH2", "PADI4", "PAX7"))

	res =
		test_gene_overrepresentation(
			df_entrez,
			.sample = .sample,
			.entrez = entrez,
			.do_test = do_test,
			species="Homo sapiens"
		)

	res |>
	  pull(pvalue) |>
	  magrittr::extract2(1) |>
	  expect_equal(0.0004572092, tolerance = 0.0001	)



})

test_that("as_SummarizedExperiment",{
  input_df |>
    as_SummarizedExperiment(
      .sample = c(a, condition),
      .transcript = c(b, entrez),
      .abundance = c
    ) |>
    nrow() |>
    expect_equal(527)

})

# test_that("bibliography",{
#
# 		tidybulk(
# 			input_df,
# 			.sample = a,
# 			.transcript = b,
# 			.abundance = c
# 		) |>
# 		scale_abundance() |>
# 		get_bibliography()
#
# })
#

