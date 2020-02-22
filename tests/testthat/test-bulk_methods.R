context('Bulk methods')

input_df = setNames(tidyBulk::counts_mini, c("a", "b", "Cell type", "c",  "time" , "condition"))

input_df_breast = setNames(tidyBulk::breast_tcga_mini, c("a", "b", "c norm", "call", "c"))

test_that("Creating tt object from tibble, number of parameters, methods",{

	expect_equal(

		length(
			attr(
				tidyBulk(
					input_df,
					.sample = a,
					.transcript = b,
					.abundance = c
				) ,
				"tt_internals"
			)$tt_columns
		),
		3
	)

})

test_that("Test class identity of tt object",{

	expect_equal(
		class(
			tidyBulk(
				input_df,
				.sample = a,
				.transcript = b,
				.abundance = c
			)
		)[1],
		"tidyBulk"
	)

})

test_that("Getting scaled counts - no object",{

	res =
		scale_abundance(
			input_df,
			.sample = a,
			.transcript = b,
			.abundance = c,
			action = "get"
		)

	expect_equal(
		unique(res$multiplier),
		c(1.124713, 1.011405, 1.511470, 0.828865, 1.191472),
		tolerance=1e-6
	)

	expect_equal(
		ncol(res),
		6
	)

	internals = attr(scale_abundance(tidyBulk(input_df, a, b, c)), "tt_internals")

	expect_equal(length(internals$tt_columns), 4 )

	expect_equal(quo_name(internals$tt_columns[[4]]), "c_scaled" )

})


test_that("Adding scaled counts - no object",{

	res =
		scale_abundance(
			input_df,
			.sample = a,
			.transcript = b,
			.abundance = c,
			action = "add"
		)

	expect_equal(
		unique(res$multiplier),
		c(1.124713 ,1.011405, 1.511470, 0.828865, 1.191472),
		tolerance=1e-6
	)

	expect_equal(
		ncol(res),
		10
	)

})

test_that("filter variable - no object",{

	res =
		filter_variable(
			input_df,
			.sample = a,
			.transcript = b,
			.abundance = c,
			top = 5
		)

	expect_equal(
		nrow(res),
		25,
		tolerance=1e-6
	)

	expect_equal(
		sort(unique(res$b)),
		c("FCN1",  "IGHD",  "IGHM",  "IGKC",  "TCL1A")
	)

})

test_that("Get differential trancript abundance - no object",{

	res =
		test_differential_abundance(
			input_df,
			~ condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			action="get"
		)

	expect_equal(
		unique(res$logFC)[1:4],
		c(-12.48201, -12.10269 ,-11.48896, -13.44406),
		tolerance=1e-6
	)

	expect_equal(
		ncol(res),
		8
	)

	expect_equal(	class(attr(res, "tt_internals")$edgeR)[1], 	"DGEGLM"  )

	# Continuous covariate
	sam = distinct(input_df, a)
	sam = mutate(sam, condition_cont = c(-0.4943428,  0.2428346,  0.7500223, -1.2440371,  1.4582024))

	res =
		test_differential_abundance(
			left_join(input_df, sam),
			~ condition_cont,
			.sample = a,
			.transcript = b,
			.abundance = c,
			action="get"
		)

	expect_equal(
		unique(res$logFC)[1:4],
		c(-4.403626, -4.390714, -4.746695, -4.467572),
		tolerance=1e-6
	)

	expect_equal(
		ncol(res),
		8
	)

	expect_equal(	class(attr(res, "tt_internals")$edgeR)[1], 	"DGEGLM"  )

	# Just one covariate error
	expect_error(
		test_differential_abundance(
			filter(input_df, condition),
			~ condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			action="get"
		),
		"Design matrix not of full rank"
	)

	# Just one sample per covariate error
	expect_error(
		test_differential_abundance(
			filter(input_df, a %in% c("SRR1740034", "SRR1740035", "SRR1740043", "SRR1740058")),
			~ condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			action="get"
		),
		"You need at least two replicated"
	)

	# Setting filtering manually
	res =
		test_differential_abundance(
			input_df,
			~ condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			action="get",
			prop_threshold = 0.5,
			cpm_threshold = 30
		)

	expect_equal(
		unique(res$logFC)[1:4],
		c(-12.33879, -12.66131, -11.70319, -13.62226),
		tolerance=1e-6
	)

	expect_equal(	class(attr(res, "tt_internals")$edgeR)[1], 	"DGEGLM"  )

	# Change scaling method
	res =
		test_differential_abundance(
			input_df,
			~ condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			scaling_method = "TMM",
			action="get"
		)
})

test_that("Get differential trancript abundance - no object - with contrasts",{

	res =
		test_differential_abundance(
			input_df,
			~ 0 + condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			.contrasts = c( "conditionTRUE - conditionFALSE",  "conditionFALSE - conditionTRUE"),
			action="get"
		)

	expect_equal(
		unique(res$`logFC_conditionTRUE - conditionFALSE`)[1:4],
		c(-12.48201, -12.10269, -11.48896, -13.44406),
		tolerance=1e-6
	)

	expect_equal(
		ncol(res),
		14
	)

	expect_equal(	class(attr(res, "tt_internals")$edgeR)[1], 	"DGEGLM"  )

})


test_that("Add differential trancript abundance - no object",{

	res =
		test_differential_abundance(
			input_df,
			~ condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			action="add"
		)

	expect_equal(
		pull(filter(distinct(res, b, logFC), b %in% c("HK3", "FCN1", "CLEC7A", "FAM198B")) , "logFC"),
		c(-12.10269, -12.48201 ,-11.48896, -13.44406),
		tolerance=1e-6
	)

	expect_equal(
		ncol(res),
		13
	)

	expect_equal(	class(attr(res, "tt_internals")$edgeR)[1], 	"DGEGLM"  )

})

test_that("Get entrez from symbol - no object",{

	res =
		symbol_to_entrez(input_df, .transcript = b, .sample = a)

	expect_equal(
		res$entrez[1:4],
		c( "7293",  "9651",  "23569" ,"5081" )
	)

})

# test_that("Get gene enrichment - no object",{
#
# 	if ("EGSEA" %in% rownames(installed.packages()) == FALSE) {
# 		writeLines("Installing EGSEA needed for differential transcript abundance analyses")
# 		if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
# 		BiocManager::install("EGSEA")
# 	}
#
# 	library(EGSEA)
#
# 	res =
# 		analise_gene_enrichment(
# 			aggregate_duplicates(
# 				dplyr::rename(symbol_to_entrez(
# 					#dplyr::filter(input_df, grepl("^B", b)),
# 					input_df,
# 					.transcript = b, .sample = a), d = entrez
# 				),
# 				.transcript = d,
# 				.sample = a,
# 				.abundance = c
# 			),
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

test_that("Get adjusted counts - no object",{

	cm = input_df
	cm$batch = 0
	cm$batch[cm$a %in% c("SRR1740035", "SRR1740043")] = 1

	res =
		adjust_abundance(
			cm,
			~ condition + batch,
			.sample = a,
			.transcript = b,
			.abundance = c,
			action="get"
		)

	expect_equal(
		unique(res$`c_adjusted`)[c(1, 2, 3, 5)],
		c( 693,   47, 7941, 7068),
		tolerance=1e-6
	)

	expect_equal(
		ncol(res),
		3
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
		2
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
		7
	)

})

test_that("Get cluster lables based on SNN - no object",{

	res =
		cluster_elements(
			input_df_breast,
			.element = a,
			.feature = b,
			.abundance = `c norm`,
			method="SNN",
			resolution=0.8,
			action="get"
		)

	expect_equal(
		typeof(res$`cluster SNN`),
		"integer"
	)

	expect_equal(
		ncol(res),
		2
	)

})

test_that("Add cluster lables based on SNN - no object",{

	res =
		cluster_elements(
			input_df_breast,
			.element = a,
			.feature = b,
			.abundance = `c norm`,
			method="SNN",
			resolution=0.8,
			action="add"
		)

	expect_equal(
		typeof(res$`cluster SNN`),
		"integer"
	)

	expect_equal(
		ncol(res),
		6
	)

})

test_that("Get reduced dimensions MDS - no object",{

	res =
		reduce_dimensions(
			input_df,
			method = "MDS",
			.abundance = c,
			.element = a,
			.feature = b,
			action="get"
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

	expect_equal(	class(attr(res, "tt_internals")$MDS)[1], 	"MDS"  )

})

test_that("Add reduced dimensions MDS - no object",{

	res =
		reduce_dimensions(
			input_df,
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
		8
	)

	expect_equal(	class(attr(res, "tt_internals")$MDS)[1], 	"MDS"  )
})

test_that("Get reduced dimensions PCA - no object",{

	res =
		reduce_dimensions(
			input_df,
			method="PCA",
			.abundance = c,
			.element = a,
			.feature = b,
			action="get"
		)

	expect_equal(
		res$PC1,
		c( -0.5070361 ,-0.5086096, -0.4070171, -0.3378837, -0.4521088),
		tolerance=1e-1
	)

	expect_equal(
		ncol(res),
		3
	)

	expect_equal(	class(attr(res, "tt_internals")$PCA), 	"prcomp"  )
})

test_that("Add reduced dimensions PCA - no object",{

	res =
		reduce_dimensions(
			input_df,
			method="PCA",
			.abundance = c,
			.element = a,
			.feature = b,
			action="add"
		)

	expect_equal(
		typeof(res$`PC1`),
		"double"
	)

	expect_equal(
		ncol(res),
		8
	)

	expect_equal(	class(attr(res, "tt_internals")$PCA), 	"prcomp"  )

})

test_that("Add reduced dimensions tSNE - no object",{

	set.seed(132)

	res =
		reduce_dimensions(
			setNames(tidyBulk::counts, c("a", "b", "Cell type", "c",  "time" , "condition", "batch", "factor_of_interest")) ,
			method="tSNE",
			.abundance = c,
			.element = a,
			.feature = b,
			action="add"
		)

	expect_equal(
		typeof(res$`tSNE1`),
		"double",
		tolerance=1e-1
	)

	expect_equal(
		ncol(res),
		10
	)

})

test_that("Get rotated dimensions - no object",{

	res.pca =
		reduce_dimensions(
			input_df,
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
			action="get"
		)

	expect_equal(
		res$`PC1 rotated 45`,
		c(-0.08299217 ,-0.08765521 ,-0.71713866 ,-0.03872173 ,-0.68530405),
		tolerance=1e-1
	)

	expect_equal(
		ncol(res),
		3
	)

})

test_that("Add rotated dimensions - no object",{

	res.pca =
		reduce_dimensions(
			input_df,
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
		res$`PC1 rotated 45`[1:4],
		c( -0.08299217, -0.08299217, -0.08299217, -0.08299217),
		tolerance=1e-1
	)

	expect_equal(
		ncol(res),
		10
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
		c( "TNFRSF4", "PLCH2" ,  "PADI4" ,  "PAX7"   )
	)

	expect_equal(
		ncol(res),
		7
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
		6
	)

})


test_that("Get symbol from ensambl - no object",{

	res =
		annotate_symbol(
			tidyBulk::counts_ensembl,
			.ensembl = ens,
			action="get"
		)

	expect_equal(
		as.character(res$transcript),
		"TSPAN6"
	)

	expect_equal(
		ncol(res),
		3
	)

})

test_that("Add symbol from ensambl - no object",{

	res =
		annotate_symbol(
			tidyBulk::counts_ensembl,
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

test_that("Get cell type proportions - no object",{

	# Cibersort
	res =
		deconvolve_cellularity(
			input_df,
			.sample = a,
			.transcript = b,
			.abundance = c,
			action="get", cores=1
		)

	expect_equal(
		as.numeric(res[1,2:5]),
		c(0.6223514, 0.2378625, 0.0000000 ,0.0000000),
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
			action="get", cores=1
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

})

test_that("Add cell type proportions - no object",{

	res =
		deconvolve_cellularity(
			input_df,
			.sample = a,
			.transcript = b,
			.abundance = c,
			action="add", cores=1
		)

	expect_equal(
		as.numeric(res[1,7:10]),
		c(0.6223514, 0.2378625, 0.0000000 ,0.0000000),
		tolerance=1e-3
	)

	expect_equal(
		ncol(res),
		28
	)

})

test_that("filter abundant - no object",{

	res =
		filter_abundant(
			input_df,
			.sample = a,
			.transcript = b,
			.abundance = c
		)

	expect_equal(
		res$b[1:4],
		c("TNFRSF4", "PLCH2" ,  "PADI4" ,  "CDA"    )
	)

	expect_equal(	ncol(res),	6	)

	expect_equal(	nrow(res),	2220	)

})

test_that("nest - no object",{

	expect_equal(	class(nest(tidyBulk(input_df, a, b, c), data = a))[1],	"tbl_df"	)

})
