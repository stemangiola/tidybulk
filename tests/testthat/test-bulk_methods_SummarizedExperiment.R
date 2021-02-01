context('Bulk methods SummarizedExperiment')

input_df = setNames(tidybulk::counts_mini, c("a", "b", "Cell type", "c",  "time" , "condition"))

input_df_breast = setNames(tidybulk::breast_tcga_mini, c("a", "b", "c norm", "call", "c"))

test_that("tidybulk SummarizedExperiment conversion",{

	res = tidybulk(tidybulk::se)

	expect_equal(	class(res)[1],	"tidybulk"	)

	expect_equal(	nrow(res),	800	)

	expect_equal(	ncol(res),	12	)

	res = res %>% tidybulk:::tidybulk_to_SummarizedExperiment()

	expect_equal(	class(res)[1],	"SummarizedExperiment"	)

	expect_equal(	nrow(res),	100	)

	expect_equal(	ncol(res),	8	)

})

test_that("tidybulk SummarizedExperiment normalisation manual",{

	res = tidybulk(tidybulk:::tidybulk_to_SummarizedExperiment(scale_abundance(tidybulk(se) %>% identify_abundant())))

	expect_equal(
		res[1:4,]$`counts_scaled`,
		c(1796.091258, 1162.818960, 1216.046589,    2.729622),
		tolerance=1e-3
	)

	expect_equal(	nrow(res),	800	)

	expect_equal(	ncol(res),	16	)


	res = rlang::quo_name(attr(res, "internals")$tt_columns[[4]])

	expect_equal( res,	"counts_scaled"	)

})

test_that("tidybulk SummarizedExperiment normalisation",{

	res = scale_abundance(se %>% identify_abundant())

	expect_equal(
		names(SummarizedExperiment::assays(res)),
		c("counts" ,"counts_scaled")
	)

})

test_that("tidybulk SummarizedExperiment clustering",{

	res = cluster_elements(se, method="kmeans", centers = 2)

	expect_equal(
		tail(names(SummarizedExperiment::colData(res)), 1),
		"cluster_kmeans"
	)

	expect_equal(
		levels(SummarizedExperiment::colData(res)$cluster_kmeans),
		c("1", "2")
	)

})

test_that("tidybulk SummarizedExperiment clustering",{

	res = se %>% identify_abundant() %>% reduce_dimensions(method="PCA")

	expect_equal(
		tail(names(SummarizedExperiment::colData(res)), 1),
		"PC2"
	)

})

test_that("Get rotated dimensions - SummarizedExperiment",{

	res.pca =
		reduce_dimensions(se %>% identify_abundant(),		method="PCA"	)

	res =
		rotate_dimensions(
			res.pca,
			dimension_1_column = PC1,
			dimension_2_column = PC2,
			rotation_degrees = 45
		)

	expect_equal(
		tail(names(SummarizedExperiment::colData(res)), 1),
		"PC2.rotated.45"
	)

})

test_that("Drop redundant correlated - SummarizedExperiment",{

	res =
		remove_redundancy(
			se,
			method = "correlation", correlation_threshold = 0.99	)

	expect_equal(
		nrow(res),
		100
	)

})

test_that("Get adjusted counts - SummarizedExperiment",{

	cm = input_df
	cm$batch = 0
	cm$batch[cm$a %in% c("SRR1740035", "SRR1740043")] = 1

	res =
		adjust_abundance(
			tidybulk:::tidybulk_to_SummarizedExperiment(cm, a, b, c) %>% identify_abundant(),
			~ condition + batch
		)

	expect_equal(nrow(res),	527	)

	expect_equal(	names(SummarizedExperiment::assays(res)),	c("c" ,"c_adjusted")	)


})

test_that("Aggregate duplicated transcript - SummarizedExperiment",{

	res =	aggregate_duplicates(	se)

	expect_equal(	dim(res),	c( 100,   8  )	)


})

test_that("Add cell type proportions - SummarizedExperiment",{

	res =		deconvolve_cellularity(tidybulk:::tidybulk_to_SummarizedExperiment(input_df, a, b, c), cores=1	)

	expect_equal(
		as.numeric(as.data.frame(res@colData[1, 4:7])),
		c( 0.6221895, 0.2380869, 0.0000000, 0.0000000),
		tolerance=1e-3
	)

})

test_that("differential trancript abundance - SummarizedExperiment",{

	res =		test_differential_abundance(
		tidybulk:::tidybulk_to_SummarizedExperiment(input_df, a, b, c) %>% 
			identify_abundant(factor_of_interest = condition),	
		~ condition	
	)
	
	w = match(  c("CLEC7A" , "FAM198B", "FCN1"  ,  "HK3"   ), rownames(res) )
	
	# Quasi likelihood
	res_tibble =		test_differential_abundance(
		input_df %>% identify_abundant(a, b, c, factor_of_interest = condition),
		~ condition	,
		a, b, c
	)

	expect_equal(
		res@elementMetadata[w,]$logFC,
		c(-11.58385, -13.53406, -12.58204, -12.19271),
		tolerance=1e-4
	)
	
	expect_equal(
		res@elementMetadata[w,]$logFC,
		res_tibble %>% 
			pivot_transcript(b) %>% 
			filter(b %in% rownames(res)[w]) %>% 
			dplyr::arrange(b) %>%
			dplyr::pull(logFC),
		tolerance=1e-4
	)
	
	# Likelihood ratio
	res2 =		test_differential_abundance(
		tidybulk:::tidybulk_to_SummarizedExperiment(input_df, a, b, c) %>%
			identify_abundant(factor_of_interest = condition),	
		~ condition, method = "edgeR_likelihood_ratio"	)
	
	res2_tibble =		test_differential_abundance(
		input_df %>% identify_abundant(a, b, c, factor_of_interest = condition),
		~ condition	,
		a, b, c, method = "edgeR_likelihood_ratio"
	)

	expect_equal(
		res2@elementMetadata[w,]$logFC,
		c(-11.57989, -13.53476, -12.57969, -12.19303),
		tolerance=1e-4
	)

	expect_equal(
		res2@elementMetadata[w,]$logFC,
		res2_tibble %>% 
			pivot_transcript(b) %>% 
			filter(b %in% rownames(res)[w]) %>% 
			dplyr::arrange(b) %>%
			dplyr::pull(logFC),
		tolerance=1e-4
	)

})



test_that("filter abundant - SummarizedExperiment",{

	res =		keep_abundant(		se	)

	expect_equal(		nrow(res),		23	)

})

test_that("filter variable - no object",{

	res =		keep_variable(se, top = 5		)

	expect_equal(	nrow(res),5	)

	res =
		keep_variable(
			tidybulk:::tidybulk_to_SummarizedExperiment(input_df, a, b, c),
			top = 5
		)

	expect_equal(	nrow(res),5	)

	expect_equivalent(
		sort(rownames(res)),
		c("FCN1",  "IGHD",  "IGHM",  "IGKC",  "TCL1A")
	)

})
