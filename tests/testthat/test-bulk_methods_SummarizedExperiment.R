context('Bulk methods SummarizedExperiment')

input_df = setNames(ttBulk::counts_mini, c("a", "b", "Cell type", "c",  "time" , "condition"))

input_df_breast = setNames(ttBulk::breast_tcga_mini, c("a", "b", "c norm", "call", "c"))

test_that("ttBulk SummarizedExperiment conversion",{

	res = ttBulk(ttBulk::se)

	expect_equal(	class(res)[1],	"ttBulk"	)

	expect_equal(	nrow(res),	800	)

	expect_equal(	ncol(res),	12	)

	res = res %>% ttBulk:::ttBulk_to_SummarizedExperiment()

	expect_equal(	class(res)[1],	"SummarizedExperiment"	)

	expect_equal(	nrow(res),	100	)

	expect_equal(	ncol(res),	8	)

})

test_that("ttBulk SummarizedExperiment normalisation manual",{

	res = ttBulk(ttBulk:::ttBulk_to_SummarizedExperiment(scale_abundance(ttBulk(se))))

	expect_equal(
		res[1:4,]$`counts scaled`,
		c(1327.286584 , 859.307120 , 898.641600  ,  2.017153),
		tolerance=1e-6
	)

	expect_equal(	nrow(res),	800	)

	expect_equal(	ncol(res),	16	)


	res = rlang::quo_name(attr(res, "parameters")[[4]])

	expect_equal( res,	"counts scaled"	)

})

test_that("ttBulk SummarizedExperiment normalisation",{

	res = scale_abundance(se)

	expect_equal(
		names(SummarizedExperiment::assays(res)),
		c("counts" ,"counts scaled")
	)

})

test_that("ttBulk SummarizedExperiment clustering",{

	res = cluster_elements(se, method="kmeans", centers = 2)

	expect_equal(
		tail(names(SummarizedExperiment::colData(res)), 1),
		"cluster.kmeans"
	)

	expect_equal(
		levels(SummarizedExperiment::colData(res)$cluster.kmeans),
		c("1", "2")
	)

})

test_that("ttBulk SummarizedExperiment clustering",{

	res = reduce_dimensions(se, method="PCA")

	expect_equal(
		tail(names(SummarizedExperiment::colData(res)), 1),
		"PC2"
	)

})

test_that("Get rotated dimensions - SummarizedExperiment",{

	res.pca =
		reduce_dimensions(se,		method="PCA"	)

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
			ttBulk:::ttBulk_to_SummarizedExperiment(cm, a, b, c),
			~ condition + batch
		)

	expect_equal(nrow(res),	527	)

	expect_equal(	names(SummarizedExperiment::assays(res)),	c("c" ,"c adjusted")	)


})

test_that("Aggregate duplicated transcript - SummarizedExperiment",{

	res =	aggregate_duplicates(	se)

	expect_equal(	dim(res),	c( 100,   8  )	)


})

test_that("Add cell type proportions - SummarizedExperiment",{

	res =		deconvolve_cellularity(ttBulk:::ttBulk_to_SummarizedExperiment(input_df, a, b, c), cores=1	)

	expect_equal(
		as.numeric(as.data.frame(res@colData[1, 4:7])),
		c( 0.6221895, 0.2380869, 0.0000000, 0.0000000),
		tolerance=1e-3
	)

})

test_that("Add differential trancript abundance - SummarizedExperiment",{

	res =		test_differential_abundance(ttBulk:::ttBulk_to_SummarizedExperiment(input_df, a, b, c),	~ condition	)

	w = match(  c("HK3", "FCN1", "CLEC7A", "FAM198B"), rownames(res) )

	expect_equal(
		res@elementMetadata[w,]$logFC,
		c(-12.10269, -12.48201 ,-11.48896, -13.44406),
		tolerance=1e-6
	)


})

test_that("filter abundant - SummarizedExperiment",{

	res =		filter_abundant(		se	)

	expect_equal(		nrow(res),		31	)

})

test_that("filter variable - no object",{

	res =		filter_variable(se, top = 5		)

	expect_equal(	nrow(res),5	)

	res =
		filter_variable(
			ttBulk:::ttBulk_to_SummarizedExperiment(input_df, a, b, c),
			top = 5
		)

	expect_equal(	nrow(res),5	)

	expect_equivalent(
		sort(rownames(res)),
		c("FCN1",  "IGHD",  "IGHM",  "IGKC",  "TCL1A")
	)

})
