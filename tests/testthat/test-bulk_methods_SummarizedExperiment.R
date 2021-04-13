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

	res2 = tidybulk(se) %>% identify_abundant() %>% scale_abundance()
	
	res %>% distinct(sample, multiplier) %>% pull(multiplier)
	res2 %>% distinct(sample, multiplier) %>% pull(multiplier)
	
	
	expect_equal(
		res %>% distinct(sample, multiplier) %>% pull(multiplier),
		res2 %>% distinct(sample, multiplier) %>% pull(multiplier) %>% as.numeric(),
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
		"PC2_rotated_45"
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
	
	# Treat
	input_df %>% 
		tidybulk:::tidybulk_to_SummarizedExperiment(a, b, c) %>%
		identify_abundant(a, b, c, factor_of_interest = condition) %>%
		test_differential_abundance(
			~ condition,
			.sample = a,
			.transcript = b,
			.abundance = c,
			scaling_method = "TMM",
			method = "edgeR_likelihood_ratio",
			test_above_log2_fold_change = 1,
			action="only"
		) %>%
		`@` (elementMetadata) %>%
		as_tibble() %>%
		filter(FDR<0.05) %>%
		nrow %>%
		expect_equal(169)

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

test_that("impute missing",{
	
	res =
		input_df %>%
		dplyr::slice(-1) %>%
		tidybulk:::tidybulk_to_SummarizedExperiment(a, b, c) %>%
		impute_missing_abundance(	~ condition	)
	
	expect_equal(	assays(res) %>% as.list() %>% .[[1]] %>% .["TNFRSF4", "SRR1740034"],	203.5	)
	

	expect_equal(	nrow(res)*ncol(res),	nrow(input_df)	)
	
})

test_that("differential composition",{
	
	# Cibersort
	input_df %>%
		tidybulk:::tidybulk_to_SummarizedExperiment(a, b, c) %>%
		test_differential_cellularity(. ~ condition	, cores = 1	) %>% 
		pull(`estimate_(Intercept)`) %>%
		.[[1]] %>%
		as.integer %>%
		expect_equal(	-2, 	tollerance =1e-3)
	
	# llsr
	input_df %>%
		tidybulk:::tidybulk_to_SummarizedExperiment(a, b, c) %>%
	test_differential_cellularity(
		. ~ condition,
		method="llsr",
		cores = 1
	) %>% 
		pull(`estimate_(Intercept)`) %>%
		.[[1]] %>%
		as.integer %>%
		expect_equal(	-2, 	tollerance =1e-3)
	
	# Survival analyses
	input_df %>%
		dplyr::select(a, b, c) %>%
		nest(data = -a) %>%
		mutate(
			days = c(1, 10, 500, 1000, 2000),
			dead = c(1, 1, 1, 0, 1)
		) %>%
		unnest(data) %>%
		tidybulk:::tidybulk_to_SummarizedExperiment(a, b, c) %>%
		test_differential_cellularity(
			survival::Surv(days, dead) ~ .,
			cores = 1
		) %>%
		pull(estimate) %>%
		.[[1]] %>%
		#expect_equal(97)
		round() %in% c(
			97,  # 97 is the github action MacOS that has different value
			112, # 112 is the github action UBUNTU that has different value
			93 # 93 is the github action Windows that has different value
		) %>% 
		expect_true()
	
})

test_that("test_stratification_cellularity",{
	
	# Cibersort
	input_df %>%
		select(a, b, c) %>%
		nest(data = -a) %>%
		mutate(
			days = c(1, 10, 500, 1000, 2000),
			dead = c(1, 1, 1, 0, 1)
		) %>%
		unnest(data) %>%
		tidybulk:::tidybulk_to_SummarizedExperiment(a, b, c) %>%
		test_stratification_cellularity(
			survival::Surv(days, dead) ~ .,
			cores = 1
		) %>%
		pull(.low_cellularity_expected) %>%
		.[[1]] %>%
		expect_equal(3.35, tolerance  =1e-1)
	
	# llsr
	input_df %>%
		select(a, b, c) %>%
		nest(data = -a) %>%
		mutate(
			days = c(1, 10, 500, 1000, 2000),
			dead = c(1, 1, 1, 0, 1)
		) %>%
		unnest(data) %>%
		test_stratification_cellularity(
			survival::Surv(days, dead) ~ .,
			.sample = a,
			.transcript = b,
			.abundance = c,
			cores = 1, method = "llsr"
		) %>%
		pull(.low_cellularity_expected) %>%
		.[[1]] %>%
		expect_equal(3.35, tolerance  =1e-1)
})


# test_that("Get gene enrichment - no object",{
# 	
# 	if (find.package("EGSEA", quiet = TRUE) %>% length %>% equals(0)) {
# 		message("Installing EGSEA needed for differential transcript abundance analyses")
# 		if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
# 		BiocManager::install("EGSEA")
# 	}
# 	
# 	library(EGSEA)
# 	
# 	res =
# 		aggregate_duplicates(
# 			dplyr::rename(symbol_to_entrez(
# 				#dplyr::filter(input_df, grepl("^B", b)),
# 				input_df,
# 				.transcript = b, .sample = a), d = entrez
# 			),
# 			.transcript = d,
# 			.sample = a,
# 			.abundance = c
# 		) %>% identify_abundant(a, b, c, factor_of_interest = condition) %>%
# 		tidybulk:::tidybulk_to_SummarizedExperiment(a, b, c) %>%
# 		test_gene_enrichment(
# 			~ condition,
# 			.entrez = d,
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

test_that("pivot",{
	
	expect_equal(	ncol(pivot_sample(tidybulk:::tidybulk_to_SummarizedExperiment(tidybulk(input_df, a, b, c)))	), 4)
	
	expect_equal(	ncol(pivot_transcript(tidybulk:::tidybulk_to_SummarizedExperiment(tidybulk(input_df, a, b, c)))	), 1)
	
})

test_that("gene over representation",{
	
	df_entrez = symbol_to_entrez(tidybulk::counts_mini, .transcript = transcript, .sample = sample)
	df_entrez = aggregate_duplicates(df_entrez, aggregation_function = sum, .sample = sample, .transcript = entrez, .abundance = count)
	df_entrez = mutate(df_entrez, do_test = transcript %in% c("TNFRSF4", "PLCH2", "PADI4", "PAX7"))
	
	res =
		df_entrez %>%
		tidybulk:::tidybulk_to_SummarizedExperiment(sample, transcript, count) %>%
		test_gene_overrepresentation(
			.entrez = entrez,
			.do_test = do_test,
			species="Homo sapiens"
		)
	
	expect_equal(	ncol(res),	10	)
	
	
	
})


test_that("Only reduced dimensions MDS - no object",{
	
	
	

	res =
		tidybulk:::tidybulk_to_SummarizedExperiment(tidybulk(input_df, a, b, c)) %>%
		reduce_dimensions(
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
	
	# Duplicate genes/samples
	expect_error(
		reduce_dimensions(
			input_df %>% identify_abundant(a, b, c) %>% bind_rows( (.) %>% dplyr::slice(1) %>% mutate(c = c+1) ),
			method = "MDS",
			.abundance = c,
			.element = a,
			.feature = b, action="only"
		),
		"Your dataset include duplicated "
	)
})
