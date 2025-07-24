context('Bulk methods SummarizedExperiment')

data("se_mini")
data("breast_tcga_mini_SE")

library(dplyr)
library(SummarizedExperiment)



test_that("tidybulk SummarizedExperiment normalisation",{
  
  res = se |> identify_abundant() |> scale_abundance()
  
  expect_equal(
    names(SummarizedExperiment::assays(res)),
    c("counts" ,"counts_scaled")
  )
  
})

test_that("quantile normalisation",{
  
  res = se_mini |> quantile_normalise_abundance()
  

  

  # preprocessCore
  res = se_mini |> quantile_normalise_abundance(method = "preprocesscore_normalize_quantiles_use_target")
  
 
  

  
  target_distribution = 
    se_mini |> 
    assay( "count") |> 
    as.matrix() |> 
    preprocessCore::normalize.quantiles.determine.target() 
  
  se_mini |> 
    quantile_normalise_abundance(
      method = "preprocesscore_normalize_quantiles_use_target", 
      target_distribution = target_distribution
    ) |> 
    expect_no_error()
  
  
})

test_that("tidybulk SummarizedExperiment normalisation subset",{
  
  res = se |> identify_abundant() |> scale_abundance(
    .subset_for_scaling = .abundant & grepl("^E", .feature)
  )
  
  expect_equal(
    unique(SummarizedExperiment::colData(res)$multiplier),
    c(1.425486, 1.645707, 1.216928, 2.222842, 1.264932, 1.000000, 1.507632, 1.556014),
    tolerance=1e-3
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
  
  res = se |> identify_abundant() |> reduce_dimensions(method="PCA")
  
  expect_equal(
    tail(names(SummarizedExperiment::colData(res)), 1),
    "PC2"
  )
  
})

test_that("Get rotated dimensions - SummarizedExperiment",{
  
  res.pca =
    reduce_dimensions(se |> identify_abundant(),		method="PCA"	)
  
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
  
  cm = se_mini
  cm$batch = 0
  cm$batch[colnames(cm) %in% c("SRR1740035", "SRR1740043")] = 1
  
  res =
    cm |>
    identify_abundant() |>
    adjust_abundance(
      ~ condition + batch,
      method = "combat"
    )
  
  expect_equal(nrow(res),	527	)
  
  expect_equal(	names(SummarizedExperiment::assays(res)),	c("count" ,"count_adjusted")	)
  
  
})

test_that("Get adjusted counts multiple factors - SummarizedExperiment",{
  
  cm = se_mini
  cm$batch = 0
  cm$batch[colnames(cm) %in% c("SRR1740035", "SRR1740043")] = 1
  cm =
    cm |>
    identify_abundant() |>
    scale_abundance()
  cm@assays@data$count_scaled  = apply(cm@assays@data$count_scaled, 2, as.integer)
  
  
  res =
    cm |>
    adjust_abundance(.factor_unwanted = c(batch),
                     .factor_of_interest =  time,
                     .abundance = count_scaled,
                     method = "combat_seq",
                     shrink.disp = TRUE,
                     shrink = TRUE,
                     gene.subset.n = 100
    )
  
  expect_equal(nrow(res),	527	)
  
  expect_equal(	names(SummarizedExperiment::assays(res)),	c("count" , "count_scaled", "count_scaled_adjusted")	)
  
  
})

test_that("Aggregate duplicated transcript - SummarizedExperiment",{
  library(dplyr)
  
  se = tidybulk::se
  SummarizedExperiment::rowData(se)$bla =
    rownames(se) |> purrr::map_chr(~ {
      if(.x  %in% c("LRG_239", "LRG_405")) "BLAAA"
      else .x
    })
  
  
  rowData(se)$bla = if_else(rowData(se)$bla  %in% c("LRG_239", "LRG_405"), "BLAAA", rowData(se)$bla )
  
  res =  aggregate_duplicates(se,	.transcript = bla )
  
  
  expect_equal(
    as.data.frame(rowRanges(se["ENSG00000272397","SRR1039508"])) |>
      mutate_if(is.factor, as.character),
    as.data.frame(rowRanges(res["ENSG00000272397","SRR1039508"]))|>
      mutate_if(is.factor, as.character),
  )
  
  expect_equal(	dim(res),	c( 99,   8  )	)
  
  
})

test_that("Add cell type proportions - SummarizedExperiment",{
  
  res =		deconvolve_cellularity(se_mini, cores=1	)
  
  expect_equal(
    as.numeric(as.data.frame(res@colData[1, 6:9])),
    c(  0.619662   ,     0.25256  ,          0    ,       0),
    tolerance=1e-3
  )
  
})

test_that("differential trancript abundance - SummarizedExperiment",{
  
  res =		test_differential_abundance(
    se_mini |>
      identify_abundant(formula_design = ~  condition),
    ~ condition
  )
  
  w = match(  c("CLEC7A" , "FAM198B", "FCN1"  ,  "HK3"   ), rownames(res) )
  

  expect_equal(
    res@elementMetadata[w,]$logFC,
    c(-11.58385, -13.53406, -12.58204, -12.19271),
    tolerance=1e-3
  )
  

  
  # Likelihood ratio
  res2 =		test_differential_abundance(
    se_mini |>
      identify_abundant(formula_design = ~  condition),
    ~ condition, method = "edgeR_likelihood_ratio"	)
  

  expect_equal(
    res2@elementMetadata[w,]$logFC,
    c(-11.57989, -13.53476, -12.57969, -12.19303),
    tolerance=1e-3
  )
  

  
  # Treat
  se_mini |>
    identify_abundant(formula_design = ~  condition) |>
    test_differential_abundance(
      ~ condition,
      scaling_method = "TMM",
      method = "edgeR_likelihood_ratio",
      test_above_log2_fold_change = 1,
      action="only"
    ) |> SummarizedExperiment::rowData() |>
    as_tibble() |>
    filter(FDR<0.05) |>
    nrow() |>
    expect_equal(171)
  
})

test_that("differential trancript abundance - SummarizedExperiment - alternative .abundance",{
  
  library(SummarizedExperiment)
  
  assays(se_mini) = list(counts = assay(se_mini), bla = assay(se_mini))
  
  
  res =	 se_mini |>
    identify_abundant(formula_design = ~  condition) |>
    test_differential_abundance(   ~ condition , .abundance =  bla)
  
  w = match(  c("CLEC7A" , "FAM198B", "FCN1"  ,  "HK3"   ), rownames(res) )
  

  expect_equal(
    res@elementMetadata[w,]$logFC,
    c(-11.58385, -13.53406, -12.58204, -12.19271),
    tolerance=1e-3
  )
  

  # Likelihood ratio
  res2 =		test_differential_abundance(
    se_mini |>
      identify_abundant(formula_design = ~  condition),
    ~ condition, .abundance =  bla, method = "edgeR_likelihood_ratio"	)
  

  
  expect_equal(
    res2@elementMetadata[w,]$logFC,
    c(-11.57989, -13.53476, -12.57969, -12.19303),
    tolerance=1e-3
  )
  

  
  # Treat
  se_mini |>
    identify_abundant( formula_design = ~  condition) |>
    test_differential_abundance(
      ~ condition, .abundance =  bla,
      scaling_method = "TMM",
      method = "edgeR_likelihood_ratio",
      test_above_log2_fold_change = 1,
      action="only"
    ) |> SummarizedExperiment::rowData() |>
    as_tibble() |>
    filter(FDR<0.05) |>
    nrow() |>
    expect_equal(171)
  
})

test_that("DE interaction effects", {
  
  # Att interaction factor
  col_data = colData(breast_tcga_mini_SE)
  set.seed(42)
  col_data$interaction_term = sample(c(0,1), size = nrow(col_data), replace = TRUE)
  colData(breast_tcga_mini_SE) = col_data
  
  breast_tcga_mini_SE |> 
    identify_abundant(formula_design = ~  Call) |>
    test_differential_abundance( 
      ~ Call * interaction_term, 
      contrasts = "CallHer2___interaction_term",
      method="edgeR_quasi_likelihood" 
    ) |> 
    expect_no_error()
  
})

test_that("Voom with treat method",{
  
  se_mini |>
    identify_abundant(formula_design = ~ condition) |>
    test_differential_abundance(
      ~ condition,
      method = "limma_voom",
      test_above_log2_fold_change = 1,
      action="only"
    ) |> SummarizedExperiment::rowData() |>
    as_tibble() |>
    filter(adj.P.Val<0.05) |>
    nrow() |>
    expect_equal(97)
  
  # with multiple contrasts
  res <-
    se_mini |>
    identify_abundant( formula_design = ~ Cell.type) |>
    test_differential_abundance(
      ~ 0 + Cell.type,
      .contrasts = c("Cell.typeb_cell-Cell.typemonocyte", "Cell.typeb_cell-Cell.typet_cell"),
      method = "limma_voom",
      test_above_log2_fold_change = 1,
      action="only"
    ) |> SummarizedExperiment::rowData() |>
    as_tibble()
  
  res |>
    filter(adj.P.Val___Cell.typeb_cell.Cell.typemonocyte < 0.05) |>
    nrow() |>
    expect_equal(294)
  
  res |>
    filter(adj.P.Val___Cell.typeb_cell.Cell.typet_cell < 0.05) |>
    nrow() |>
    expect_equal(246)
  
})

test_that("differential trancript abundance - random effects SE",{
  
  set.seed(42)
  res =
    se_mini[1:10,] |>
    identify_abundant(formula_design = ~  condition) |>
    #mutate(time = time |> stringr::str_replace_all(" ", "_")) |>
    test_differential_abundance(
      ~ condition + (1 + condition | time),
      method = "glmmseq_lme4",
      cores = 1
    )
  
  rowData(res)[,"P_condition_adjusted"] |>
    head(4) |>
    expect_equal(
      c(0.1578695, 0.1221392, 0.1221392, 0.2262688),
      tolerance=1e-2
    )
  
  # Custom dispersion
  se_mini =
    se_mini |>
    identify_abundant(formula_design = ~  condition)
  
  rowData(se_mini)$disp_ = rep(2,nrow(se_mini))
  
  res =
    se_mini[1:10,] |>
    #mutate(time = time |> stringr::str_replace_all(" ", "_")) |>
    test_differential_abundance(
      ~ condition + (1 + condition | time),
      method = "glmmseq_lme4",
      .dispersion = disp_,
      cores = 1
    )
  
  rowData(res)[,"P_condition_adjusted"] |>
    head(4) |>
    expect_equal(
      c(0.2633982, 0.2633982, 0.2633982, 0.5028348),
      tolerance=1e-2
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
      se_mini,
      top = 5
    )
  
  expect_equal(	nrow(res),5	)
  
  expect_equivalent(
    sort(rownames(res)),
    c("FCN1",  "IGHD",  "IGHM",  "IGKC",  "TCL1A")
  )
  
})



test_that("differential composition",{
  
  # Cibersort
  se_mini |>
    test_differential_cellularity(. ~ condition	, cores = 1	) |>
    pull(`estimate_(Intercept)`) |>
    magrittr::extract2(1) |>
    as.integer() |>
    expect_equal(	-2, 	tollerance =1e-3)
  
  # llsr
  se_mini |>
    test_differential_cellularity(
      . ~ condition,
      method="llsr"
    ) |>
    pull(`estimate_(Intercept)`) |>
    magrittr::extract2(1) |>
    as.integer() |>
    expect_equal(	-2, 	tollerance =1e-3)
  
 
  
})




# test_that("Get gene enrichment - no object",{
#
# 	if (find.package("EGSEA", quiet = TRUE) |> length |> equals(0)) {
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
# 		) |> identify_abundant(a, b, c, formula_design = ~  condition) |>
# 		tidybulk:::tidybulk_to_SummarizedExperiment(a, b, c) |>
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
  
  expect_equal(	ncol(pivot_sample(se_mini)	), 6)
  
  expect_equal(	ncol(pivot_transcript(se_mini)	), 2)
  
})



test_that("Only reduced dimensions MDS - no object",{
  
  
  
  
  res =
    breast_tcga_mini_SE |>
    reduce_dimensions(method = "MDS")
  
  expect_equal(
    res$`Dim1`[1:4],
    c(-0.2723808836, -0.1105770207, -0.3034092668, -0.0064569358),
    tolerance=10
  )
  
  expect_equal(
    ncol(SummarizedExperiment::colData(res)),
    3
  )
  
  expect_equal(	class(attr(res, "internals")$MDS[[1]])[1], 	"MDS"  )
  
})

test_that("Only reduced dimensions PCA - no object",{
  
  
  
  
  res =
    breast_tcga_mini_SE |>
    reduce_dimensions(  method = "PCA"  )
  
  expect_equal(
    res$`PC1`[1:4],
    c(-8.5145575 ,  4.2451190 , -6.8991306 , -6.0787597 ),
    tolerance=10
  )
  
  expect_equal(
    ncol(SummarizedExperiment::colData(res)),
    3
  )
  
  expect_equal(	class(attr(res, "internals")$PCA[[1]])[1], 	"numeric"  )
  
})

test_that("Only reduced dimensions tSNE - no object",{
  
  
  
  
  res =
    breast_tcga_mini_SE |>
    reduce_dimensions(  method = "tSNE"  )
  
  expect_equal(
    res$`tSNE1`[1:4],
    c(-2.4677708, -0.3447403, -5.8059667,  1.9685007 ),
    tolerance=10
  )
  
  expect_equal(
    ncol(SummarizedExperiment::colData(res)),
    3
  )
  
  expect_equal(	class(attr(res, "internals")$tSNE[[1]])[1], 	"integer"  )
  
})

test_that("Only reduced dimensions UMAP - no object",{
  
  
  
  
  res =
    breast_tcga_mini_SE |>
    reduce_dimensions(  method = "UMAP"  )
  
  expect_equal(
    res$`UMAP1`[1:4],
    c(-1.9784263, -0.8536380, -0.6537955, -1.8840860 ),
    tolerance=10
  )
  
  expect_equal(
    ncol(SummarizedExperiment::colData(res)),
    3
  )
  
  expect_equal(	class(attr(res, "internals")$UMAP[[1]])[1], 	"numeric"  )
  
})


test_that("resolve_complete_confounders_of_non_interest",{
  
  
  #library(tidySummarizedExperiment)
  library(SummarizedExperiment)
  library(dplyr)
  
  # Sample annotations
  sample_annotations <- data.frame(
    sample_id = paste0("Sample", seq(1, 9)),
    factor_of_interest =  c(rep("treated", 4), rep("untreated", 5)),
    A = c("a1", "a2", "a1", "a2", "a1", "a2", "a1", "a2", "a3"),
    B = c("b1", "b1", "b2", "b1", "b1", "b1", "b2", "b1", "b3"),
    C = c("c1", "c1", "c1", "c1", "c1", "c1", "c1", "c1", "c3"),
    stringsAsFactors = FALSE
  )
  
  # Simulated assay data (e.g., gene expression data)
  # Let's assume we have 100 genes (rows) and 9 samples (columns)
  assay_data <- matrix(rnorm(100 * 9), nrow = 100, ncol = 9)
  
  # Row data (e.g., gene annotations)
  # For simplicity, we'll just use a sequence of gene IDs
  row_data <- data.frame(gene_id = paste0("Gene", seq_len(100)))
  
  # Create SummarizedExperiment object
  se <- SummarizedExperiment(assays = list(counts = assay_data),
                             rowData = row_data,
                             colData = DataFrame(sample_annotations))
  
  expected_tibble =
    tibble(
      .sample = c("1", "5", "2", "4", "6", "8", "9", "3", "7"),
      A = c("a1", "a1", "a2", "a2", "a2", "a2", "a3", "a1", "a1"),
      B = c("b1", "b1", "b1", "b1", "b1", "b1", "b1", "b2", "b2"),
      C = rep("c1", 9)
    ) |> arrange(.sample)
  
  se_deconfounded = 
    se |>
    resolve_complete_confounders_of_non_interest(A, B, C) 
  
  se_deconfounded |>
    colData() |> 
    _[, c("A___altered", "B___altered", "C___altered")] |>
    as_tibble(rownames = ".sample") |> 
    setNames(c(".sample", "A", "B", "C")) |> 
    expect_identical(expected_tibble )
  
  
  se_deconfounded |>
    colData() |> 
    ncol() |> 
    expect_equal(se |> colData() |> ncol() + 3 )
  
})
