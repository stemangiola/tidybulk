context('Bulk methods SummarizedExperiment')

data("se_mini")
data("breast_tcga_mini_SE")

input_df =  se_mini |> tidybulk() |> as_tibble() |> setNames( c( "b","a",  "c", "Cell type",  "time" , "condition", "days",  "dead", "entrez"))
input_df_breast = setNames( breast_tcga_mini_SE |> tidybulk() |> as_tibble(), c( "b", "a","c", "c norm", "call" ))

test_that("tidybulk SummarizedExperiment conversion",{

  res = tidybulk(tidybulk::se)

  expect_equal(	class(res)[1],	"tidybulk"	)

  expect_equal(	nrow(res),	800	)

  expect_equal(	ncol(res),	13	)

  res = res |> tidybulk:::tidybulk_to_SummarizedExperiment()

  expect_equal(	class(res)[1],	"SummarizedExperiment"	)

  expect_equal(	nrow(res),	100	)

  expect_equal(	ncol(res),	8	)

})

test_that("tidybulk SummarizedExperiment normalisation manual",{

  res = tidybulk(tidybulk:::tidybulk_to_SummarizedExperiment(scale_abundance(tidybulk(se) |> identify_abundant())))

  res2 = tidybulk(se) |> identify_abundant() |> scale_abundance()

  res |> distinct(.sample, multiplier) |> pull(multiplier)
  res2 |> distinct(.sample, multiplier) |> pull(multiplier)


  expect_equal(
    res |> distinct(.sample, multiplier) |> pull(multiplier),
    res2 |> distinct(.sample, multiplier) |> pull(multiplier) |> as.numeric(),
    tolerance=1e-3
  )

  expect_equal(	nrow(res),	800	)

  expect_equal(	ncol(res),	17	)


  res = rlang::quo_name(attr(res, "internals")$tt_columns[[4]])

  expect_equal( res,	"counts_scaled"	)

})

test_that("tidybulk SummarizedExperiment normalisation",{

  res = se |> identify_abundant() |> scale_abundance()

  expect_equal(
    names(SummarizedExperiment::assays(res)),
    c("counts" ,"counts_scaled")
  )

})


test_that("quantile normalisation",{

  res = se_mini |> quantile_normalise_abundance()

  res_tibble =
    input_df |>
    quantile_normalise_abundance(
      .sample = a,
      .transcript = b,
      .abundance = c,
      action = "get"
    )


    SummarizedExperiment::assay(res, "count_scaled")["ABCB9","SRR1740035"] |>
  expect_equal(
    res_tibble |>
      filter(a=="SRR1740035" & b=="ABCB9") |>
      dplyr::pull(c_scaled)
  )

    # preprocessCore
    res = se_mini |> quantile_normalise_abundance(method = "preprocesscore_normalize_quantiles_use_target")

    res_tibble =
      input_df |>
      quantile_normalise_abundance(
        .sample = a,
        .transcript = b,
        .abundance = c,
        method = "preprocesscore_normalize_quantiles_use_target",
        action = "get"
      )


    SummarizedExperiment::assay(res, "count_scaled")["ABCB9","SRR1740035"] |>
      expect_equal(
        res_tibble |>
          filter(a=="SRR1740035" & b=="ABCB9") |>
          dplyr::pull(c_scaled)
      )

})


test_that("tidybulk SummarizedExperiment normalisation subset",{

  res = se |> identify_abundant() |> scale_abundance(
    .subset_for_scaling = .abundant & grepl("^E", .feature)
  )

  expect_equal(
    unique(SummarizedExperiment::colData(res)$multiplier),
    c(1.3648110, 1.5756592, 1.1651309, 2.1282288, 1.2110911, 0.9574359, 1.4434610, 1.4897840),
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
                     .factor_of_interest = time,
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
      identify_abundant(factor_of_interest = condition),
    ~ condition
  )

  w = match(  c("CLEC7A" , "FAM198B", "FCN1"  ,  "HK3"   ), rownames(res) )

  # Quasi likelihood
  res_tibble =		test_differential_abundance(
    input_df |> identify_abundant(a, b, c, factor_of_interest = condition),
    ~ condition	,
    a, b, c
  )

  expect_equal(
    res@elementMetadata[w,]$logFC,
    c(-11.58385, -13.53406, -12.58204, -12.19271),
    tolerance=1e-3
  )

  expect_equal(
    res@elementMetadata[w,]$logFC,
    res_tibble |>
      pivot_transcript(b) |>
      filter(b %in% rownames(res)[w]) |>
      dplyr::arrange(b) |>
      dplyr::pull(logFC),
    tolerance=1e-3
  )

  # Likelihood ratio
  res2 =		test_differential_abundance(
    se_mini |>
      identify_abundant(factor_of_interest = condition),
    ~ condition, method = "edgeR_likelihood_ratio"	)

  res2_tibble =		test_differential_abundance(
    input_df |> identify_abundant(a, b, c, factor_of_interest = condition),
    ~ condition	,
    a, b, c, method = "edgeR_likelihood_ratio"
  )

  expect_equal(
    res2@elementMetadata[w,]$logFC,
    c(-11.57989, -13.53476, -12.57969, -12.19303),
    tolerance=1e-3
  )

  expect_equal(
    res2@elementMetadata[w,]$logFC,
    res2_tibble |>
      pivot_transcript(b) |>
      filter(b %in% rownames(res)[w]) |>
      dplyr::arrange(b) |>
      dplyr::pull(logFC),
    tolerance=1e-3
  )

  # Treat
  se_mini |>
    identify_abundant(factor_of_interest = condition) |>
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
    expect_equal(169)

})


test_that("differential trancript abundance - SummarizedExperiment - alternative .abundance",{

  assays(se_mini) = list(counts = assay(se_mini), bla = assay(se_mini))


  res =	 se_mini |>
    identify_abundant(factor_of_interest = condition) |>
    test_differential_abundance(   ~ condition , .abundance =  bla)

  w = match(  c("CLEC7A" , "FAM198B", "FCN1"  ,  "HK3"   ), rownames(res) )

  # Quasi likelihood
  res_tibble =		test_differential_abundance(
    input_df |> identify_abundant(a, b, c, factor_of_interest = condition),
    ~ condition	,
    a, b, c
  )

  expect_equal(
    res@elementMetadata[w,]$logFC,
    c(-11.58385, -13.53406, -12.58204, -12.19271),
    tolerance=1e-3
  )

  expect_equal(
    res@elementMetadata[w,]$logFC,
    res_tibble |>
      pivot_transcript(b) |>
      filter(b %in% rownames(res)[w]) |>
      dplyr::arrange(b) |>
      dplyr::pull(logFC),
    tolerance=1e-3
  )

  # Likelihood ratio
  res2 =		test_differential_abundance(
    se_mini |>
      identify_abundant(factor_of_interest = condition),
    ~ condition, .abundance =  bla, method = "edgeR_likelihood_ratio"	)

  res2_tibble =		test_differential_abundance(
    input_df |> identify_abundant(a, b, c, factor_of_interest = condition),
    ~ condition	,
    a, b, c, method = "edgeR_likelihood_ratio"
  )

  expect_equal(
    res2@elementMetadata[w,]$logFC,
    c(-11.57989, -13.53476, -12.57969, -12.19303),
    tolerance=1e-3
  )

  expect_equal(
    res2@elementMetadata[w,]$logFC,
    res2_tibble |>
      pivot_transcript(b) |>
      filter(b %in% rownames(res)[w]) |>
      dplyr::arrange(b) |>
      dplyr::pull(logFC),
    tolerance=1e-3
  )

  # Treat
  se_mini |>
    identify_abundant( factor_of_interest = condition) |>
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
    expect_equal(169)

})


test_that("Voom with treat method",{

  se_mini |>
    identify_abundant(factor_of_interest = condition) |>
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
    identify_abundant( factor_of_interest = Cell.type) |>
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
    expect_equal(293)

  res |>
    filter(adj.P.Val___Cell.typeb_cell.Cell.typet_cell < 0.05) |>
    nrow() |>
    expect_equal(246)

})

test_that("differential trancript abundance - random effects SE",{

 res =
   se_mini[1:10,] |>
    identify_abundant(factor_of_interest = condition) |>
    #mutate(time = time |> stringr::str_replace_all(" ", "_")) |>
    test_differential_abundance(
      ~ condition + (1 + condition | time),
      method = "glmmseq_lme4",
      cores = 1
    )

 rowData(res)[,"P_condition_adjusted"] |>
    head(4) |>
    expect_equal(
      c(0.03394914, 0.03394914, 0.03394914,  NA),
      tolerance=1e-2
    )

 # Custom dispersion
 se_mini =
   se_mini |>
   identify_abundant(factor_of_interest = condition)

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
     c(0.1153254, 0.1668555, 0.1668555 ,       NA),
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

test_that("impute missing",{

  res =
    input_df |>
    dplyr::slice(-1) |>
    tidybulk:::tidybulk_to_SummarizedExperiment(a, b, c) |>
    impute_missing_abundance(	~ condition	)

  list_SE = SummarizedExperiment::assays(res) |> as.list()

  list_SE[[1]]["TNFRSF4", "SRR1740034"] |>
    expect_equal(6)


  expect_equal(	nrow(res)*ncol(res),	nrow(input_df)	)

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

  # Survival analyses
  input_df |>
    dplyr::select(a, b, c) |>
    nest(data = -a) |>
    mutate(
      days = c(1, 10, 500, 1000, 2000),
      dead = c(1, 1, 1, 0, 1)
    ) |>
    unnest(data) |>
    tidybulk:::tidybulk_to_SummarizedExperiment(a, b, c) |>
    test_differential_cellularity(
      survival::Surv(days, dead) ~ .,
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
    tidybulk:::tidybulk_to_SummarizedExperiment(a, b, c) |>
    test_stratification_cellularity(
      survival::Surv(days, dead) ~ .,
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
# 		) |> identify_abundant(a, b, c, factor_of_interest = condition) |>
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

test_that("gene over representation",{

  df_entrez = aggregate_duplicates(input_df, aggregation_function = sum, .sample = a, .transcript = entrez, .abundance = c)
  df_entrez = mutate(df_entrez, do_test = b %in% c("TNFRSF4", "PLCH2", "PADI4", "PAX7"))

  res =
    df_entrez |>
    tidybulk:::tidybulk_to_SummarizedExperiment(a, b, c) |>
    test_gene_overrepresentation(
      .entrez = entrez,
      .do_test = do_test,
      species="Homo sapiens"
    )

  expect_equal(	ncol(res),	10	)



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
