context('Bulk methods')


test_that("Creating tt object from tibble, number of parameters, methods",{

  expect_equal(

    length(
      attr(
        create_ttBulk(
          ttBulk::counts_mini,
          sample_column = sample,
          transcript_column = transcript,
          counts_column = `count`
        ) ,
        "parameters"
      )
    ),
    3
  )

})

test_that("Test class identity of tt object",{

  expect_equal(
    class(
      create_ttBulk(
        ttBulk::counts_mini,
        sample_column = sample,
        transcript_column = transcript,
        counts_column = `count`
      )
    )[1],
    "ttBulk"
  )

})

test_that("Getting normalised counts - no object",{

  res =
    normalise_counts(
      ttBulk::counts_mini,
      sample_column = sample,
      transcript_column = transcript,
      counts_column = `count`,
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

})


test_that("Adding normalised counts - no object",{

  res =
    normalise_counts(
      ttBulk::counts_mini,
      sample_column = sample,
      transcript_column = transcript,
      counts_column = `count`,
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

test_that("Get differential trancript abundance - no object",{

  res =
    annotate_differential_transcription(
      ttBulk::counts_mini,
      ~ condition,
      sample_column = sample,
      transcript_column = transcript,
      counts_column = `count`,
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

})

test_that("Add differential trancript abundance - no object",{

  res =
    annotate_differential_transcription(
      ttBulk::counts_mini,
      ~ condition,
      sample_column = sample,
      transcript_column = transcript,
      counts_column = `count`,
      action="add"
    )

  expect_equal(
    unique(res$logFC)[1:4],
    c(-12.10269, -12.48201, -11.48896, -13.44406),
    tolerance=1e-6
  )

  expect_equal(
    ncol(res),
    13
  )

})

test_that("Get adjusted counts - no object",{

  cm = ttBulk::counts_mini
  cm$batch = 0
  cm$batch[cm$sample %in% c("SRR1740035", "SRR1740043")] = 1

  res =
    adjust_counts(
      cm,
      ~ condition + batch,
      sample_column = sample,
      transcript_column = transcript,
      counts_column = `count`,
      action="get"
    )

  expect_equal(
    unique(res$`count adjusted`)[c(1, 2, 3, 5)],
    c( 6, 1014,   25 ,   0),
    tolerance=1e-6
  )

  expect_equal(
    ncol(res),
    4
  )

})


test_that("Get cluster lables based on Kmeans - no object",{

  res =
    annotate_clusters(
      ttBulk::counts_mini,
      value_column = `count`,
      elements_column = sample,
      feature_column = transcript,
      number_of_clusters = 2,
      action="get"
    )

  expect_equal(
    typeof(res$cluster),
    "integer"
  )

  expect_equal(
    ncol(res),
    2
  )

})

test_that("Add cluster lables based on Kmeans - no object",{

  res =
    annotate_clusters(
      ttBulk::counts_mini,
      value_column = `count`,
      elements_column = sample,
      feature_column = transcript,
      number_of_clusters = 2,
      action="add"
    )

  expect_equal(
    typeof(res$cluster),
    "integer"
  )

  expect_equal(
    ncol(res),
    7
  )

})

test_that("Get reduced dimensions MDS - no object",{

  res =
    reduce_dimensions(
      ttBulk::counts_mini,
      method = "MDS",
      value_column = `count`,
      elements_column = sample,
      feature_column = transcript,
      action="get"
    )

  expect_equal(
    res$`Dim 1`,
    c(1.4048441,  1.3933490, -2.0138120 , 0.8832354, -1.6676164),
    tolerance=1e-6
  )

  expect_equal(
    ncol(res),
    3
  )

})

test_that("Add reduced dimensions MDS - no object",{

  res =
    reduce_dimensions(
      ttBulk::counts_mini,
      method = "MDS",
      value_column = `count`,
      elements_column = sample,
      feature_column = transcript,
      action="add"
    )

  expect_equal(
    (res$`Dim 1`)[1:4],
    c( 1.404844, 1.404844, 1.404844, 1.404844),
    tolerance=1e-6
  )

  expect_equal(
    ncol(res),
    8
  )

})

test_that("Get reduced dimensions PCA - no object",{

  res =
    reduce_dimensions(
      ttBulk::counts_mini,
      method="PCA",
      value_column = `count`,
      elements_column = sample,
      feature_column = transcript,
      action="get"
    )

  expect_equal(
    res$PC1,
    c( -0.4959163, -0.4977283, -0.4145928 ,-0.3582299, -0.4540019),
    tolerance=1e-6
  )

  expect_equal(
    ncol(res),
    3
  )

})

test_that("Add reduced dimensions PCA - no object",{

  res =
    reduce_dimensions(
      ttBulk::counts_mini,
      method="PCA",
      value_column = `count`,
      elements_column = sample,
      feature_column = transcript,
      action="add"
    )

  expect_equal(
    res$PC1[1:4],
    c( -0.4959163, -0.4959163, -0.4959163, -0.4959163),
    tolerance=1e-6
  )

  expect_equal(
    ncol(res),
    8
  )

})

test_that("Get rotated dimensions - no object",{

  res.pca =
    reduce_dimensions(
      ttBulk::counts_mini,
      method="PCA",
      value_column = `count`,
      elements_column = sample,
      feature_column = transcript,
      action="add"
    )

  res =
    rotate_dimensions(
      res.pca,
      dimension_1_column = PC1,
      dimension_2_column = PC2,
      rotation_degrees = 45,
      elements_column = sample,
      action="get"
    )

  expect_equal(
    res$`PC1 rotated 45`,
    c(-0.08299217 ,-0.08765521 ,-0.71713866 ,-0.03872173 ,-0.68530405),
    tolerance=1e-6
  )

  expect_equal(
    ncol(res),
    3
  )

})

test_that("Add rotated dimensions - no object",{

  res.pca =
    reduce_dimensions(
      ttBulk::counts_mini,
      method="PCA",
      value_column = `count`,
      elements_column = sample,
      feature_column = transcript,
      action="add"
    )

  res =
    rotate_dimensions(
      res.pca,
      dimension_1_column = PC1,
      dimension_2_column = PC2,
      rotation_degrees = 45,
      elements_column = sample,
      action="add"
    )

  expect_equal(
    res$`PC1 rotated 45`[1:4],
    c( -0.08299217, -0.08299217, -0.08299217, -0.08299217),
    tolerance=1e-6
  )

  expect_equal(
    ncol(res),
    10
  )

})

test_that("Aggregate duplicated transcript - no object",{

  res =
    aggregate_duplicates(
      ttBulk::counts_mini,
      sample_column = sample,
      transcript_column = transcript,
      counts_column = `count`
    )

  expect_equal(
    res$transcript[1:4],
    c( "TNFRSF4", "PLCH2" ,  "PADI4" ,  "PAX7"   )
  )

  expect_equal(
    ncol(res),
    7
  )

})

test_that("Drop redundant correlated - no object",{

  res =
    drop_redundant(
      ttBulk::counts_mini,
      method = "correlation",
      value_column = `count`,
      elements_column = sample,
      feature_column = transcript
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
      ttBulk::counts_ensembl,
      ensembl_transcript_column = ens,
      action="get"
    )

  expect_equal(
    res$hgnc_symbol,
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
      ttBulk::counts_ensembl,
      ensembl_transcript_column = ens,
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

  res =
    annotate_cell_type(
      ttBulk::counts_mini,
      sample_column = sample,
      transcript_column = transcript,
      counts_column = `count`,
      action="get", cores=2
    )

  expect_equal(
    as.numeric(res[1,2:5]),
    c(0.6223514, 0.2378625, 0.0000000 ,0.0000000),
    tolerance=1e-6
  )

  expect_equal(
    ncol(res),
    23
  )

})
#
test_that("Add cell type proportions - no object",{

  res =
    annotate_cell_type(
      ttBulk::counts_mini,
      sample_column = sample,
      transcript_column = transcript,
      counts_column = `count`,
      action="add", cores=2
    )

  expect_equal(
    as.numeric(res[1,7:10]),
    c(0.6223514, 0.2378625, 0.0000000 ,0.0000000),
    tolerance=1e-6
  )

  expect_equal(
    ncol(res),
    28
  )

})
