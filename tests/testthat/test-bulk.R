context('Bulk functions')

test_that("Test data frame",{ expect_equal( ncol(ttBulk::counts_mini), 6 ) })

test_that("Creating tt object from tibble, number of parameters",{

  expect_equal(

    length(
      attr(
        create_tt_from_tibble_bulk(
          ttBulk::counts_mini,
          .sample = sample,
          .transcript = transcript,
          .abundance = `count`
        ) ,
      "parameters"
      )
    ),
    3
  )

})

test_that("Getting normalised counts - no object",{

  res =
    get_normalised_counts_bulk(
      ttBulk::counts_mini,
      .sample = sample,
      .transcript = transcript,
      .abundance = `count`
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
    add_normalised_counts_bulk(
      ttBulk::counts_mini,
      .sample = sample,
      .transcript = transcript,
      .abundance = `count`
    )

  expect_equal(
    unique(res$multiplier),
    c(1.124713, 1.011405, 1.511470, 0.828865, 1.191472),
    tolerance=1e-6
  )

  expect_equal(
    ncol(res),
    10
  )

})

test_that("Get differential trancript abundance - no object",{

  res =
    get_differential_transcript_abundance_bulk(
      ttBulk::counts_mini,
      ~ condition,
      .sample = sample,
      .transcript = transcript,
      .abundance = `count`
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

test_that("Get cluster lables based on Kmeans - no object",{

  res =
    get_clusters_kmeans_bulk(
      ttBulk::counts_mini,
      .abundance = `count`,
      .element = sample,
      .feature = transcript,
      centers = 2
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
    add_clusters_kmeans_bulk(
      ttBulk::counts_mini,
      .abundance = `count`,
      .element = sample,
      .feature = transcript,
      centers = 2
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



test_that("Get reduced dimensions MDS - no object",{

  res =
    get_reduced_dimensions_MDS_bulk(
      ttBulk::counts_mini,
      .abundance = `count`,
      .element = sample,
      .feature = transcript
    )

  expect_equal(
    res$`Dim 1`,
    c(1.4048441 , 1.3933490, -2.0138120 , 0.8832354 ,-1.6676164),
    tolerance=1e-6
  )

  expect_equal(
    ncol(res),
    3
  )

})

test_that("Add reduced dimensions MDS - no object",{

  res =
    add_reduced_dimensions_MDS_bulk(
      ttBulk::counts_mini,
      .abundance = `count`,
      .element = sample,
      .feature = transcript
    )

  expect_equal(
    (res$`Dim 1`)[1:4],
    c(1.404844, 1.404844, 1.404844, 1.404844),
    tolerance=1e-6
  )

  expect_equal(
    ncol(res),
    8
  )

})

test_that("Get reduced dimensions PCA - no object",{

  res =
    get_reduced_dimensions_PCA_bulk(
      ttBulk::counts_mini,
      .abundance = `count`,
      .element = sample,
      .feature = transcript
    )

  expect_equal(
    res$PC1,
    c( -0.4959163, -0.4977283 ,-0.4145928, -0.3582299, -0.4540019),
    tolerance=1e-1
  )

  expect_equal(
    ncol(res),
    3
  )

})

test_that("Add reduced dimensions PCA - no object",{

  res =
    add_reduced_dimensions_PCA_bulk(
      ttBulk::counts_mini,
      .abundance = `count`,
      .element = sample,
      .feature = transcript
    )

  expect_equal(
    res$PC1[1:4],
    c(-0.4959163 ,-0.4959163, -0.4959163, -0.4959163),
    tolerance=1e-1
  )

  expect_equal(
    ncol(res),
    8
  )

})

test_that("Get rotated dimensions - no object",{

  res.pca =
    add_reduced_dimensions_PCA_bulk(
      ttBulk::counts_mini,
      .abundance = `count`,
      .element = sample,
      .feature = transcript
    )

  res =
    get_rotated_dimensions(
      res.pca,
      dimension_1_column = PC1,
      dimension_2_column = PC2,
      rotation_degrees = 45,
      .element = sample
    )

  expect_equal(
    res$`PC1 rotated 45`,
    c(-0.08299217 ,-0.08765521, -0.71713866, -0.03872173, -0.68530405),
    tolerance=1e-1
  )

  expect_equal(
    ncol(res),
    3
  )

})

test_that("Add rotated dimensions - no object",{

  res.pca =
    add_reduced_dimensions_PCA_bulk(
      ttBulk::counts_mini,
      .abundance = `count`,
      .element = sample,
      .feature = transcript
    )

  res =
    add_rotated_dimensions(
      res.pca,
      dimension_1_column = PC1,
      dimension_2_column = PC2,
      rotation_degrees = 45,
      .element = sample
    )

  expect_equal(
    res$`PC1 rotated 45`[1:4],
    c(-0.08299217, -0.08299217, -0.08299217, -0.08299217),
    tolerance=1e-1
  )

  expect_equal(
    ncol(res),
    10
  )

})

test_that("Aggregate duplicated transcript - no object",{

  res =
    aggregate_duplicated_transcripts_bulk(
      ttBulk::counts_mini,
      .sample = sample,
      .transcript = transcript,
      .abundance = `count`
    )

  expect_equal(
    res$transcript[1:4],
    c("TNFRSF4", "PLCH2" ,  "PADI4" ,  "PAX7"   )
  )

  expect_equal(
    ncol(res),
    7
  )

})

test_that("Drop redundant correlated - no object",{

  res =
    remove_redundancy_elements_through_correlation(
      ttBulk::counts_mini,
      .abundance = `count`,
      .element = sample,
      .feature = transcript
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
    get_symbol_from_ensembl(
      ttBulk::counts_ensembl,
      .ensembl = ens
    )

  expect_equal(
    res$transcript,
    "TSPAN6"
  )

  expect_equal(
    ncol(res),
    3
  )

})

test_that("Add symbol from ensambl - no object",{

  res =
    add_symbol_from_ensembl(
      ttBulk::counts_ensembl,
      .ensembl = ens
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

# test_that("Get cell type proportions - no object",{
#
#   expect_error(
#     get_cell_type_proportions(
#       ttBulk::counts_mini,
#       .sample = sample,
#       .transcript = transcript,
#       .abundance = `count`,
#       cores = 2
#     ),
#     "You have less than 50 genes that overlap the Cibersort signature"
#   )
#
#   res =
#     get_cell_type_proportions(
#       ttBulk::counts[
#         ttBulk::counts$sample %in%
#           unique(ttBulk::counts_mini$sample),
#         ],
#       .sample = sample,
#       .transcript = transcript,
#       .abundance = `count`,
#       cores = 2
#     )
#
#   expect_equal(
#     res$proportion[1:4],
#     c(0.6273726, 0.6004113, 0.5931743, 0.5811784),
#     tolerance=1e-6
#   )
#
#   expect_equal(
#     ncol(res),
#     3
#   )
#
# })
#
# test_that("Add cell type proportions - no object",{
#
#   expect_error(
#     add_cell_type_proportions(
#       ttBulk::counts_mini,
#       .sample = sample,
#       .transcript = transcript,
#       .abundance = `count`,
#       cores = 2
#     ),
#     "You have less than 50 genes that overlap the Cibersort signature"
#   )
#
#   res =
#     add_cell_type_proportions(
#       ttBulk::counts[
#         ttBulk::counts$sample %in%
#           unique(ttBulk::counts_mini$sample),
#         ],
#       .sample = sample,
#       .transcript = transcript,
#       .abundance = `count`,
#       cores = 2
#     )
#
#   expect_equal(
#     res$proportion[1:4],
#     c(0.6273726, 0.2553588, 0.0000000, 0.0000000),
#     tolerance=1e-6
#   )
#
#   expect_equal(
#     ncol(res),
#     7
#   )
#
# })
