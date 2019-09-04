context('Bulk methods')


test_that("Creating tt object from tibble, number of parameters, methods",{

  expect_equal(

    length(
      attr(
        create_ttBulk(
          ttBulk::counts_mini,
          sample_column = sample,
          transcript_column = transcript,
          counts_column = `read count`
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
        counts_column = `read count`
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
      counts_column = `read count`,
      action = "get"
    )

  expect_equal(
    unique(res$multiplier),
    c(1.0057758, 0.8828939, 0.8369487, 1.1423585),
    tolerance=1e-7
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
      counts_column = `read count`,
      action = "add"
    )

  expect_equal(
    unique(res$multiplier),
    c(1.0057758, 0.8828939, 0.8369487, 1.1423585),
    tolerance=1e-7
  )

  expect_equal(
    ncol(res),
    11
  )

})

test_that("Get differential trancript abundance - no object",{

  res =
    annotate_differential_transcription(
      ttBulk::counts_mini,
      ~ condition,
      sample_column = sample,
      transcript_column = transcript,
      counts_column = `read count`,
      action="get"
    )

  expect_equal(
    unique(res$logFC)[1:4],
    c(-0.5642906,  0.6875711,  0.1718000, -0.4769712),
    tolerance=1e-7
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
      counts_column = `read count`,
      action="add"
    )

  expect_equal(
    unique(res$logFC)[1:4],
    c(-0.56429058,  0.17179998, -0.47697124,  0.08888069),
    tolerance=1e-7
  )

  expect_equal(
    ncol(res),
    14
  )

})

# test_that("Get adjusted counts - no object",{
#
#   res =
#     adjust_counts(
#       ttBulk::counts_mini,
#       ~ condition + batch,
#       sample_column = sample,
#       transcript_column = transcript,
#       counts_column = `read count`,
#       action="get"
#     )
#
#   expect_equal(
#     unique(res$`read count adjusted`)[c(1, 3, 4, 5)],
#     c(372 ,  19 ,  14, 1209),
#     tolerance=1e-7
#   )
#
#   expect_equal(
#     ncol(res),
#     4
#   )
#
# })


test_that("Get cluster lables based on Kmeans - no object",{

  res =
    annotate_clusters(
      ttBulk::counts_mini,
      value_column = `read count`,
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
      value_column = `read count`,
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
    8
  )

})

test_that("Get reduced dimensions MDS - no object",{

  res =
    reduce_dimensions(
      ttBulk::counts_mini,
      method = "MDS",
      value_column = `read count`,
      elements_column = sample,
      feature_column = transcript,
      action="get"
    )

  expect_equal(
    res$`Dimension 1`,
    c(0.2772910, -0.2172255, -0.1806345,  0.1205691),
    tolerance=1e-7
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
      value_column = `read count`,
      elements_column = sample,
      feature_column = transcript,
      action="add"
    )

  expect_equal(
    (res$`Dimension 1`)[1:4],
    c(0.277291, 0.277291, 0.277291, 0.277291),
    tolerance=1e-7
  )

  expect_equal(
    ncol(res),
    9
  )

})

test_that("Get reduced dimensions PCA - no object",{

  res =
    reduce_dimensions(
      ttBulk::counts_mini,
      method="PCA",
      value_column = `read count`,
      elements_column = sample,
      feature_column = transcript,
      action="get"
    )

  expect_equal(
    res$PC1,
    c(0.4993684, 0.4994505, 0.5004880, 0.5006918),
    tolerance=1e-7
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
      value_column = `read count`,
      elements_column = sample,
      feature_column = transcript,
      action="add"
    )

  expect_equal(
    res$PC1[1:4],
    c(0.4993684, 0.4993684, 0.4993684, 0.4993684),
    tolerance=1e-7
  )

  expect_equal(
    ncol(res),
    9
  )

})

test_that("Get rotated dimensions - no object",{

  res.pca =
    reduce_dimensions(
      ttBulk::counts_mini,
      method="PCA",
      value_column = `read count`,
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
    c(-0.08831853,  0.77379931 , 0.61726229 , 0.11145282),
    tolerance=1e-7
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
      value_column = `read count`,
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
    c( -0.08831853, -0.08831853, -0.08831853, -0.08831853),
    tolerance=1e-7
  )

  expect_equal(
    ncol(res),
    11
  )

})

test_that("Aggregate duplicated transcript - no object",{

  res =
    aggregate_duplicates(
      ttBulk::counts_mini,
      sample_column = sample,
      transcript_column = transcript,
      counts_column = `read count`
    )

  expect_equal(
    res$transcript[1:4],
    c("KLHL17" ,     "TRN-GTT12-1", "RAB25"  , "DLG5")
  )

  expect_equal(
    ncol(res),
    8
  )

})

test_that("Drop redundant correlated - no object",{

  res =
    drop_redundant(
      ttBulk::counts_mini,
      method = "correlation",
      value_column = `read count`,
      elements_column = sample,
      feature_column = transcript
    )

  expect_equal(
    nrow(res),
    30
  )

  expect_equal(
    ncol(res),
    7
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

# test_that("Get cell type proportions - no object",{
#
#   expect_error(
#     annotate_cell_type(
#       ttBulk::counts_mini,
#       sample_column = sample,
#       transcript_column = transcript,
#       counts_column = `read count`
#     ),
#     "You have less than 50 genes that overlap the Cibersort signature"
#   )
#
#   res =
#     annotate_cell_type(
#       ttBulk::counts[
#         ttBulk::counts$sample %in%
#           unique(ttBulk::counts_mini$sample),
#         ],
#       sample_column = sample,
#       transcript_column = transcript,
#       counts_column = `read count`,
#       action="get"
#     )
#
#   expect_equal(
#     res$proportion[1:4],
#     c(0.6273726, 0.6004113, 0.5931743, 0.5811784),
#     tolerance=1e-7
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
#     annotate_cell_type(
#       ttBulk::counts_mini,
#       sample_column = sample,
#       transcript_column = transcript,
#       counts_column = `read count`,
#       cores = 2
#     ),
#     "You have less than 50 genes that overlap the Cibersort signature"
#   )
#
#   res =
#     annotate_cell_type(
#       ttBulk::counts[
#         ttBulk::counts$sample %in%
#           unique(ttBulk::counts_mini$sample),
#         ],
#       sample_column = sample,
#       transcript_column = transcript,
#       counts_column = `read count`,
#       action="add",
#       cores = 2
#     )
#
#   expect_equal(
#     res$proportion[1:4],
#     c(0.6273726, 0.2553588, 0.0000000, 0.0000000),
#     tolerance=1e-7
#   )
#
#   expect_equal(
#     ncol(res),
#     7
#   )
#
# })
