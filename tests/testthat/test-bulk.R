context('Bulk functions')

test_that("Test data frame",{ expect_equal( ncol(tidyTranscriptomics::counts_mini), 6 ) })

test_that("Creating tt object from tibble, number of parameters",{

  expect_equal(

    length(
      attr(
        create_tt_from_tibble_bulk(
          tidyTranscriptomics::counts_mini,
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

test_that("Getting normalised counts - no object",{

  res =
    get_normalised_counts_bulk(
      tidyTranscriptomics::counts_mini,
      sample_column = sample,
      transcript_column = transcript,
      counts_column = `read count`
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
    add_normalised_counts_bulk(
      tidyTranscriptomics::counts_mini,
      sample_column = sample,
      transcript_column = transcript,
      counts_column = `read count`
    )

  expect_equal(
    unique(res$multiplier),
    c(1.0057758, 0.8828939, 0.8369487, 1.1423585),
    tolerance=1e-7
  )

  expect_equal(
    ncol(res),
    10
  )

})

test_that("Get differential trancript abundance - no object",{

  res =
    get_differential_transcript_abundance_bulk(
      tidyTranscriptomics::counts_mini,
      ~ condition,
      sample_column = sample,
      transcript_column = transcript,
      counts_column = `read count`
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

test_that("Get cluster lables based on Kmeans - no object",{

  res =
    get_clusters_kmeans_bulk(
      tidyTranscriptomics::counts_mini,
      value_column = `read count`,
      elements_column = sample,
      feature_column = transcript,
      number_of_clusters = 2
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
    add_clusters_kmeans_bulk(
      tidyTranscriptomics::counts_mini,
      value_column = `read count`,
      elements_column = sample,
      feature_column = transcript,
      number_of_clusters = 2
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
    get_reduced_dimensions_MDS_bulk(
      tidyTranscriptomics::counts_mini,
      value_column = `read count`,
      elements_column = sample,
      feature_column = transcript
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
    add_reduced_dimensions_MDS_bulk(
      tidyTranscriptomics::counts_mini,
      value_column = `read count`,
      elements_column = sample,
      feature_column = transcript
    )

  expect_equal(
    (res$`Dimension 1`)[1:4],
    c(0.277291, 0.277291, 0.277291, 0.277291),
    tolerance=1e-7
  )

  expect_equal(
    ncol(res),
    8
  )

})

test_that("Get reduced dimensions PCA - no object",{

  res =
    get_reduced_dimensions_PCA_bulk(
      tidyTranscriptomics::counts_mini,
      value_column = `read count`,
      elements_column = sample,
      feature_column = transcript
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
    add_reduced_dimensions_PCA_bulk(
      tidyTranscriptomics::counts_mini,
      value_column = `read count`,
      elements_column = sample,
      feature_column = transcript
    )

  expect_equal(
    res$PC1[1:4],
    c(0.4993684, 0.4993684, 0.4993684, 0.4993684),
    tolerance=1e-7
  )

  expect_equal(
    ncol(res),
    8
  )

})

test_that("Get rotated dimensions - no object",{

  res.pca =
    add_reduced_dimensions_PCA_bulk(
      tidyTranscriptomics::counts_mini,
      value_column = `read count`,
      elements_column = sample,
      feature_column = transcript
    )

  res =
    get_rotated_dimensions(
      res.pca,
      dimension_1_column = PC1,
      dimension_2_column = PC2,
      rotation_degrees = 45,
      elements_column = sample
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
    add_reduced_dimensions_PCA_bulk(
      tidyTranscriptomics::counts_mini,
      value_column = `read count`,
      elements_column = sample,
      feature_column = transcript
    )

  res =
    add_rotated_dimensions(
      res.pca,
      dimension_1_column = PC1,
      dimension_2_column = PC2,
      rotation_degrees = 45,
      elements_column = sample
    )

  expect_equal(
    res$`PC1 rotated 45`[1:4],
    c( -0.08831853, -0.08831853, -0.08831853, -0.08831853),
    tolerance=1e-7
  )

  expect_equal(
    ncol(res),
    10
  )

})

test_that("Aggregate duplicated transcript - no object",{

  res =
    aggregate_duplicated_transcripts_bulk(
      tidyTranscriptomics::counts_mini,
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
    7
  )

})

test_that("Drop redundant correlated - no object",{

  res =
    drop_redundant_elements_through_correlation(
      tidyTranscriptomics::counts_mini,
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
    6
  )

})


test_that("Get symbol from ensambl - no object",{

  res =
    get_symbol_from_ensembl(
      tidyTranscriptomics::counts_ensembl,
      ensembl_transcript_column = ens
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
    add_symbol_from_ensembl(
      tidyTranscriptomics::counts_ensembl,
      ensembl_transcript_column = ens
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

  expect_error(
    get_cell_type_proportions(
      tidyTranscriptomics::counts_mini,
      sample_column = sample,
      transcript_column = transcript,
      counts_column = `read count`
    ),
    "You have less than 50 genes that overlap the Cibersort signature"
  )

  res =
    get_cell_type_proportions(
      tidyTranscriptomics::counts[
        tidyTranscriptomics::counts$sample %in%
          unique(tidyTranscriptomics::counts_mini$sample),
        ],
      sample_column = sample,
      transcript_column = transcript,
      counts_column = `read count`
    )

  expect_equal(
    res$proportion[1:4],
    c(0.6273726, 0.6004113, 0.5931743, 0.5811784),
    tolerance=1e-7
  )

  expect_equal(
    ncol(res),
    3
  )

})

test_that("Add cell type proportions - no object",{

  expect_error(
    add_cell_type_proportions(
      tidyTranscriptomics::counts_mini,
      sample_column = sample,
      transcript_column = transcript,
      counts_column = `read count`
    ),
    "You have less than 50 genes that overlap the Cibersort signature"
  )

  res =
    add_cell_type_proportions(
      tidyTranscriptomics::counts[
        tidyTranscriptomics::counts$sample %in%
          unique(tidyTranscriptomics::counts_mini$sample),
        ],
      sample_column = sample,
      transcript_column = transcript,
      counts_column = `read count`
    )

  expect_equal(
    res$proportion[1:4],
    c(0.6273726, 0.2553588, 0.0000000, 0.0000000),
    tolerance=1e-7
  )

  expect_equal(
    ncol(res),
    7
  )

})
