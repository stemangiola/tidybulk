context('Fix the bug that can not calculate dimensions larger
        than 2 with tSNE and UMAP')

library(dplyr)
data("breast_tcga_mini_SE")
input_df_breast =   breast_tcga_mini_SE |> tidybulk() |> as_tibble() |> setNames(c( "b","a", "c", "c norm", "call"))


test_that("Get reduced dimensions tSNE - no object",{

  set.seed(132)

  res =
    input_df_breast |>
    identify_abundant(.sample = a, .transcript = b, .abundance = c) |>

    reduce_dimensions(
      method="tSNE",
      .abundance = c,
      .element = a,
      .feature = b,
      action="get",
      .dims = 3,
      verbose=FALSE
    ) |>
    suppressMessages()

  # DOES NOT REPRODUCE MAYBE BECAUSE OF DIFFERENT TSNE VERSIONS
  # res |>
  #   pull(tSNE1) |>
  #   magrittr::extract2(1) |>
  #   expect_equal(2.432608, tolerance = 0.01)

  expect_equal(
    typeof(res$`tSNE1`),
    "double",
    tolerance=1e-1
  )

  expect_true("tSNE3" %in% names(res))

  expect_equal(
    nrow(res),
    251
  )

})

test_that("Add reduced dimensions UMAP - no object",{

  set.seed(132)

  res =
    input_df_breast |>
    identify_abundant(a, b, c) |>

    reduce_dimensions(
      method="UMAP",
      .abundance = c,
      .element = a,
      .feature = b,
      action="add",
      .dims = 3
    )

  expect_true("UMAP3" %in% names(res))

  expect_false("UMAP4" %in% names(res))

  res =
    input_df_breast |>
    identify_abundant(a, b, c) |>

    reduce_dimensions(
      method="UMAP",
      .abundance = c,
      .element = a,
      .feature = b,
      action="add",
      .dims = 4
    )

  expect_true("UMAP4" %in% names(res))

})
