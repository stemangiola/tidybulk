#' Global variable declarations to avoid R CMD check warnings
#'
#' This script declares global variables used throughout the package to avoid R CMD check NOTES.
#'
#' @noRd
if(getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "transcript", "read count", "n", "m", ".", ".feature", ".abundance_scaled", 
    "tt_columns", "seurat_clusters", "cluster", "tagwise.dispersion", "Component", 
    "Component value", "sdev", "name", "value", "x", "Y", "x",
    "item1", "rowid", "n1", "n2", "se_data", "sample 1", "sample 2",
    ".proportion", ".cell_type", "cell_type_proportions", "surv_test",
    "X_cibersort", "P-value", "Correlation", "RMSE", "cell_type",
    "pathway", "my_do_test", "sample_idx", "(Intercept)",
    "group_name", "group", "EPIC", "rotation",
    "variable", "nulldist", "my_n", "gs_collection", "test", "geneID", 
    "ncbi_gene", "fit", "FDR", "constrast", "adj.P.Val", "group_id", 
    "parameter", "CI", "lower", "upper", "term", "sample a", "sample b",
    "abundance", "element", "feature", "rc", "correlation",
    "buildCustomIdx", "buildIdx", "egsea", "normalize.quantiles"
  ))
} 