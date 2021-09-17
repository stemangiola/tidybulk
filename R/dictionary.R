# Dictionary

scaled_string = "_scaled"
adjusted_string = "_adjusted"

warning_for_scaling_with_few_genes = "tidybulk says: There are < 100 features/genes that are present in all you samples. Because edgeR::calcNormFactors does not allow NAs, the scaling is performed on that limited set of features.genes. The scaling could not be accurate, it is adivasble to perform impute_missing_abundance() before scaling. It is possible to filter the imputed counts after scaling."
