





library(tidySummarizedExperiment)
library(tidybulk)

# Sample annotations
sample_annotations <- data.frame(
	sample_id = paste0("Sample", seq(1, 9)),
	factor_of_interest = c(rep("treated", 4), rep("untreated", 5)),
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

se |>
  resolve_complete_confounders_of_non_interest(A, B, C) |>
  distinct(.sample, A, B, C)



