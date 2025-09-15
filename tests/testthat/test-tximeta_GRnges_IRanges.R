context('tximeta and Granges')

library(tidybulk)
# Load required data
data(tximeta_summarizeToGene_object)

test_that("tximeta 1",{

  duplicate = tximeta_summarizeToGene_object[1,]
  rownames(duplicate) = "dup"


  SummarizedExperiment::rbind(duplicate, tximeta_summarizeToGene_object)  %>%
  aggregate_duplicates(feature = "gene_id")

  tximeta_summarizeToGene_object  %>%
    aggregate_duplicates(feature = "gene_id")

})





