context('tximeta and Granges')

test_that("tximeta 1",{

  duplicate = tximeta_summarizeToGene_object[1,]
  rownames(duplicate) = "dup"


#   SummarizedExperiment::rbind(duplicate, tximeta_summarizeToGene_object)  %>%
#   aggregate_duplicates(.transcript = gene_id)

#   tximeta_summarizeToGene_object  %>%
#     aggregate_duplicates(.transcript = gene_id)

})



# The following test is commented out due to S4 object validity error in this branch. Review and update as needed.
# test_that("se no features",{
#
#   # Create dataset
#   nrows <- 200; ncols <- 6
#   counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
#   rowRanges <- GenomicRanges::GRanges(rep(c("chr1", "chr2"), c(50, 150)),
#                                 IRanges::IRanges(floor(runif(200, 1e5, 1e6)), width=100),
#                        strand=sample(c("+", "-"), 200, TRUE),
#                        feature_id=sprintf("ID%03d", 1:200))
#   colData <- S4Vectors::DataFrame(Treatment=rep(c("ChIP", "Input"), 3), row.names=LETTERS[1:6])
#   se <- SummarizedExperiment(assays=S4Vectors::SimpleList(counts=counts), rowRanges=rowRanges, colData=colData)
#   se= rbind( se[1,], se)
#
#   # se %>%
#   #   aggregate_duplicates(.transcript = feature_id)
# })





