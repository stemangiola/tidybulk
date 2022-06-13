context('tximeta and Granges')

test_that("tximeta 1",{

  duplicate = tximeta_summarizeToGene_object[1,]
  rownames(duplicate) = "dup"


  rbind(duplicate, tximeta_summarizeToGene_object)  %>%
  aggregate_duplicates(.transcript = gene_id)

  tximeta_summarizeToGene_object  %>%
    aggregate_duplicates(.transcript = gene_id)

})



test_that("se no features",{

  # Create dataset
  nrows <- 200; ncols <- 6
  counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
  rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                       IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                       strand=sample(c("+", "-"), 200, TRUE),
                       feature_id=sprintf("ID%03d", 1:200))
  colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3), row.names=LETTERS[1:6])
  se <- SummarizedExperiment(assays=SimpleList(counts=counts), rowRanges=rowRanges, colData=colData)
  se= rbind( se[1,], se)

  se %>%
    aggregate_duplicates(.transcript = feature_id)
})





