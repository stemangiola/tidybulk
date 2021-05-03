\name{NEWS}
\title{News for Package \pkg{tidybulk}}

\section{Changes in version 1.2.0, Bioconductor 3.12 Release}{
\itemize{
    \item Make gene filtering functionality `identify_abundance` explicit, a warning will be given if this has not been performed before the majority of workflow steps (e.g. `test_differential_abundance`).
    \item Add Automatic bibliography `get_bibliography`.
    \item Add DESeq2 and limma-voom to the methods for `test_differential_abundance` (method="DESeq2").
    \item Add prefix to test_differential_abundance for multi-methods analyses.
    \item Add other cell-type signature to `deconvolve_cellularity`.
    \item Add differential cellularity analyses `test_differential_cellularity`.
    \item Add gene descrption annotation utility `describe_transcript`.
    \item Add `nest` functionality for functional-coding transcriptomic analyses.
    \item Add gene overrepresentation functionality `test_gene_overrepresentation`.
    \item Add github website.
    \item Seep up data frame vadidation.
    \item Several bug fixes.
}}

\section{Changes in version 1.3.2, Bioconductor 3.13 Release}{
\itemize{
    \item Tidybulk now operates natively with SummarizedExperment data container, in a seamless way thanks to tidySummarisedExperiment 10.18129/B9.bioc.tidySummarizedExperiment     
    \item Added robust edgeR as it outperforms many other methods as shown here doi.org/10.1093/nargab/lqab005
    \item Added test stratifiction cellularity, to easily calculate Kaplan-Meier curves
    \item Production of SummarizedExperiment from BAM or SAM files
    \item Added treat method to edgeR and voom differential transcription tests doi.org/10.1093/bioinformatics/btp053
    \item Added the method as_SummarizedExperiment
    \item Vastly improved test_gene_enrichment
    \item Added test_gene_rank, based on GSEA
    \item Several bug fixes.
}}
