\name{NEWS}
\title{News for Package \pkg{tidybulk}}

\section{Changes in version 1.2.0, Bioconductor 3.12 Release}{
\itemize{
		\item Make gene filtering functionality `identify_abundance` explicit, a warning will be given if this has not been performed before the majority of workflow steps (e.g. `test_differential_abundance`).
    \item Add DESeq2 and limma-voom to the methods for `test_differential_abundance` (method="DESeq2").
    \item Add other cell-type signature to `deconvolve_cellularity`.
    \item Add differential cellularity analyses `test_differential_cellularity`.
    \item Add gene descrption annotation utility `describe_transcript`.
    \item Add `nest` functionality for functional-coding transcriptomic analyses.
    \item Add gene overrepresentation functionality `test_gene_overrepresentation`.
    \item Add github website.
    \item Several bug fixes.
}}

