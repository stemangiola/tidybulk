## ---- echo=FALSE, include=FALSE-----------------------------------------------
library(knitr)
#library(kableExtra)
knitr::opts_chunk$set(cache = TRUE, warning = FALSE,
                      message = FALSE, cache.lazy = FALSE)
#options(width = 120)
options(pillar.min_title_chars = Inf)


library(tibble)
library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)
# library(widyr)
library(rlang)
library(purrr)
library(tidybulk)

my_theme = 	
	theme_bw() +
	theme(
		panel.border = element_blank(),
		axis.line = element_line(),
		panel.grid.major = element_line(size = 0.2),
		panel.grid.minor = element_line(size = 0.1),
		text = element_text(size=12),
		legend.position="bottom",
		aspect.ratio=1,
		strip.background = element_blank(),
		axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
		axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
	)

# counts_mini =
# 	tidybulk::counts %>%
# 	filter(transcript %in% (tidybulk::X_cibersort %>% rownames)) %>%
# 	filter(sample %in% c("SRR1740034", "SRR1740035", "SRR1740058", "SRR1740043", "SRR1740067")) %>%
# 	mutate(condition = ifelse(sample %in% c("SRR1740034", "SRR1740035", "SRR1740058"), TRUE, FALSE))

# se_mini
se_mini = tidybulk:::tidybulk_to_SummarizedExperiment(tidybulk::counts_mini, sample, transcript, count)
se_breast_tcga_mini = tidybulk:::tidybulk_to_SummarizedExperiment( tidybulk::breast_tcga_mini, sample, ens, `count`)
se.cibersort =
	tidybulk:::tidybulk_to_SummarizedExperiment(tidybulk::counts,  sample ,  transcript, count)
se.norm.batch =
	tidybulk:::tidybulk_to_SummarizedExperiment(tidybulk::counts,  sample ,  transcript, count) %>%
	scale_abundance()

## ----adjust se, cache=TRUE----------------------------------------------------
se.norm.batch %>%
  adjust_abundance(
  	~ factor_of_interest + batch,
  	.abundance = count_scaled
  )

## ----cibersort, cache=TRUE----------------------------------------------------
tt.cibersort =
	tt %>%
	deconvolve_cellularity(action="get", cores=2)

tt.cibersort %>% select(sample, contains("cibersort:")) 

## ----plot_cibersort, cache=TRUE-----------------------------------------------
tt.cibersort %>%
	gather(`Cell type inferred`, `proportion`, 5:26) %>%
  distinct(sample, `Cell type`, `Cell type inferred`, proportion) %>%
  ggplot(aes(x=`Cell type inferred`, y=proportion, fill=`Cell type`)) +
  geom_boxplot() +
  facet_wrap(~`Cell type`) +
  my_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), aspect.ratio=1/5)

## ----cibersort se, cache=TRUE-------------------------------------------------

se.cibersort %>% deconvolve_cellularity(cores=2)


## ----cluster, cache=TRUE------------------------------------------------------
tt.norm.cluster = tt.norm %>%
  cluster_elements(.abundance = count_scaled, method="kmeans",	centers = 2 )

tt.norm.cluster

## ----plot_cluster, cache=TRUE-------------------------------------------------
 tt.norm.MDS %>%
  cluster_elements(
  	.abundance = count_scaled,
  	method="kmeans",
  	centers = 2,
  	action="get"
  ) %>%
	ggplot(aes(x=`Dim1`, y=`Dim2`, color=`cluster kmeans`)) +
  geom_point() +
  my_theme

## ----cluster se, cache=TRUE---------------------------------------------------
se.norm %>%
  cluster_elements(.abundance = count_scaled, method="kmeans",	centers = 2 )


## ----SNN, cache=TRUE----------------------------------------------------------
tt.norm.SNN =	tt.norm.tSNE %>%	cluster_elements(.abundance= count_scaled, method = "SNN")

tt.norm.SNN %>%
	pivot_sample()

tt.norm.SNN %>%
	select(contains("tSNE", ignore.case = FALSE), `cluster SNN`, sample, Call) %>%
	gather(source, Call, c("cluster SNN", "Call")) %>%
	distinct() %>%
	ggplot(aes(x = `tSNE1`, y = `tSNE2`, color=Call)) + geom_point() + facet_grid(~source) + my_theme


# Do differential transcription between clusters
tt.norm.SNN %>%
	mutate(factor_of_interest = `cluster SNN` == 3) %>%
	test_differential_abundance(
    ~ factor_of_interest,
    action="only"
   )

## ----SNN se, cache=TRUE-------------------------------------------------------
se.norm.tSNE %>%	cluster_elements(.abundance= count_scaled, method = "SNN")

## ----drop, cache=TRUE---------------------------------------------------------
tt.norm.non_redundant = tt.norm.MDS %>%  remove_redundancy(	method = "correlation" )

## ----plot_drop, cache=TRUE----------------------------------------------------
tt.norm.non_redundant %>%
	pivot_sample() %>%
	ggplot(aes(x=`Dim1`, y=`Dim2`, color=`Cell type`)) +
  geom_point() +
  my_theme


## ----drop se, cache=TRUE------------------------------------------------------
se.norm.MDS %>%  remove_redundancy(	method = "correlation" )

## ----drop2, cache=TRUE--------------------------------------------------------
tt.norm.non_redundant =
	tt.norm.MDS %>%
  remove_redundancy(
  	method = "reduced_dimensions",
  	.element = sample,
  	.feature = transcript,
  	Dim_a_column = `Dim1`,
  	Dim_b_column = `Dim2`
  )

## ----plot_drop2, cache=TRUE---------------------------------------------------
tt.norm.non_redundant %>%
	pivot_sample() %>%
	ggplot(aes(x=`Dim1`, y=`Dim2`, color=`Cell type`)) +
  geom_point() +
  my_theme


## ----drop2 se, cache=TRUE-----------------------------------------------------

se.norm.MDS %>%
remove_redundancy(
	method = "reduced_dimensions",
	.element = sample,
	.feature = transcript,
	Dim_a_column = `Dim1`,
	Dim_b_column = `Dim2`
)

## ----eval=FALSE---------------------------------------------------------------
#  counts = tidybulk_SAM_BAM(
#  	file_names,
#  	genome = "hg38",
#  	isPairedEnd = TRUE,
#  	requireBothEndsMapped = TRUE,
#  	checkFragLength = FALSE,
#  	useMetaFeatures = TRUE
#  )

## ----ensembl, cache=TRUE------------------------------------------------------
counts_ensembl %>% ensembl_to_symbol(ens)

## ---- cache=TRUE--------------------------------------------------------------
  tt.norm

## ---- cache=TRUE--------------------------------------------------------------
  tt.norm %>%
    reduce_dimensions(
    	.abundance = count_scaled,
    	method="MDS" ,
    	.element = sample,
    	.feature = transcript,
    	.dims = 3,
    	action="add"
    )

## ---- cache=TRUE--------------------------------------------------------------
  tt.norm %>%
    reduce_dimensions(
    	.abundance = count_scaled,
    	method="MDS" ,
    	.element = sample,
    	.feature = transcript,
    	.dims = 3,
    	action="get"
    )

## ---- cache=TRUE--------------------------------------------------------------
  tt.norm %>%
    reduce_dimensions(
    	.abundance = count_scaled,
    	method="MDS" ,
    	.element = sample,
    	.feature = transcript,
    	.dims = 3,
    	action="only"
    )

## -----------------------------------------------------------------------------
sessionInfo()

