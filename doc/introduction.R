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

## -----------------------------------------------------------------------------
tt = counts_mini %>% tidybulk(sample, transcript, count)

## ----aggregate, cache=TRUE----------------------------------------------------
tt.aggr =  tt %>% aggregate_duplicates( 	aggregation_function = sum )

tt.aggr

## ----aggregate se, cache=TRUE-------------------------------------------------
se.aggr =  se_mini %>% aggregate_duplicates( 	aggregation_function = sum )

se.aggr

## ----normalise, cache=TRUE----------------------------------------------------
tt.norm =  tt.aggr %>% scale_abundance(method="TMM")

tt.norm %>% select(`count`, count_scaled, lowly_abundant, everything())

## ----plot_normalise, cache=TRUE-----------------------------------------------
tt.norm %>%
	ggplot(aes(count_scaled + 1, group=sample, color=`Cell type`)) +
	geom_density() +
	scale_x_log10() +
	my_theme

## ----normalise se, cache=TRUE-------------------------------------------------
se.norm =  se.aggr %>% scale_abundance(method="TMM")

se.norm

## ----filter variable, cache=TRUE----------------------------------------------
tt.norm.variable = tt.norm %>% keep_variable()

## ----mds, cache=TRUE----------------------------------------------------------
tt.norm.MDS =  tt.norm %>% reduce_dimensions(.abundance = count_scaled, method="MDS", .dims = 3)

tt.norm.MDS %>% select(sample, contains("Dim"), `Cell type`, time ) %>% distinct()

## ----plot_mds, cache=TRUE, eval=FALSE-----------------------------------------
#  tt.norm.MDS %>%
#  	select(contains("Dim"), sample, `Cell type`) %>%
#    distinct() %>%
#    GGally::ggpairs(columns = 1:3, ggplot2::aes(colour=`Cell type`))

## ----mds se, cache=TRUE-------------------------------------------------------
se.norm.MDS =  se.norm %>% reduce_dimensions(.abundance = count_scaled, method="MDS", .dims = 3)

se.norm.MDS

## ----pca, cache=TRUE----------------------------------------------------------
tt.norm.PCA = tt.norm %>% reduce_dimensions(.abundance = count_scaled, method="PCA" ,  .dims = 3)

tt.norm.PCA %>% select(sample, contains("PC"), `Cell type`, time ) %>% distinct()

## ----plot_pca, cache=TRUE, eval=FALSE-----------------------------------------
#  tt.norm.PCA %>%
#  	select(contains("PC"), sample, `Cell type`) %>%
#    distinct() %>%
#    GGally::ggpairs(columns = 1:3, ggplot2::aes(colour=`Cell type`))

## ----pca se, cache=TRUE-------------------------------------------------------
se.norm.PCA = se.norm %>% reduce_dimensions(.abundance = count_scaled, method="PCA" ,  .dims = 3)

se.norm.PCA

## ---- echo=FALSE, include=FALSE-----------------------------------------------
tt_tcga_breast =
	tidybulk::breast_tcga_mini %>%
	tidybulk(sample, ens, `count`)

## ----tsne, cache=TRUE---------------------------------------------------------
tt.norm.tSNE =
	tt_tcga_breast %>%
	reduce_dimensions(
		.abundance = count_scaled,
		method = "tSNE",
		top = 500,
		perplexity=10,
		pca_scale =TRUE
	)

tt.norm.tSNE %>%
	select(contains("tSNE", ignore.case = FALSE), sample, Call) %>%
	distinct()

tt.norm.tSNE %>%
	pivot_sample() %>%
	ggplot(aes(x = `tSNE1`, y = `tSNE2`, color=Call)) + geom_point() + my_theme

## ----tsne se, cache=TRUE------------------------------------------------------
se.norm.tSNE =
	se_breast_tcga_mini %>%
	reduce_dimensions(
		.abundance = count_scaled,
		method = "tSNE",
		top = 500,
		perplexity=10,
		pca_scale =TRUE
	)
se.norm.tSNE

## ----rotate, cache=TRUE-------------------------------------------------------
tt.norm.MDS.rotated =
  tt.norm.MDS %>%
	rotate_dimensions(`Dim1`, `Dim2`, rotation_degrees = 45, .element = sample)

## ----plot_rotate_1, cache=TRUE------------------------------------------------
tt.norm.MDS.rotated %>%
	pivot_sample() %>%
	ggplot(aes(x=`Dim1`, y=`Dim2`, color=`Cell type` )) +
  geom_point() +
  my_theme

## ----plot_rotate_2, cache=TRUE------------------------------------------------
tt.norm.MDS.rotated %>%
	pivot_sample() %>%
	ggplot(aes(x=`Dim1 rotated 45`, y=`Dim2 rotated 45`, color=`Cell type` )) +
  geom_point() +
  my_theme

## ----rotate se, cache=TRUE----------------------------------------------------
se.norm.MDS %>%
rotate_dimensions(`Dim1`, `Dim2`, rotation_degrees = 45, .element = sample)

## ----de, cache=TRUE-----------------------------------------------------------
tt %>%	test_differential_abundance(  ~ condition,  action="only")

## ----de se, cache=TRUE--------------------------------------------------------
se_mini %>%	test_differential_abundance(  ~ condition)

## ---- echo=FALSE, include=FALSE-----------------------------------------------
tt.norm.batch =
	tt.norm %>%

	  # Add fake batch and factor of interest
	  left_join(
	  	(.) %>%
	  		distinct(sample) %>%
	  		mutate(batch = c(0,1,0,1,1))
	  ) %>%
	 	mutate(factor_of_interest = `Cell type` == "b_cell")


## ----adjust, cache=TRUE-------------------------------------------------------
tt.norm.adj =
	tt.norm.batch %>%
	  adjust_abundance(
	  	~ factor_of_interest + batch,
	  	.abundance = count_scaled,
	  	action = "only"
	  )

tt.norm.adj

## ----adjust se, cache=TRUE----------------------------------------------------
se.norm.batch %>%
  adjust_abundance(
  	~ factor_of_interest + batch,
  	.abundance = count_scaled
  )

## ----cibersort, cache=TRUE----------------------------------------------------
tt.cibersort =
	tt %>%
	deconvolve_cellularity(action="get", cores=1)

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

se.cibersort %>% deconvolve_cellularity(cores=1)


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

