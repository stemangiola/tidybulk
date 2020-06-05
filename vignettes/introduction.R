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
library(widyr)
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

## ----pca, cache=TRUE----------------------------------------------------------
tt.norm.PCA = tt.norm %>% reduce_dimensions(.abundance = count_scaled, method="PCA" ,  .dims = 3)

tt.norm.PCA %>% select(sample, contains("PC"), `Cell type`, time ) %>% distinct()

## ----plot_pca, cache=TRUE-----------------------------------------------------
tt.norm.PCA %>%
	select(contains("PC"), sample, `Cell type`) %>%
  distinct() %>%
  GGally::ggpairs(columns = 1:3, ggplot2::aes(colour=`Cell type`))

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


