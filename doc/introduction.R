## ---- echo=FALSE, include=FALSE-----------------------------------------------
library(knitr)
knitr::opts_chunk$set(cache = TRUE, warning = FALSE,
                      message = FALSE, cache.lazy = FALSE)

library(dplyr)
library(tidyr)
library(tibble)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(tidybulk)
library(tidySummarizedExperiment)

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


## ----eval=FALSE---------------------------------------------------------------
#  BiocManager::install("tidybulk")

## ---- eval=FALSE--------------------------------------------------------------
#  devtools::install_github("stemangiola/tidybulk")

## -----------------------------------------------------------------------------
counts_SE

## -----------------------------------------------------------------------------
class(counts_SE)

## ----eval=FALSE---------------------------------------------------------------
#  counts_SE %>%	get_bibliography()

## ----aggregate, message=FALSE, warning=FALSE, results='hide', class.source='yellow'----
counts_SE.aggr = counts_SE %>% aggregate_duplicates()

## ----normalise----------------------------------------------------------------
counts_SE.norm = counts_SE.aggr %>% identify_abundant(factor_of_interest = condition) %>% scale_abundance()

## ---- include=FALSE-----------------------------------------------------------
counts_SE.norm %>% select(`count`, count_scaled, .abundant, everything())

## ----plot_normalise-----------------------------------------------------------
counts_SE.norm %>%
	ggplot(aes(count_scaled + 1, group=sample, color=`Cell.type`)) +
	geom_density() +
	scale_x_log10() +
	my_theme

## ----filter variable----------------------------------------------------------
counts_SE.norm.variable = counts_SE.norm %>% keep_variable()

## ----filter variable long, eval=FALSE-----------------------------------------
#  library(edgeR)
#  
#  x = norm_counts.table
#  
#  s <- rowMeans((x-rowMeans(x))^2)
#  o <- order(s,decreasing=TRUE)
#  x <- x[o[1L:top],,drop=FALSE]
#  
#  norm_counts.table = norm_counts.table[rownames(x)]
#  
#  norm_counts.table$cell_type = tidybulk::counts[
#  	match(
#  		tidybulk::counts$sample,
#  		rownames(norm_counts.table)
#  	),
#  	"Cell.type"
#  ]

## ----mds----------------------------------------------------------------------
counts_SE.norm.MDS =
  counts_SE.norm %>%
  reduce_dimensions(method="MDS", .dims = 6)


## ---- eval = FALSE------------------------------------------------------------
#  library(limma)
#  
#  count_m_log = log(count_m + 1)
#  cmds = limma::plotMDS(ndim = .dims, plot = FALSE)
#  
#  cmds = cmds %$%	
#  	cmdscale.out %>%
#  	setNames(sprintf("Dim%s", 1:6))
#  
#  cmds$cell_type = tidybulk::counts[
#  	match(tidybulk::counts$sample, rownames(cmds)),
#  	"Cell.type"
#  ]

## ----plot_mds-----------------------------------------------------------------
counts_SE.norm.MDS %>% pivot_sample()  %>% select(contains("Dim"), everything())

counts_SE.norm.MDS %>%
	pivot_sample() %>%
  GGally::ggpairs(columns = 10:15, ggplot2::aes(colour=`Cell.type`))



## ----pca, message=FALSE, warning=FALSE, results='hide'------------------------
counts_SE.norm.PCA =
  counts_SE.norm %>%
  reduce_dimensions(method="PCA", .dims = 6)

## ----eval=FALSE---------------------------------------------------------------
#  count_m_log = log(count_m + 1)
#  pc = count_m_log %>% prcomp(scale = TRUE)
#  variance = pc$sdev^2
#  variance = (variance / sum(variance))[1:6]
#  pc$cell_type = counts[
#  	match(counts$sample, rownames(pc)),
#  	"Cell.type"
#  ]

## ----plot_pca-----------------------------------------------------------------

counts_SE.norm.PCA %>% pivot_sample() %>% select(contains("PC"), everything())

counts_SE.norm.PCA %>%
	 pivot_sample() %>%
  GGally::ggpairs(columns = 11:13, ggplot2::aes(colour=`Cell.type`))

## ----tsne, message=FALSE, warning=FALSE, results='hide'-----------------------
counts_SE.norm.tSNE =
	breast_tcga_mini_SE %>%
	identify_abundant() %>%
	reduce_dimensions(
		method = "tSNE",
		perplexity=10,
		pca_scale =TRUE
	)

## ---- eval=FALSE--------------------------------------------------------------
#  count_m_log = log(count_m + 1)
#  
#  tsne = Rtsne::Rtsne(
#  	t(count_m_log),
#  	perplexity=10,
#  		pca_scale =TRUE
#  )$Y
#  tsne$cell_type = tidybulk::counts[
#  	match(tidybulk::counts$sample, rownames(tsne)),
#  	"Cell.type"
#  ]

## -----------------------------------------------------------------------------
counts_SE.norm.tSNE %>%
	pivot_sample() %>%
	select(contains("tSNE"), everything()) 

counts_SE.norm.tSNE %>%
	pivot_sample() %>%
	ggplot(aes(x = `tSNE1`, y = `tSNE2`, color=Call)) + geom_point() + my_theme

## ----rotate-------------------------------------------------------------------
counts_SE.norm.MDS.rotated =
  counts_SE.norm.MDS %>%
	rotate_dimensions(`Dim1`, `Dim2`, rotation_degrees = 45, action="get")

## ---- eval=FALSE--------------------------------------------------------------
#  rotation = function(m, d) {
#  	r = d * pi / 180
#  	((bind_rows(
#  		c(`1` = cos(r), `2` = -sin(r)),
#  		c(`1` = sin(r), `2` = cos(r))
#  	) %>% as_matrix) %*% m)
#  }
#  mds_r = pca %>% rotation(rotation_degrees)
#  mds_r$cell_type = counts[
#  	match(counts$sample, rownames(mds_r)),
#  	"Cell.type"
#  ]

## ----plot_rotate_1------------------------------------------------------------
counts_SE.norm.MDS.rotated %>%
	ggplot(aes(x=`Dim1`, y=`Dim2`, color=`Cell.type` )) +
  geom_point() +
  my_theme

## ----plot_rotate_2------------------------------------------------------------
counts_SE.norm.MDS.rotated %>%
	pivot_sample() %>%
	ggplot(aes(x=`Dim1_rotated_45`, y=`Dim2_rotated_45`, color=`Cell.type` )) +
  geom_point() +
  my_theme

## ----de, message=FALSE, warning=FALSE, results='hide'-------------------------
counts_SE.de =
	counts_SE %>%
	test_differential_abundance( ~ condition, action="get")
counts_SE.de

## ---- eval=FALSE--------------------------------------------------------------
#  library(edgeR)
#  
#  dgList <- DGEList(counts=counts_m,group=group)
#  keep <- filterByExpr(dgList)
#  dgList <- dgList[keep,,keep.lib.sizes=FALSE]
#  dgList <- calcNormFactors(dgList)
#  design <- model.matrix(~group)
#  dgList <- estimateDisp(dgList,design)
#  fit <- glmQLFit(dgList,design)
#  qlf <- glmQLFTest(fit,coef=2)
#  topTags(qlf, n=Inf)

## ----de contrast, message=FALSE, warning=FALSE, results='hide', eval=FALSE----
#  counts_SE.de =
#  	counts_SE %>%
#  	identify_abundant(factor_of_interest = condition) %>%
#  	test_differential_abundance(
#  		~ 0 + condition,
#  		.contrasts = c( "conditionTRUE - conditionFALSE"),
#  		action="get"
#  	)

## ----adjust, message=FALSE, warning=FALSE, results='hide'---------------------
counts_SE.norm.adj =
	counts_SE.norm %>% adjust_abundance(	~ factor_of_interest + batch)


## ---- eval=FALSE--------------------------------------------------------------
#  library(sva)
#  
#  count_m_log = log(count_m + 1)
#  
#  design =
#  		model.matrix(
#  			object = ~ factor_of_interest + batch,
#  			data = annotation
#  		)
#  
#  count_m_log.sva =
#  	ComBat(
#  			batch =	design[,2],
#  			mod = design,
#  			...
#  		)
#  
#  count_m_log.sva = ceiling(exp(count_m_log.sva) -1)
#  count_m_log.sva$cell_type = counts[
#  	match(counts$sample, rownames(count_m_log.sva)),
#  	"Cell.type"
#  ]
#  

## ----cibersort----------------------------------------------------------------
counts_SE.cibersort =
	counts_SE %>%
	deconvolve_cellularity(action="get", cores=1, prefix = "cibersort__") 


## ---- eval=FALSE--------------------------------------------------------------
#  
#  source(‘CIBERSORT.R’)
#  count_m %>% write.table("mixture_file.txt")
#  results <- CIBERSORT(
#  	"sig_matrix_file.txt",
#  	"mixture_file.txt",
#  	perm=100, QN=TRUE
#  )
#  results$cell_type = tidybulk::counts[
#  	match(tidybulk::counts$sample, rownames(results)),
#  	"Cell.type"
#  ]
#  

## ----plot_cibersort-----------------------------------------------------------
counts_SE.cibersort %>%
	pivot_longer(
		names_to= "Cell_type_inferred", 
		values_to = "proportion", 
		names_prefix ="cibersort__", 
		cols=contains("cibersort__")
	) %>%
  ggplot(aes(x=`Cell_type_inferred`, y=proportion, fill=`Cell.type`)) +
  geom_boxplot() +
  facet_wrap(~`Cell.type`) +
  my_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), aspect.ratio=1/5)

## ----DC-----------------------------------------------------------------------

	counts_SE %>%
	test_differential_cellularity(. ~ condition )


## ----DC_censored, eval=FALSE--------------------------------------------------
#  
#  	counts_SE %>%
#  	test_differential_cellularity(survival::Surv(time, dead) ~ .)
#  

## ----cluster------------------------------------------------------------------
counts_SE.norm.cluster = counts_SE.norm.MDS %>%
  cluster_elements(method="kmeans",	centers = 2, action="get" )

## ---- eval=FALSE--------------------------------------------------------------
#  count_m_log = log(count_m + 1)
#  
#  k = kmeans(count_m_log, iter.max = 1000, ...)
#  cluster = k$cluster
#  
#  cluster$cell_type = tidybulk::counts[
#  	match(tidybulk::counts$sample, rownames(cluster)),
#  	c("Cell.type", "Dim1", "Dim2")
#  ]
#  

## ----plot_cluster-------------------------------------------------------------
 counts_SE.norm.cluster %>%
	ggplot(aes(x=`Dim1`, y=`Dim2`, color=`cluster_kmeans`)) +
  geom_point() +
  my_theme

## ----SNN, message=FALSE, warning=FALSE, results='hide'------------------------
counts_SE.norm.SNN =
	counts_SE.norm.tSNE %>%
	cluster_elements(method = "SNN")

## ---- eval=FALSE--------------------------------------------------------------
#  library(Seurat)
#  
#  snn = CreateSeuratObject(count_m)
#  snn = ScaleData(
#  	snn, display.progress = TRUE,
#  	num.cores=4, do.par = TRUE
#  )
#  snn = FindVariableFeatures(snn, selection.method = "vst")
#  snn = FindVariableFeatures(snn, selection.method = "vst")
#  snn = RunPCA(snn, npcs = 30)
#  snn = FindNeighbors(snn)
#  snn = FindClusters(snn, method = "igraph", ...)
#  snn = snn[["seurat_clusters"]]
#  
#  snn$cell_type = tidybulk::counts[
#  	match(tidybulk::counts$sample, rownames(snn)),
#  	c("Cell.type", "Dim1", "Dim2")
#  ]
#  

## ----SNN_plot-----------------------------------------------------------------
counts_SE.norm.SNN %>%
	pivot_sample() %>%
	select(contains("tSNE"), everything()) 

counts_SE.norm.SNN %>%
	pivot_sample() %>%
	gather(source, Call, c("cluster_SNN", "Call")) %>%
	distinct() %>%
	ggplot(aes(x = `tSNE1`, y = `tSNE2`, color=Call)) + geom_point() + facet_grid(~source) + my_theme


# Do differential transcription between clusters
counts_SE.norm.SNN %>%
	mutate(factor_of_interest = `cluster_SNN` == 3) %>%
	test_differential_abundance(
    ~ factor_of_interest,
    action="get"
   )

## ----drop---------------------------------------------------------------------
counts_SE.norm.non_redundant =
	counts_SE.norm.MDS %>%
  remove_redundancy(	method = "correlation" )

## ---- eval=FALSE--------------------------------------------------------------
#  library(widyr)
#  
#  .data.correlated =
#  	pairwise_cor(
#  		counts,
#  		sample,
#  		transcript,
#  		rc,
#  		sort = TRUE,
#  		diag = FALSE,
#  		upper = FALSE
#  	) %>%
#  	filter(correlation > correlation_threshold) %>%
#  	distinct(item1) %>%
#  	rename(!!.element := item1)
#  
#  # Return non redudant data frame
#  counts %>% anti_join(.data.correlated) %>%
#  	spread(sample, rc, - transcript) %>%
#  	left_join(annotation)
#  
#  
#  

## ----plot_drop----------------------------------------------------------------
counts_SE.norm.non_redundant %>%
	pivot_sample() %>%
	ggplot(aes(x=`Dim1`, y=`Dim2`, color=`Cell.type`)) +
  geom_point() +
  my_theme


## ----drop2--------------------------------------------------------------------
counts_SE.norm.non_redundant =
	counts_SE.norm.MDS %>%
  remove_redundancy(
  	method = "reduced_dimensions",
  	Dim_a_column = `Dim1`,
  	Dim_b_column = `Dim2`
  )

## ----plot_drop2---------------------------------------------------------------
counts_SE.norm.non_redundant %>%
	pivot_sample() %>%
	ggplot(aes(x=`Dim1`, y=`Dim2`, color=`Cell.type`)) +
  geom_point() +
  my_theme


## ----eval=FALSE---------------------------------------------------------------
#  counts = tidybulk_SAM_BAM(
#  	file_names,
#  	genome = "hg38",
#  	isPairedEnd = TRUE,
#  	requireBothEndsMapped = TRUE,
#  	checkFragLength = FALSE,
#  	useMetaFeatures = TRUE
#  )

## ----ensembl------------------------------------------------------------------
counts_ensembl %>% ensembl_to_symbol(ens)

## ----description--------------------------------------------------------------
counts_SE %>% describe_transcript() %>% select(transcript, description, everything())

## -----------------------------------------------------------------------------
sessionInfo()

