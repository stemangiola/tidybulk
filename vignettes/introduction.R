## ----echo=FALSE, out.width = "800px"------------------------------------------
knitr::include_graphics("../man/figures/new_SE_usage-01.png")

## ----echo=FALSE, include=FALSE------------------------------------------------
library(knitr)
# knitr::opts_chunk$set(cache = TRUE, warning = FALSE,
#                       message = FALSE, cache.lazy = FALSE)

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

data(se_mini)
tibble_counts = tidybulk::se_mini |> tidybulk() |> as_tibble()


## ----eval=FALSE---------------------------------------------------------------
#  BiocManager::install("tidybulk")

## ----eval=FALSE---------------------------------------------------------------
#  devtools::install_github("stemangiola/tidybulk")

## -----------------------------------------------------------------------------
se_mini

## -----------------------------------------------------------------------------
class(se_mini)

## ----eval=FALSE---------------------------------------------------------------
#  se_mini |>	get_bibliography()

## ----aggregate, message=FALSE, warning=FALSE, results='hide', class.source='yellow'----
rowData(se_mini)$gene_name = rownames(se_mini)
se_mini.aggr = se_mini |> aggregate_duplicates(.transcript = gene_name)

## ----aggregate long, eval=FALSE-----------------------------------------------
#  temp = data.frame(
#  	symbol = dge_list$genes$symbol,
#  	dge_list$counts
#  )
#  dge_list.nr <- by(temp,	temp$symbol,
#  	function(df)
#  		if(length(df[1,1])>0)
#  			matrixStats:::colSums(as.matrix(df[,-1]))
#  )
#  dge_list.nr <- do.call("rbind", dge_list.nr)
#  colnames(dge_list.nr) <- colnames(dge_list)

## ----normalise----------------------------------------------------------------
se_mini.norm = se_mini.aggr |> identify_abundant(factor_of_interest = condition) |> scale_abundance()

## ----normalise long, eval=FALSE-----------------------------------------------
#  library(edgeR)
#  
#  dgList <- DGEList(count_m=x,group=group)
#  keep <- filterByExpr(dgList)
#  dgList <- dgList[keep,,keep.lib.sizes=FALSE]
#  [...]
#  dgList <- calcNormFactors(dgList, method="TMM")
#  norm_counts.table <- cpm(dgList)

## ----include=FALSE------------------------------------------------------------
se_mini.norm |> select(`count`, count_scaled, .abundant, everything())

