library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(tidybulk)
library(tidySummarizedExperiment)

tcga_joined  <- readRDS("dev/counts_tt.rds")


counts_dupsrem <- counts_tt %>% aggregate_duplicates()
