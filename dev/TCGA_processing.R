library(TCGAbiolinks)
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(tidybulk)
library(tidySummarizedExperiment)


query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  platform = "Illumina HiSeq",
                  file.type  = "normalized_results",
                  experimental.strategy = "RNA-Seq",
                  legacy = TRUE)


GDCdownload(query, method = "api", files.per.chunk = 10)
counts_se <- GDCprepare(query)

query_clin <- GDCquery(project = "TCGA-BRCA",
                       data.category = "Clinical",
                       data.type = "Clinical data",
                       file.type = "txt",
                       legacy = TRUE)


GDCdownload(query_clin)

clinical_list <- GDCprepare(query_clin)

clinical_patient_brca <- clinical_list$clinical_patient_brca

groups <- c("HER2_pos", "HER2_low", "HER2_neg")

tcga_her2 <- clinical_patient_brca %>%
  mutate(her2_group = case_when((her2_ihc_score == 3 | her2_fish_status == "Positive") ~ "HER2_pos",
                                ((her2_ihc_score == 1 | her2_ihc_score == 2) & her2_fish_status == "Negative") ~ "HER2_low",
                                her2_ihc_score == 0 ~ "HER2_neg"
  )) %>%
  filter(her2_group %in% groups) %>%
  select(bcr_patient_barcode, her2_group)

tcga_her2

tcga_joined <- inner_join(counts_se, tcga_her2, by = c("barcode" = "bcr_patient_barcode"))

counts_se %>% select(barcode)

df = readRDS("dev/TCGA_breast.rds") %>% as_tibble()

df %>% aggregate_duplicates(.sample, .feature, normalized_count)
