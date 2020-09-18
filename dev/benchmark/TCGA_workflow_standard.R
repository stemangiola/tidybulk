library(reshape)
library(tidyverse)
library(tictoc)
library(ComplexHeatmap)
library(edgeR)
library(GGally)
library(sva)
library(tidybulk)

### Reading data and sample annotation
tt = readRDS("dev/PRAD_TCGA_counts.rds") %>% tidybulk(patient, transcript, count) %>% aggregate_duplicates() %>% filter(transcript %>% is.na %>% `!`)

cts = tt %>% distinct(patient, transcript, count) %>% spread(patient, count) %>% nanny::as_matrix(rownames = "transcript")

coldata = tt %>% pivot_sample() %>%
  
  # Add fake batch
  nest(data = -patient) %>%
  mutate(batch = sample(0:1, size = n(), replace = T)) %>%
  unnest(data) 

time_df = tibble(step = "", time = list(), lines = NA, assignments = NA)

# START WORKFLOW
plot_densities = function(){

### Remving lowly expressed genes
keep1 = filterByExpr(cts)
sum(keep1)
keep2 = apply(cts, 1, function(x)  length(x[x > 5]) > 2)
sum(keep2) # 7846
cts = cts[keep2,]

### Ploting distribution of samples
col.type = c('red', 'black')[coldata$PFI.2]
plot(density(log2(cts[, 1] + 1)), type = 'n', ylim = c(0, .25))
for (i in 1:ncol(cts))
  lines(density(log2(cts[, i] + 1)), col = col.type[i])

### TMM normalization

dge = DGEList(counts = cts,
              sample = coldata,
              group = coldata$PFI.2)
dge = calcNormFactors(dge, method = "TMM")
logCPM = cpm(dge, log = TRUE, prior.count = 0.5)
plot(density(logCPM[, 1]), type = 'n', ylim = c(0, .25))
for (i in 1:ncol(cts))
  lines(density(logCPM[, i]), col = col.type[i])

list(dge = dge, logCPM = logCPM)

}
plot_MDS = function(){

### dimensionality reduction

mds = plotMDS(logCPM, ndim = 3)
d = data.frame(   'type' = coldata$PFI.2,    'dim1' = mds$cmdscale.out[, 1],  'dim2' = mds$cmdscale.out[, 2],  'dim3' = mds$cmdscale.out[, 3])
p = ggpairs(d, columns = 2:ncol(d), ggplot2::aes(colour = as.factor(PFI.2)))

d
}
plot_adjusted_MDS = function(){
  ### ComBat

    
  batch = coldata$batch
  mod.combat = model.matrix( ~ 1, data = coldata)
  mod.condition = model.matrix( ~ PFI.2, data = coldata)
  combat.corrected = ComBat(  dat = logCPM,  batch = batch,  mod = mod.condition,  par.prior = TRUE,  prior.plots = FALSE)
  mds.combat = plotMDS(combat.corrected, ndim = 3)
  d2 = data.frame(  'cond' = coldata$PFI.2,  'type' = coldata$batch,   'dim1' = mds.combat$cmdscale.out[, 1],  'dim2' = mds.combat$cmdscale.out[, 2], 'dim3' = mds.combat$cmdscale.out[, 3])
  final.d = bind_rows(d, d2)
  final.d = gather(final.d, dim, dist, dim1:dim3, factor_key = TRUE)
  final.d2 = gather(final.d, cond, batch, cond:batch, factor_key = TRUE)
  final.d$new = paste0(final.d$cond, final.d$batch)
  p = ggplot(final.d2, aes(x = cond, y = dist, fill = batch)) +
    geom_boxplot() +
    facet_wrap( ~  dim)
  
  combat.corrected
  
}
test_abundance = function(){
    # DE (comparison 1)
design = model.matrix( ~ coldata$PFI.2, data = coldata)
dge = estimateGLMCommonDisp(dge, design)
dge = estimateGLMTagwiseDisp(dge, design)
fit = glmQLFit(dge, design)
lrt = glmQLFTest(fit, coef = 2)
de = topTags(lrt, n = nrow(dge$counts))
#hist(de.table$PValue)

de.table  = de$table

list(
  de.table = de.table,
  de.genes = de.table[abs(de.table$logFC) >= 2,],
  de.genes.lable = de.table[abs(de.table$logFC) >= 3,]
)
}
plot_MA = function(){
  ### MA plot
  n.genes = nrow(dge$counts)
  gene.de = rep(NA, n.genes)
  gene.de[which(row.names(de.table) %in% row.names(de.genes))] = row.names(de.genes)
  gene.de.color = rep('black', n.genes)
  gene.de.color[which(row.names(de.table) %in% row.names(de.genes))] = 'red'
  size.point = ifelse(gene.de.color == 'black', .1, .2)
  gene.lable = rep(NA, n.genes)
  gene.lable[which(row.names(de.table) %in% row.names(de.genes.lable))] = row.names(de.genes.lable)
  p = ggplot(de.table, aes(x = logCPM, y = logFC, label = gene.lable)) +
    geom_point(aes(
      color = gene.de.color,
      size = size.point,
      alpha = size.point
    )) +
    ggrepel::geom_text_repel()
}
plot_DE_comparative = function(){
  ### Boxplot of 6 DE genes
  de.genes = row.names(de.genes.lable)[1:6]
  ### Raw
  count.df = log2(dge$counts[de.genes ,] + 1)
  colnames(count.df) = coldata$PFI.2
  count.df = melt(count.df)
  count.df$data = 'count'
  ### cpm
  cpm.df = logCPM[de.genes,]
  colnames(cpm.df) = coldata$PFI.2
  cpm.df = melt(cpm.df)
  cpm.df$data = 'cpm'
  ### combat
  combat.df = combat.corrected[de.genes,]
  colnames(combat.df) = coldata$PFI.2
  combat.df = melt(combat.df)
  combat.df$data = 'combat'
  ### Boxplot of all data
  final = rbind(count.df, cpm.df, combat.df)
  final$data = factor(final$data, levels = c('count', 'cpm', 'combat'))
  p = ggplot(final, aes(x = data, y = value, fill = X2)) +
    geom_boxplot() +
    facet_wrap( ~ X1)
  
  de.genes
}
plot_heatmap = function(){
  ######## complex heatmap
  de.data = logCPM[de.genes ,]
  
  gene.labels = c(rep('AB', floor(length(de.genes)/2)), rep('BA', ceiling(length(de.genes)/2)))
  h1 = Heatmap(t(de.data), top_annotation = HeatmapAnnotation(labels = gene.labels))
  h2 = Heatmap(coldata$PFI.2)
  h3 = Heatmap(coldata$type)
  p = draw(h1 + h2 + h3)
}

tic()
pd_res = plot_densities()
time_df = time_df %>% bind_rows(tibble(step = "Normalisation", time = list(toc()), lines = 19, assignments = 8))

logCPM = pd_res$logCPM
dge = pd_res$dge

tic()
d  = plot_MDS()
time_df = time_df %>% bind_rows(tibble(step = "Reduce dimensionality", time = list(toc()), lines = 3, assignments = 2))

tic()
combat.corrected= plot_adjusted_MDS()
time_df = time_df %>% bind_rows(tibble(step = "Removal unwanted variation", time = list(toc()), lines = 13, assignments = 10))

tic()
de_list = test_abundance()
time_df = time_df %>% bind_rows(tibble(step = "Test differential abundance", time = list(toc()), lines = 7, assignments = 6))

de.table = de_list$de.table
de.genes = de_list$de.genes
de.genes.lable = de_list$de.genes.lable

tic()
plot_MA()
time_df = time_df %>% bind_rows(tibble(step = "Plot MA", time = list(toc()), lines = 11, assignments = 8))

tic()
de.genes = plot_DE_comparative()
time_df = time_df %>% bind_rows(tibble(step = "Plot results across stages", time = list(toc()), lines = 18, assignments = 15))

tic()
plot_heatmap()
time_df = time_df %>% bind_rows(tibble(step = "Plot heatmap", time = list(toc()), lines = 6, assignments = 5))

saveRDS(time_df, "dev/stats_TCGA_standard.rds")

# 
# > time_df = time_df %>% bind_rows(tibble(step = "Normalisation", time = list(toc()), lines = 19, assignments = 8))
# 18.725 sec elapsed
# > 
#   > logCPM = pd_res$logCPM
# > dge = pd_res$dge
# > 
#   > tic()
# > d  = plot_MDS()
# > time_df = time_df %>% bind_rows(tibble(step = "Reduce dimensionality", time = list(toc()), lines = 3, assignments = 2))
# 212.277 sec elapsed
# > 
#   > tic()
# > de_list = test_abundance()
# > time_df = time_df %>% bind_rows(tibble(step = "Test differential abundance", time = list(toc()), lines = 7, assignments = 6))
# 119.962 sec elapsed