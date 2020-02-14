library("pasilla")
library(edgeR)
library(reshape)
### Reading data and sample annotation
pasCts = system.file("extdata",
                     "pasilla_gene_counts.tsv",
                     package = "pasilla",
                     mustWork = TRUE)
pasAnno = system.file(
  "extdata",
  "pasilla_sample_annotation.csv",
  package = "pasilla",
  mustWork = TRUE
)
cts = as.matrix(read.csv(pasCts, sep = "\t", row.names = "gene_id"))
dim(cts) # 14599     7
coldata = read.csv(pasAnno, row.names = 1)
coldata = coldata[, c("condition", "type")]
coldata$new.annot = row.names(coldata)
coldata$new.annot = gsub('fb', '', coldata$new.annot)
cts = cts[, match(coldata$new.annot, colnames(cts))]
# START WORKFLOW
plot_densities = function(){
### Remving lowly expressed genes
keep1 = filterByExpr(cts, group = factor(coldata$condition))
sum(keep1)
keep2 = apply(cts, 1, function(x)  length(x[x > 5]) > 2)
sum(keep2) # 7846
cts = cts[keep2,]
### Ploting distribution of samples
col.type = c('red', 'black')[coldata$type]
col.conditions = c('blue', 'cyan')[coldata$condition]
plot(density(log2(cts[, 1] + 1)), type = 'n', ylim = c(0, .25))
for (i in 1:ncol(cts))
  lines(density(log2(cts[, i] + 1)), col = col.type[i])
### TMM normalization
library(edgeR)
dge = DGEList(counts = cts,
              sample = coldata$condition,
              group = coldata$type)
dge = calcNormFactors(dge, method = "TMM")
logCPM = cpm(dge, log = TRUE, prior.count = 0.5)
plot(density(logCPM[, 1]), type = 'n', ylim = c(0, .25))
for (i in 1:ncol(cts))
  lines(density(logCPM[, i]), col = col.type[i])
}
plot_densities()
plot_MDS = function(){
  ### dimensionality reduction
library(GGally)
mds = plotMDS(logCPM, ndim = 3)
d = data.frame( 'cond' = coldata$condition,  'type' = coldata$type,  'data' = rep('CPM', 7),  'dim1' = mds$cmdscale.out[, 1],  'dim2' = mds$cmdscale.out[, 2],  'dim3' = mds$cmdscale.out[, 3])
ggpairs(d, columns = 4:ncol(d), ggplot2::aes(colour = type))
}
plot_MDS()
plot_adjusted_MDS = function(){
  ### ComBat
library(sva)
batch = coldata$type
mod.combat = model.matrix( ~ 1, data = coldata)
mod.condition = model.matrix( ~ condition, data = coldata)
combat.corrected = ComBat(  dat = logCPM,  batch = batch,  mod = mod.condition,  par.prior = TRUE,  prior.plots = FALSE)
mds.combat = plotMDS(combat.corrected, ndim = 3)
d2 = data.frame(
  'cond' = coldata$condition,
  'type' = coldata$type,
  'data' = rep('ComBat', 7),
  'dim1' = mds.combat$cmdscale.out[, 1],
  'dim2' = mds.combat$cmdscale.out[, 2],
  'dim3' = mds.combat$cmdscale.out[, 3]
)
final.d = rbind(d, d2)
library(tidyr)
final.d = gather(final.d, dim, dist, dim1:dim3, factor_key = TRUE)
final.d2 = gather(final.d, cond, type, cond:type, factor_key = TRUE)
final.d$new = paste0(final.d$cond, final.d$type)
ggplot(final.d2, aes(x = cond, y = dist, fill = type)) +
  geom_boxplot() +
  facet_wrap( ~ data + dim)
}
plot_adjusted_MDS()
test_abundance = function(){
    # DE (comparison 1)
design = model.matrix( ~ coldata$condition + coldata$type, data = coldata$condition)
dge = estimateGLMCommonDisp(dge, design)
dge = estimateGLMTagwiseDisp(dge, design)
fit = glmFit(dge, design)
lrt = glmLRT(fit, coef = 2)
de = topTags(lrt, n = nrow(dge$counts))
#hist(de.table$PValue)
list(
  de.table = de$table,
  de.genes = de.table[abs(de.table$logFC) >= 2,],
  de.genes.lable = de.table[abs(de.table$logFC) >= 3,]
)
  }
de_list = test_abundance()
de.table = de_list$de.table
de.genes = de_list$de.genes
de.genes.lable = de_list$de.genes.lable
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
ggplot(de.table, aes(x = logCPM, y = logFC, label = gene.lable)) +
  geom_point(aes(
    color = gene.de.color,
    size = size.point,
    alpha = size.point
  )) +
  ggrepel::geom_text_repel()
}
plot_MA()
plot_DE_comparative = function(){
  ### Boxplot of 6 DE genes
de.genes = row.names(de.genes.lable)[1:6]
### Raw
count.df = log2(dge$counts[de.genes ,] + 1)
colnames(count.df) = coldata$condition
count.df = melt(count.df)
count.df$data = 'count'
### cpm
cpm.df = logCPM[de.genes,]
colnames(cpm.df) = coldata$condition
cpm.df = melt(cpm.df)
cpm.df$data = 'cpm'
### combat
combat.df = combat.corrected[de.genes,]
colnames(combat.df) = coldata$condition
combat.df = melt(combat.df)
combat.df$data = 'combat'
### Boxplot of all data
final = rbind(count.df, cpm.df, combat.df)
final$data = factor(final$data, levels = c('count', 'cpm', 'combat'))
ggplot(final, aes(x = data, y = value, fill = Var2)) +
  geom_boxplot() +
  facet_wrap( ~ Var1)
}
plot_DE_comparative()
plot_heatmap = function(){
  ######## complex heatmap
  de.genes = row.names(de.genes)
  de.data = logCPM[de.genes ,]
  library(ComplexHeatmap)
  gene.labels = c(rep('AB', 28), rep('BA', 28))
  h1 = Heatmap(t(de.data), top_annotation = HeatmapAnnotation(labels = gene.labels))
  h2 = Heatmap(coldata$condition)
  h3 = Heatmap(coldata$type)
  draw(h1 + h2 + h3)
}
plot_heatmap()
