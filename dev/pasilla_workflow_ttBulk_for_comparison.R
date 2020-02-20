library(tidyverse)
library(ttBulk)
library("pasilla")
library(tictoc)
library(edgeR)
library(limma)

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
coldata = read.csv(pasAnno, row.names = 1)
coldata = coldata[, c("condition", "type")]
# Create ttBulk object
tt =
	cts %>%
	as_tibble(rownames = "transcript") %>%
	pivot_longer(names_to = "sample",
							 values_to = "count",
							 cols = -transcript) %>%
	left_join(coldata %>%
							as_tibble(rownames = "sample") %>%
							mutate(sample = gsub("fb", "", sample))) %>%
	# # attach symbol
	# left_join(tibble(
	# 	transcript = (.) %>% distinct(transcript) %>% pull(1),
	# 	symbol = flybaseR::id.converter((.) %>% distinct(transcript) %>% pull(1),
	# 																	symbols = TRUE,
	# 																	convert.into = "genes"
	# 	)
	# )) %>%
	# ttBulk
	mutate_if(is.character, as.factor) %>%
	ttBulk(sample, transcript, count)

time_df = tibble(step = "", time = list(), lines = NA, assignments = NA)

# START WORKFLOW
plot_densities = function(){
	# Normalise
tt_scaled = tt %>% scale_abundance()
# Plot densities
 p = tt_scaled %>%
	pivot_longer(	values_to = "count",	names_to = "Normalisation",	cols = c(count, `count scaled`)	) %>%
	ggplot(aes(count + 1, group = sample, color = type)) +
	facet_grid( ~ Normalisation) +
	geom_density() +
	scale_x_log10()

	tt_scaled
}
plot_MDS = function(){
	tt_mds = tt_scaled %>% 	reduce_dimensions(method = "MDS", .dims = 3)
# Visualise MDS for cell types
p = tt_mds %>%
	select(contains("Dim"), sample, type,  condition) %>%
	distinct() %>%
	GGally::ggpairs(	columns = 1:3,	ggplot2::aes(colour = type),	upper = list(continuous = GGally::wrap("cor", size = 3))	)

	tt_mds
}
plot_adjusted_MDS = function(){
# Adjust abundance for type, and check MDS densities
 p = tt_mds %>%
	adjust_abundance( ~ condition + type) %>%
	filter(`count scaled adjusted` %>% is.na %>% `!`) %>%
	reduce_dimensions(.abundance = `count scaled adjusted`,	method = "MDS",	.dims = 3) %>%

	# Plot
	select(contains("Dim"), sample, type,  condition) %>%
	distinct() %>%
	pivot_longer(names_to = "Dim", values_to = ".value", cols = contains("Dim")) %>%
	separate(Dim, c("Dim", "Adj"), sep = "\\.") %>%
	pivot_longer(names_to = "variation",	values_to = "Type of variation",	cols = c(type, condition)	) %>%
	ggplot(aes(y = .value, x = variation, fill = `Type of variation`)) +
	geom_boxplot() +
	facet_grid(Adj ~ Dim)

}
test_abundance = function(){
	# DE analyses
 tt %>% test_differential_abundance( ~ condition + type)
}
plot_MA = function(){
	# MA plot
p = tt_test %>%
	filter(!`filter out low counts`) %>%
	mutate(de = FDR < 0.05 & abs(logFC) >= 2) %>%
	distinct(transcript, logCPM, logFC, de) %>%
	mutate(transcript = ifelse(de &	abs(logFC) > 3, as.character(transcript), NA)) %>%
	ggplot(aes(x = logCPM, y = logFC, label = transcript)) +
	geom_point(aes(	color = de,	size = de,	alpha = de)) +
	ggrepel::geom_text_repel()
}
plot_DE_comparative = function(){
	# Plot top genes for raw, scaled and ajusted counts
 p = tt_test %>%
	inner_join((.) %>% distinct(transcript, PValue) %>% arrange(PValue) %>% head(6)) %>%
	pivot_longer(names_to = "Stage",	 values_to = "count",	 cols = starts_with("count")) %>%
	ggplot(aes(x = Stage, y = count + 1, fill = condition)) +
	#ggpol::geom_boxjitter(jitter.width = 0.1,  jitter.size = 0.2) +
	facet_wrap( ~ transcript) +
	scale_y_log10()
}
plot_heatmap = function(){
	# Heatmap
 p = tt_test %>%
	filter(FDR < 0.05 & abs(logFC) > 2) %>%
	tidyHeatmap::heatmap(	.horizontal = sample,		.vertical = transcript,		.abundance = `count`,annotation = c(condition, type),	log_transform = TRUE	)
}


tic()
tt_scaled = plot_densities()
time_df = time_df %>% bind_rows(tibble(step = "Normalisation", time = list(toc()), lines = 7, assignments = 1))

tic()
tt_mds = plot_MDS()
time_df = time_df %>% bind_rows(tibble(step = "Reduce dimensionality", time = list(toc()), lines = 5, assignments = 1))

tic()
plot_adjusted_MDS()
time_df = time_df %>% bind_rows(tibble(step = "Removal unwanted variation", time = list(toc()), lines = 12, assignments = 0))

tic()
tt_test = test_abundance()
time_df = time_df %>% bind_rows(tibble(step = "Test differential abundance", time = list(toc()), lines = 1, assignments = 0))

tic()
plot_MA()
time_df = time_df %>% bind_rows(tibble(step = "Plot MA", time = list(toc()), lines = 8, assignments = 0))

tic()
plot_DE_comparative()
time_df = time_df %>% bind_rows(tibble(step = "Plot results across stages", time = list(toc()), lines = 7, assignments = 0))

tic()
plot_heatmap()
time_df = time_df %>% bind_rows(tibble(step = "Plot heatmap", time = list(toc()), lines = 3, assignments = 0))



time_df  %>% mutate(step = factor(step, levels = unique(step))) %>%  saveRDS("dev/stats_pasilla_ttBulk.rds")

