library(tidyverse)
library(tidyBulk)

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

library("pasilla")
pasCts <- system.file("extdata",
											"pasilla_gene_counts.tsv",
											package="pasilla", mustWork=TRUE)
pasAnno <- system.file("extdata",
											 "pasilla_sample_annotation.csv",
											 package="pasilla", mustWork=TRUE)
cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")]

# Create tidyBulk object
tt =
	cts %>%
	as_tibble(rownames = "transcript") %>%
	pivot_longer(names_to = "sample", values_to = "count", cols=-transcript) %>%
	left_join(
		coldata %>%
			as_tibble(rownames = "sample") %>%
			mutate(sample = gsub("fb", "", sample))
	) %>%

	# attach symbol
	left_join(
		tibble(
			transcript = (.) %>% distinct(transcript) %>% pull(1),
			symbol = flybaseR::id.converter((.) %>% distinct(transcript) %>% pull(1), symbols = TRUE, convert.into = "genes")
		)
	) %>%

	# tidyBulk
	mutate_if(is.character, as.factor) %>%
	tidyBulk(sample, transcript, count)

# Normalise
tt_scaled = tt %>% scale_abundance()

# Plot densities
tt_scaled %>%
	pivot_longer(values_to = "count", names_to = "Normalisation", cols = c(count, `count_scaled`)) %>%
	ggplot(aes(count + 1, group=sample, color=type)) +
	facet_grid(~Normalisation) +
	geom_density() +
	scale_x_log10() +
	my_theme

tt_mds = tt_scaled %>% 	reduce_dimensions(method="MDS", .dims = 3)

# Visualise MDS for cell types
tt_mds %>%
	select(contains("Dim"), sample, type,  condition ) %>%
	distinct() %>%
	GGally::ggpairs(columns = 1:3, ggplot2::aes(colour=condition), upper = list(continuous = GGally::wrap("cor", size = 3))) +
	theme_bw() +
	theme(
		panel.border = element_blank(),
		axis.line = element_line(),
		panel.grid.major = element_line(size = 0.2),
		panel.grid.minor = element_line(size = 0.1),
		legend.position="bottom",
		strip.background = element_blank()
	)

# Visualise MDS for cell types
tt_mds %>%
	select(contains("Dim"), sample, type,  condition ) %>%
	distinct() %>%
	GGally::ggpairs(columns = 1:3, ggplot2::aes(colour=type), upper = list(continuous = GGally::wrap("cor", size = 3))) +
	theme_bw() +
	theme(
		panel.border = element_blank(),
		axis.line = element_line(),
		panel.grid.major = element_line(size = 0.2),
		panel.grid.minor = element_line(size = 0.1),
		legend.position="bottom",
		strip.background = element_blank()
	)

# Adjust abundance for type, and check MDS densities
tt_adj = tt_mds %>%	adjust_abundance(~ condition + type)

tt_adj %>%
	filter( `count_scaled_adjusted` %>% is.na %>% `!`) %>%
	reduce_dimensions(.abundance = `count_scaled_adjusted`, method="MDS", .dims = 3) %>%

	# Plot
	select(contains("Dim"), sample, type,  condition ) %>%
	distinct() %>%
	pivot_longer(names_to = "Dim", values_to = ".value", cols = contains("Dim")) %>%
	separate(Dim, c("Dim", "Adj"), sep="\\.") %>%
	mutate(Adj = ifelse(Adj == "y", "Non-adjusted", "adjusted") %>% factor(labels = c("Non-adjusted", "adjusted"))) %>%
	pivot_longer(names_to = "variation", values_to = "Type of variation", cols = c(type, condition)) %>%
	ggplot(aes(y = .value, x = variation, fill = `Type of variation`)) +
	geom_boxplot() +
	#ggpol::geom_boxjitter(jitter.width = 0.1,  jitter.size = 0.5, jitter.stroke = 0	) +
	#geom_jitter() +
	facet_grid(Adj ~ Dim ) +
	my_theme

# DE analyses
tt_test = tt_adj %>% test_differential_abundance(~ condition + type)

# MA plot
tt_test %>%
	filter(!`lowly_abundant.x`) %>%
	mutate(de = FDR<0.05 & abs(logFC) >=2) %>%
	distinct(transcript, symbol, logCPM, logFC, de) %>%
	mutate(symbol = ifelse(de & abs(logFC)>3, as.character(symbol), NA)) %>%
	ggplot(aes(x = logCPM, y = logFC, label=symbol)) +
	geom_point(aes(color = de, size = de, alpha=de)) +
	ggrepel::geom_text_repel() +
	scale_color_manual(values=c("black", "#e11f28")) +
	scale_size_discrete(range = c(0, 2)) +
	my_theme

# Plot top genes for raw, scaled and ajusted counts
tt_test %>%
	inner_join( (.) %>% distinct(symbol, PValue) %>% arrange(PValue) %>% head(6)) %>%
	pivot_longer(names_to = "Stage", values_to = "count", cols = starts_with("count")) %>%
	ggplot(aes(x = Stage, y = count + 1, fill = condition)) +
	ggpol::geom_boxjitter(jitter.width = 0.1,  jitter.size = 0.2	) +
	facet_wrap(~symbol) +
	scale_y_log10() +
	my_theme +
	theme(axis.text.x =element_text(angle = 60, hjust = 1, vjust = 1))

# Heatmap
tt_test %>%
	filter(FDR < 0.05 & abs(logFC) > 2) %>%
	tidyHeatmap::plot_heatmap(		.horizontal = sample,		.vertical = symbol,		.abundance = `count_scaled_adjusted`,
		annotation = c(condition, type),
		log_transform = TRUE
	)

