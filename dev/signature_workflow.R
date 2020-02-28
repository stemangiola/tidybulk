library(tidyverse)
library(tidybulk)
library(furrr)
plan(multiprocess)
options(future.globals.maxSize = 50000 * 1024 ^ 2)

# Iterative version of Siberg function because fails
siberg_iterative = function(x) {
	if (x %>% unique %>% length %>% `<` (5))
		return(c(NA, NA))



	mu = NA
	max_i = ceiling(length(x) / 10)
	i = 0
	while (mu %>% is.na | i <= max_i) {
		res = SIBERG::SIBER(x, model = 'NB')

		BI = res[7]
		mu = res[1]
		x = x[-1]
		i = i + 1

	}


	if (mu %>% is.na & x %>% length %>% `<` (8))
		return(c(NA, NA))

	return(c(max(res[1], res[2]) / (min(res[1], res[2]) + 1),
					 res[7]))
}

my_theme =
	theme_bw() +
	theme(
		panel.border = element_blank(),
		axis.line = element_line(),
		panel.grid.major = element_line(size = 0.2),
		panel.grid.minor = element_line(size = 0.1),
		text = element_text(size=12),
		legend.position="bottom",
		#aspect.ratio=1,
		axis.text.x = element_text(angle = 30, hjust = 1),
		strip.background = element_blank(),
		axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
		axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
	)

data_hierarchy =
	ARMET::tree %>%
	data.tree::Clone() %>%
	ARMET::ToDataFrameTypeColFull(TRUE, "name") %>%
	pivot_longer(
		cols = -name,
		names_to = "level",
		values_to = "Cell type category",
		names_prefix = "level_"
	) %>%
	drop_na() %>%
	mutate(level = as.integer(level) - 1) %>%
	filter(level > 0) %>%
	dplyr::rename(`Cell type formatted` = name)

# Gather data
path = "/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/RNAseq-noise-model/big_data/tibble_cellType_files"

counts =
	dir(
		path = path,
		pattern="cellTypes_",
		full.names = TRUE
	) %>%
	future_map_dfr(
		~ read_csv(.x) %>%
			select(sample, `Cell type`, `Cell type formatted`, count, symbol, `Data base`)
	) %>%
	filter((!is.na(symbol)) & (symbol != "")) %>%
	filter(`Cell type formatted` %>% is.na %>% `!`) %>%

	# Select level 3
	inner_join(data_hierarchy) %>%
	filter(level==3) %>%
	filter(!`Cell type category` %in% c("t_cell", "b_cell"))

counts_proc =
	counts %>%

	# Create object
	tidybulk(sample, symbol, count) %>%

	# Sum redundant genes/isoforms
	aggregate_duplicates(aggregation_function = sum) %>%

	# Normalise
	scale_abundance()


# Remove redundancy
counts_non_red =
	counts_proc %>%
	remove_redundancy(
		method="correlation",
		correlation_threshold = 0.99,
		top=1000
	)

counts_non_red %>% saveRDS("dev/counts_non_red.rds", compress = "gzip")

counts_non_red = readRDS("dev/counts_non_red.rds")

counts_ct =
	counts_non_red %>%

	# Eliminate genes that are not in all cell types
	inner_join((.) %>%
						 	distinct(symbol, `Cell type formatted`) %>%
						 	count(symbol) %>%
						 	filter(n == max(n)))

# Calculate bimodality
bimodality =

	counts_ct %>%
	nest(data = -c(`Cell type formatted`, symbol)) %>%
	mutate(	bimodality_NB =
		map(
			data,
			~ {
				#print(.x)
				.x %>% pull(`count_scaled`) %>% as.integer %>%
				siberg_iterative() %>%
				`[` (1:2) %>%
				setNames(c("bimodality_NB_diff", "bimodality_NB")) %>%
				enframe() %>% spread(name, value)
		})
	) %>%
	select(-data) %>%
	unnest(bimodality_NB)

bimodality %>% saveRDS("dev/bimodality.rds")

bimodality = readRDS("dev/bimodality.rds")

# Assign values and filter
counts_ct_bm =
	counts_ct %>%
	left_join(bimodality) %>%

	# Too bimodal
	filter((bimodality_NB > 0.8 & bimodality_NB_diff > 20) | bimodality_NB_diff > 100)

#
# Plots and Study
counts_ct_bm %>%

	reduce_dimensions(sample, symbol, `count_scaled`, method = "tSNE") %>%
	select(contains("tSNE"), `Data base`, `Cell type formatted`) %>%
	distinct %>%
	ggplot(aes(x = `tSNE1`, y = `tSNE2`, color = `Cell type formatted`)) +
	geom_point(size =2) +
	theme_bw() +
	theme(
		panel.border = element_blank(),
		axis.line = element_line(),
		panel.grid.major = element_line(size = 0.2),
		panel.grid.minor = element_line(size = 0.1),
		text = element_text(size = 12),
		legend.position = "bottom",
		aspect.ratio = 1,
		#axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
		strip.background = element_blank(),
		axis.title.x  = element_text(margin = margin(
			t = 10,
			r = 10,
			b = 10,
			l = 10
		)),
		axis.title.y  = element_text(margin = margin(
			t = 10,
			r = 10,
			b = 10,
			l = 10
		))
	)

# Identify markers

# Setup Cell type category names
counts_light = counts_ct_bm %>% distinct(level, `Cell type category`, sample, symbol, count)
counts_light %>% saveRDS("dev/counts_light.rds")

counts_light = readRDS("dev/counts_light.rds")


markers =
	counts_light %>%
	filter(!`Cell type category` %in% c("t_cell", "b_cell")) %>%
	distinct(`Cell type category`) %>%
	pull(1) %>%
	gtools::permutations(n = length(.), r = 2, v = .) %>%
	as_tibble() %>%
	setNames(c("ct1", "ct2")) %>%
	mutate(de = map2(ct1, ct2,
									 ~ {
									 	counts_light %>%
									 		filter(`Cell type category` %in% c(.x, .y)) %>%
									 		rename(ct = `Cell type category`) %>%
									 		test_differential_abundance(
									 			~ 0 + ct,
									 			.contrasts = sprintf("ct%s - ct%s", .x, .y),
									 			fill_missing_values = TRUE,
									 			action = "only"
									 		) %>%
									 		filter(logFC > 0 & logFC > 2) %>%
									 		arrange(PValue) %>%
									 		slice(1:10) %>%
									 		select(symbol)
									 }))

markers %>% saveRDS("dev/markers_for_signature_workflow")

markers %>%
	mutate(comparison = 1:n()) %>%
	mutate(de = map(de, ~.x %>% slice(1))) %>%
	unnest(de) %>%
	filter(ct1=="monocyte") %>%
	unite(pair, c("ct1", "ct2"), remove = FALSE, sep = "\n") %>%
	gather(which, `Cell type category`, ct1, ct2) %>%
	left_join(counts_light, by = c("symbol", "Cell type category")) %>%
	ggplot(aes(y = count + 1, x = `Cell type category`)) +
	geom_boxplot() +
	facet_wrap(~pair + symbol, scales ="free_x", ncol = 5) +
	scale_y_log10() +
	my_theme

