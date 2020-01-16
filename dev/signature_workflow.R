library(tidyverse)
library(ttBulk)
library(furrr)
plan(multiprocess)
options(future.globals.maxSize = 50000 * 1024 ^ 2)

# Iterative version of Siberg function because fails
siberg_iterative = function(x) {

	if(x %>% unique %>% length %>% `<` (5)) return(c(NA, NA))



	mu = NA
	max_i = ceiling(length(x)/10)
	i = 0
	while (mu %>% is.na | i <= max_i) {
		res = SIBERG::SIBER(x, model='NB')

		BI = res[7]
		mu = res[1]
		x = x[-1]
		i = i+1

	}


	if(mu %>% is.na & x %>% length %>% `<` (8)) return(c(NA, NA))

	return(c(
		max(res[1], res[2]) / (min(res[1], res[2]) + 1),
		res[7]
	))
}

data_hierarchy =
	ARMET::tree %>%
	data.tree::Clone() %>%
	ARMET::ToDataFrameTypeColFull(F, "name") %>%
	pivot_longer(cols = -name, names_to = "level", values_to = "Cell type category", names_prefix = "level_") %>%
	drop_na() %>%
	mutate(level = as.integer(level) -1) %>%
	filter(level > 0) %>%
	dplyr::rename(`Cell type formatted` = name)

# Gather data
path = "/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/RNAseq-noise-model/big_data/tibble_cellType_files"

counts =
	dir(
		path = path,
		pattern="cellTypes_",
		full.names = T
	) %>%
	future_map_dfr(
		~ read_csv(.x) %>%
			select(sample, `Cell type`, `Cell type formatted`, count, symbol, `Data base`)
	) %>%
	filter((!is.na(symbol)) & (symbol != "")) %>%
	filter(`Cell type formatted` %>% is.na %>% `!`)

counts_proc =
	counts %>%

	# Create object
	ttBulk(sample, symbol, count) %>%

	# Sum redundant genes/isoforms
	aggregate_duplicates(aggregation_function = sum) %>%

	# Normalise
	scale_abundance() %>%

	# Remove redundant samples
	inner_join(
			filter_variable((.), top = 1000) %>%
			remove_redundancy(method="correlation", correlation_threshold = 0.95) %>%
				distinct(sample)
	)

counts_proc %>% saveRDS("dev/counts_proc.rds", compress = "gzip")

counts_ct =
	counts_proc %>%

	# Setup Cell type category names
	inner_join(data_hierarchy) %>%

	# Eliminate genes that are not in all cell types
	inner_join(
		(.) %>%
			distinct(symbol, `Cell type category`) %>%
			count(symbol) %>%
			filter(n == max(n))
	)

# Calculate bimodality
bimodality =

	counts_ct %>%
	nest(data = -c(`Cell type category`, symbol, level)) %>%
	mutate(	bimodality_NB =
		future_map(
			data,
			~ .x %>% pull(`count normalised`) %>% as.integer %>%
				siberg_iterative() %>%
				`[` (1:2) %>%
				setNames(c("bimodality_NB_diff", "bimodality_NB")) %>%
				enframe() %>% spread(name, value)
		)
	) %>%
	select(-data) %>%
	unnest(bimodality_NB)

# Assign values
counts_ct_bm = counts_ct %>% left_join(bimodality)

# Plots and Study
counts_ct_bm %>%

    reduce_dimensions(method = "tSNE") %>%
    select(contains("tSNE"), `Data base`, `Cell type formatted`) %>%
    distinct %>%
    ggplot(aes(x = `tSNE1`, y = `tSNE2`, color=`Cell type formatted`)) + geom_point(size=2) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(),
      panel.grid.major = element_line(size = 0.2),
      panel.grid.minor = element_line(size = 0.1),
      text = element_text(size=12),
      legend.position="bottom",
      aspect.ratio=1,
      #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      strip.background = element_blank(),
      axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
      axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
    )

