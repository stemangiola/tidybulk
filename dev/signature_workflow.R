library(tidyverse)
library(tidybulk)
library(furrr)
plan(multicore)
options(future.globals.maxSize = 50000 * 1024 ^ 2)
library(RColorBrewer)

# Iterative version of Siberg function because fails
siberg_iterative = function(x) {
	if (x %>% unique %>% length %>% `<` (5))
		return(c(NA, NA))



	mu = NA
	max_i = ceiling(length(x) / 10)
	#max_i = 10
	i = 0
	while (mu %>% is.na | i <= max_i) {
		res = SIBERG::SIBER(x, model = 'NB')

		BI = res[7]
		mu = res[1]
		x = x[-1]
		i = i + 1

	}


	if (mu %>% is.na & x %>% length %>% `<` (5))
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

# Expand palette
colourCount = data_hierarchy %>% filter(level ==3) %>% distinct(`Cell type category`) %>% nrow
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

# Gather data
# path = "/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/RNAseq-noise-model/big_data/tibble_cellType_files"
# counts =
# 	dir(
# 		path = path,
# 		pattern="cellTypes_",
# 		full.names = TRUE
# 	) %>%
# 	future_map_dfr(
# 		~ read_csv(.x) %>%
# 			select(sample, `Cell type`, `Cell type formatted`, count, symbol, `Data base`)
# 	) %>%
# 	filter((!is.na(symbol)) & (symbol != "")) %>%
# 	filter(`Cell type formatted` %>% is.na %>% `!`) %>%
#
# 	# Select level 3
# 	inner_join(data_hierarchy) %>%
# 	filter(level==3) %>%
# 	filter(!`Cell type category` %in% c("t_cell", "b_cell"))
#
# counts_proc =
# 	counts %>%
#
# 	# Create object
# 	tidybulk(sample, symbol, count) %>%
#
# 	# Sum redundant genes/isoforms
# 	aggregate_duplicates(aggregation_function = sum) %>%
#
# 	# Normalise
# 	scale_abundance()

# # Remove redundancy
# counts_non_red =
# 	counts_proc %>%
# 	remove_redundancy(
# 		method="correlation",
# 		correlation_threshold = 0.99,
# 		top=1000
# 	)
#
# counts_non_red %>% saveRDS("dev/counts_non_red.rds", compress = "gzip")

counts_non_red = readRDS("dev/counts_non_red.rds")

# Plots and Study
(counts_non_red %>%

		reduce_dimensions(sample, symbol, `count_scaled`, method = "tSNE") %>%
		pivot_sample() %>%
		ggplot(aes(x = `tSNE1`, y = `tSNE2`, color = `Cell type category`)) +
		geom_point(size =2) +
		scale_color_manual(values = getPalette(colourCount)) +
		my_theme + theme(aspect.ratio=1)) %>%
	ggsave(	"dev/signature_p1.pdf",	plot = .,	useDingbats=FALSE,	units = c("mm"),	width = 183/2 ,	limitsize = FALSE)


#
# counts_non_red_de =
# 	counts_non_red %>%
# 	filter(!`Cell type category` %in% c("t_cell", "b_cell")) %>%
# 	distinct(`Cell type category`) %>%
# 	pull(1) %>%
# 	gtools::permutations(n = length(.), r = 2, v = .) %>%
# 	as_tibble() %>%
# 	setNames(c("ct1", "ct2")) %>%
# 	mutate(contrast = sprintf("ct%s - ct%s", ct1, ct2)) %>%
# 	mutate(de = pmap(list(ct1, ct2, contrast),
# 									 ~ 	counts_non_red %>%
# 									 	filter(`Cell type category` %in% c(..1, ..2)) %>%
# 									 	rename(ct = `Cell type category`) %>%
#
# 									 	test_differential_abundance(
# 									 		~ 0 + ct,
# 									 		.contrasts = ..3,
# 									 		fill_missing_values = TRUE
# 									 	)
# 									 # %>%
# 									 # filter(!!(as.symbol(sprintf("logFC_%s", ..3))) > 2) %>%
# 									 # arrange(!!(as.symbol(sprintf("PValue_%s", ..3)))) %>%
# 									 # slice(1:10) %>%
# 									 # select(symbol)
# 	))
#
# counts_non_red_de %>% saveRDS("dev/counts_non_red_de.rds")

counts_non_red_de <- readRDS("dev/counts_non_red_de.rds")




# # Calculate bimodality
# bimodality =
#
# 	counts_non_red %>%
# 	#keep_variable(top = 5000) %>%
# 	tidybulk:::drop_class(c("tidybulk", "tt")) %>%
# 	tidybulk:::drop_internals() %>%
# 	nest(data = -c(`Cell type formatted`, symbol)) %>%
#
# 	#slice(1:10) %>%
# 	mutate(	bimodality_NB =
# 		map(
# 			data,
# 			~ tryCatch(
# 							.x %>% pull(`count_scaled`) %>% as.integer %>%
# 								siberg_iterative() %>%
# 								`[` (1:2) , error=function(e) c(NA, NA))		%>%
# 							setNames(c("bimodality_NB_diff", "bimodality_NB")) %>%
# 							enframe() %>% spread(name, value)
#
# 		)
# 	) %>%
# 	select(-data) %>%
# 	unnest(bimodality_NB)
#
# bimodality %>% saveRDS("dev/bimodality.rds")
#
bimodality = readRDS("dev/bimodality.rds")

non_bimodal =
	bimodality %>%
	add_count(symbol) %>%
	filter(n==max(n)) %>%
	mutate(bimodal = ((bimodality_NB > 0.8 & bimodality_NB_diff > 20) | bimodality_NB_diff > 100) ) %>%
	nest(data = -symbol) %>%
	mutate(how_many_bimod = map_int(data, ~ .x %>% pull(bimodal) %>% sum(na.rm=T))) %>%
	filter(how_many_bimod == 0)


associated_with_data_base =

	counts_non_red %>%
	tidybulk:::drop_class(c("tidybulk", "tt")) %>%
	tidybulk:::drop_internals() %>%


	nest(data = -c(`Cell type formatted`, symbol)) %>%
	filter(map_int(data, ~.x %>% distinct(`Data base`) %>% nrow) > 1) %>%

	#slice(1:10) %>%
	mutate(	anova =
						map(
							data,
							~ .x %>%
								mutate(count_scaled_log = log(count_scaled + 1)) %>%
								mutate(data_base = `Data base`) %>%
								aov(count_scaled_log ~ data_base, data = .) %>%
								broom::tidy() %>%
								filter(term=="data_base")
						)
	) %>%
	select(-data) %>%
	unnest(anova) %>%
	mutate(FDR = p.adjust(p.value, method="BH"))

associated_with_data_base %>% saveRDS("dev/associated_with_data_base.rds")

markers =
	counts_non_red_de %>%
	mutate(de =
				 	map(
				 		de,
				 		~ .x %>%
				 			filter(ct1 == `Cell type formatted`) %>%
				 			filter(symbol %in% (non_bimodal %>% pull(symbol))) %>%
				 			select(symbol, starts_with("logFC"), starts_with("PValue"), starts_with("FDR"),  `Cell type formatted`)  %>%
				 			distinct %>%
				 			setNames(c("symbol" ,   "logFC" , "PValue" ,"FDR" , "Cell type formatted"   )) %>%
				 			filter(FDR < 0.05 & logFC > 0) %>%
				 			arrange(FDR) %>%
				 			mutate(i = 1:n()) )) %>%
	unnest(de) %>%
	filter(symbol %in% (non_bimodal %>% pull(symbol))) %>%
	anti_join(associated_with_data_base %>% filter(FDR<0.00001) %>% distinct(symbol)) %>%
	filter(!`Cell type formatted` %in% c("t_cell", "b_cell"))


markers %>% saveRDS("dev/markers_for_signature_workflow.rds")


(markers %>%
		filter(ct1=="monocyte") %>%
		group_by(ct2) %>%
		arrange(i) %>%
		slice(1) %>%
		ungroup() %>%
		unite(pair, c("ct1", "ct2"), remove = FALSE, sep = "\n") %>%
		gather(which, `Cell type category`, ct1, ct2) %>%
		distinct(pair ,   contrast    ,  symbol ,   which, `Cell type category`) %>%
		left_join(counts_non_red, by = c("symbol", "Cell type category")) %>%

		# Sort labels
		#mutate(`Cell type category` = factor(`Cell type category`, levels=c("monocyte", counts_non_red %>% filter(`Cell type formatted` != "monocyte") %>% pull(`Cell type formatted`) %>% unique %>% sort))) %>%
		ggplot(aes(y = count_scaled + 1, x = `Cell type category`, fill = `Cell type category`)) +
		geom_boxplot() +
		facet_wrap(~pair+ symbol, scales ="free_x", nrow = 2) +
		scale_y_log10() +
		scale_fill_manual(values = getPalette(colourCount)) +
		my_theme +
		theme(	axis.text.x = element_text(angle = 30, hjust = 1))) %>%
	ggsave(	"dev/signature_p2.pdf",	plot = .,	useDingbats=FALSE,	units = c("mm"),	width = 183 ,	limitsize = FALSE)

# Plots and Study
(markers %>%
		unite(pair, c("ct1", "ct2"), remove = FALSE, sep = "\n") %>%
		gather(which, `Cell type category`, ct1, ct2) %>%


		filter(i < 6) %>%

		distinct(symbol) %>%

		left_join(counts_non_red, by = c("symbol"))  %>%

		# Impute missng values
	mutate(ct = `Cell type category`) %>%
		tidybulk:::fill_NA_using_formula(~ct, sample, symbol, count_scaled) %>%

		reduce_dimensions(sample, symbol, count_scaled, method = "tSNE") %>%
		pivot_sample(sample) %>%
		ggplot(aes(x = `tSNE1`, y = `tSNE2`, color = `Cell type category`)) +
		geom_point(size =2) +
		scale_color_manual(values = getPalette(colourCount)) +
		my_theme + theme(aspect.ratio=1)) %>%
	ggsave(	"dev/signature_p3.pdf",	plot = .,	useDingbats=FALSE,	units = c("mm"),	width = 183/2 ,	limitsize = FALSE)


