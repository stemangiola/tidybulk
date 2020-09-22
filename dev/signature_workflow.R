library(tidyverse)
library(tidybulk)
library(furrr)
plan(multicore, workers=15)
options(future.globals.maxSize = 50000 * 1024 ^ 2)
library(RColorBrewer)


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

ToDataFrameTypeColFull = function(tree, fill = TRUE, ...) {
	t = tree %>% data.tree::Clone()
	
	1:(t %$% Get("level") %>% max) %>%
		map_dfr(
			~ data.tree::Clone(t) %>%
				{
					data.tree::Prune(., function(x)
						x$level <= .x + 1)
					.
				} %>%
				data.tree::ToDataFrameTypeCol("name") %>%
				as_tibble
			
		) %>%
		distinct() %>%
		
		when(
			fill & ("level_3" %in% colnames(.)) ~ mutate(., level_3 = ifelse(level_3 %>% is.na, level_2, level_3)),
			fill & ("level_4" %in% colnames(.)) ~ mutate(., level_4 = ifelse(level_4 %>% is.na, level_3, level_4)),
			fill & ("level_5" %in% colnames(.)) ~ mutate(., level_5 = ifelse(level_5 %>% is.na, level_4, level_5)),
			fill & ("level_6" %in% colnames(.)) ~ mutate(., level_6 = ifelse(level_6 %>% is.na, level_5, level_6)),
			TRUE ~ (.)
		) %>%
		dplyr::select(..., everything())
	
}

data_hierarchy =
	ARMET::tree %>%
	data.tree::Clone() %>%
	ToDataFrameTypeColFull(TRUE) %>%
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
	)  %>%
	rename(cell_type = `Cell type category`, data_base = `Data base`) %>%
	mutate_if(is.character, as.factor) %>%
	tidybulk:::drop_class(c("tidybulk", "tt")) %>%
	tidybulk(sample, symbol, count, count_scaled)
#
# counts_non_red %>% saveRDS("dev/counts_non_red.rds", compress = "gzip")

counts_non_red = readRDS("dev/counts_non_red.rds")

# Plots and Study
(counts_non_red %>%
		reduce_dimensions(method = "tSNE", action="get") %>%
		ggplot(aes(x = `tSNE1`, y = `tSNE2`, color = cell_type)) +
		geom_point(size =2) +
		scale_color_manual(values = getPalette(colourCount)) +
		my_theme + theme(aspect.ratio=1)) %>%
	ggsave(	"dev/signature_p1.pdf",	plot = .,	useDingbats=FALSE,	units = c("mm"),	width = 183/2 ,	limitsize = FALSE)



# counts_non_red_de %>% saveRDS("dev/counts_non_red_de.rds")

counts_non_red_de <- readRDS("dev/counts_non_red_de.rds")




associated_with_data_base =

	counts_non_red %>%
	nest(data = -c(cell_type, symbol)) %>%

	# Eliminate one database only
	filter(map_int(data, ~.x %>% distinct(data_base) %>% nrow) > 1) %>%

	mutate(	anova =
						map(
							data,
							~ .x %>%
								aov(log(count_scaled + 1) ~ data_base, data = .) %>%
								broom::tidy() %>%
								filter(term=="data_base")
						)
	) %>%
	select(-data) %>%
	unnest(anova) %>%
	mutate(FDR = p.adjust(p.value, method="BH"))


associated_with_data_base %>% saveRDS("dev/associated_with_data_base.rds")


keep_unimodal = associated_with_data_base %>% filter(FDR>0.00001) %>% distinct(symbol) %>% pull(1)

counts_non_red_filtered = counts_non_red_common %>% filter(symbol %in% keep_unimodal)

markers =
	counts_non_red_filtered %>%
	distinct(cell_type) %>%
	pull(cell_type) %>%
	gtools::permutations(n = length(.), r = 2, v = .) %>%
	as_tibble() %>%
	setNames(c("cell_type1", "cell_type2")) %>%
	mutate(contrast = sprintf("cell_type%s - cell_type%s", cell_type1, cell_type2)) %>%
	mutate(de =
				 	pmap(
				 		list(cell_type1, cell_type2, contrast),
				 		~ 	counts_non_red_filtered %>%
				 			filter(cell_type %in% c(..1, ..2)) %>%
				 			test_differential_abundance(~ 0 + cell_type, .contrasts = ..3, fill_missing_values = TRUE, action="get", omit_contrast_in_colnames = TRUE) %>%
				 			filter(logFC > 0) %>%
				 			arrange(FDR) %>%
				 			mutate(i = 1:n())
				 	)) %>%
	unnest(de)


markers %>% saveRDS("dev/markers_for_signature_workflow.rds")


(markers %>%

		# Filter best markers for monocytes
		filter(ct1=="monocyte" & i==1) %>%

		# Prettify contrasts for plotting
		unite(pair, c("ct1", "ct2"), remove = FALSE, sep = "\n") %>%

		# Reshape
		gather(which, cell_type, ct1, ct2) %>%
		distinct(pair,  symbol,   which, cell_type) %>%

		# Attach counts
		left_join(counts_non_red, by = c("symbol", "Cell type category")) %>%

		# Plot
		ggplot(aes(y = count_scaled + 1, x = cell_type, fill = cell_type)) +
		geom_boxplot() +
		facet_wrap(~pair+ symbol, scales ="free_x", nrow = 2) +
		scale_y_log10() +
		scale_fill_manual(values = getPalette(colourCount)) +
		my_theme +
		theme(	axis.text.x = element_text(angle = 30, hjust = 1))) %>%
	ggsave(	"dev/signature_p2.pdf",	plot = .,	useDingbats=FALSE,	units = c("mm"),	width = 183 ,	limitsize = FALSE)

# Plots and Study
(markers %>%
		filter(i < 6) %>%
		unite(pair, c("ct1", "ct2"), remove = FALSE, sep = "\n") %>%
		gather(which, cell_type, ct1, ct2) %>%

		distinct(symbol) %>%
		left_join(counts_non_red, by = c("symbol"))  %>%

		# Impute missng values
		tidybulk:::fill_NA_using_formula(~cell_type, sample, symbol, count_scaled) %>%

		reduce_dimensions(sample, symbol, count_scaled, method = "tSNE") %>%
		pivot_sample(sample) %>%
		ggplot(aes(x = `tSNE1`, y = `tSNE2`, color = cell_type)) +
		geom_point(size =2) +
		scale_color_manual(values = getPalette(colourCount)) +
		my_theme + theme(aspect.ratio=1)) %>%
	ggsave(	"dev/signature_p3.pdf",	plot = .,	useDingbats=FALSE,	units = c("mm"),	width = 183/2 ,	limitsize = FALSE)


