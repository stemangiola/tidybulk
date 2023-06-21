library(mygene)
library(furrr)
plan(multiprocess)

symbol = ensembl_symbol_mapping %>% distinct(hgnc_symbol) %>% pull(1)

res <- queryMany(symbol, scopes='symbol', fields=c('entrezgene', 'go'), species='human')

saveRDS(res, file="symbol_GO_object.rds", compress = "gzip")

 BP =
 	future_map2_dfr(
		1:dim(res),
		res[[1]],
		~ {
				res[.x, ]$go.BP[[1]] %>%
					as_tibble() %>%
					mutate(gene_symbol = .y) %>%
					nest(BP = -gene_symbol)
		}
	)

 CC =
 	future_map2_dfr(
 		1:dim(res),
 		res[[1]],
 		~ {
 			res[.x, ]$go.CC[[1]] %>%
 				as_tibble() %>%
 				mutate(gene_symbol = .y) %>%
 				nest(CC = -gene_symbol)
 		}
 	)

 MF =
 	future_map2_dfr(
 		1:dim(res),
 		res[[1]],
 		~ {
 			res[.x, ]$go.MF[[1]] %>%
 				as_tibble() %>%
 				mutate(gene_symbol = .y) %>%
 				nest(MF = -gene_symbol)
 		}
 	)

GO = BP %>% full_join(CC) %>% full_join(MF)

save(GO, file="data/GO.rda", compress = "gzip")

# library(GO.db)
# as.list(GOBPCHILDREN['GO:0006412'])

res_fly <- queryMany(
	flybaseIDs,
	#flybaseR::id.converter(flybaseIDs, symbols = TRUE, convert.into = "genes"),
	scopes='symbol',
	fields=c('entrezgene', 'go'),
	species='fruitfly'
)
