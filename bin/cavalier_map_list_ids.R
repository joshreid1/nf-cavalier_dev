#!/usr/bin/env Rscript

stopifnot(
  require(tidyverse), 
  require(cavalier)
)

input  <- commandArgs(trailingOnly = TRUE)[1]
output <- commandArgs(trailingOnly = TRUE)[2]
opts   <- commandArgs(trailingOnly = TRUE)[3]

set_options_from_json(opts)

list_df <- read_tsv(input, col_types = cols(.default = "c"))
# check required columns are present
req_cols <- c('symbol', 'gene', 'hgnc_id', 'ensembl_gene_id', 'entrez_id')
stopifnot(
  length(intersect(colnames(list_df), req_cols)) > 0
)
# map gene identifiers to ensembl_gene_id
list_df <-
  list_df %>% 
  bind_rows(req_cols %>% setNames(.,.) %>% map_dfc(~ character())) %>% 
  mutate(
    entrez_id = as.integer(entrez_id),
    symbol = coalesce(symbol, gene),
    hgnc_id = coalesce(
      hgnc_id,
      hgnc_ensembl2id(ensembl_gene_id),
      hgnc_entrez2id(entrez_id),
      hgnc_sym2id(symbol),
    ),
    symbol = coalesce(symbol, hgnc_id2sym(hgnc_id)),
    ensembl_gene_id = coalesce(ensembl_gene_id, hgnc_id2ensembl(hgnc_id)),
    entrez_id = coalesce(entrez_id, hgnc_id2entrez(hgnc_id)),
  ) %>% 
  select(-gene)
  
write_tsv(list_df, output)
