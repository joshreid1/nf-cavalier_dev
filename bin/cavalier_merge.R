#!/usr/bin/env Rscript

stopifnot(require(tidyverse),
          require(docopt),
          require(cavalier))

doc <- "
Usage:
  cavalier_merge.R <output> <inputs>... [options]

Options:
  inputs                      Input filenames
  output                      Output filename
  --max-maf=<f>               Maximum MAF [default: 0.1].
  --gtex-rpkm=<f>             Path to GTEx_median_rpkm_file
  --omim-genemap2=<f>         Path to OMIM_genemap2_file
" 
opts <- docopt(doc)

qualvars <- map_df(opts$inputs, readRDS)
qualvars$homRef <- 
  qualvars %>% 
  select(contains('genotype')) %>%
  apply(1, function(x) sum(x %in% c('0/0', '0|0'), na.rm = T))

saveRDS(qualvars, opts$output)
