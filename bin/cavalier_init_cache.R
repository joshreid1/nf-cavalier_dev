#!/usr/bin/env Rscript

stopifnot(
  require(tidyverse), 
  require(cavalier)
)

opts_json      <- commandArgs(trailingOnly = TRUE)[1]
opts_json_out  <- commandArgs(trailingOnly = TRUE)[2]

opts <- jsonlite::fromJSON(opts_json)
do.call(set_cavalier_opt, opts)

build_caches(PanelApp = FALSE, Genes4Epilepsy = FALSE)

opts$hgnc_version    <- cavalier:::get_hgnc_version()
opts$hpo_version     <- cavalier:::get_hpo_version()
opts$mi_omim_version <- cavalier:::get_mi_omim_version()

cat(
  jsonlite::toJSON(opts),
  file = opts_json_out
)