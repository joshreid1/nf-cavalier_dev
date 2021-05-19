#!/usr/bin/env Rscript

stopifnot(require(tidyverse),
          require(docopt),
          require(cavalier))

doc <- "
Usage:
  cavalier_prep.R <vcf> <out_prefix> [options]

Options:
  vcf                         vcf(.gz) file with annotated variants.
  out_prefix                  Prefix of output files.
  --max-maf=<f>               Maximum MAF [default: 0.1].
  --gtex-rpkm=<f>             Path to GTEx_median_rpkm_file
  --omim-genemap2=<f>         Path to OMIM_genemap2_file
" 
opts <- docopt(doc)

MAF <- as.numeric(opts$max_maf)
# GTEx_median_rpkm_file <- opts$gtex_rpkm
# OMIM_genemap2_file <- opts$omim_genemap2

region_include <- c("exonic", "splicing", "exonic;splicing")
change_exclude <- c("synonymous SNV")

samples <- 
  colnames((vcfR::read.vcfR(opts$vcf, nrows = 1))@gt)[-1] %>% 
  { setNames(as.list(.), .) }

vars <- load_annovar_vcf(opts$vcf, samples)

qualvars <- 
  quality_filter_variants(vars) %>% 
  mutate(`MAF gnomAD genome` = suppressWarnings(as.numeric(`MAF gnomAD genome`)) %>% replace_na(0),
         `MAF gnomAD exome` = suppressWarnings(as.numeric(qualvars$`MAF gnomAD exome`)) %>% replace_na(0)) %>% 
  filter(`MAF gnomAD genome` < MAF,                       
         region %in% region_include,
         !(change %in% change_exclude))

saveRDS(qualvars, paste0(opts$out_prefix, "_qualvars_filtered_exonic.rds"))
