#!/usr/bin/env Rscript

stopifnot(require(docopt),
          require(cavalier),
          require(tidyverse))

doc <- "
Usage:
  cavalier_prep.R <vcf> <out_prefix> [options]

Options:
  vcf                         vcf(.gz) file with annotated variants
  out_prefix                  Prefix of output files
" 
opts <- docopt(doc)

MAF <- 0.1

region_include <- c("exonic", "splicing", "exonic;splicing")
change_exclude <- c("synonymous SNV")

samples <- 
  colnames((vcfR::read.vcfR(opts$vcf, nrows = 1))@gt)[-1] %>% 
  { setNames(as.list(.), .) }

vars <- (opts$vcf, samples)

qualvars <- 
  quality_filter_variants(vars) %>% 
  mutate(`MAF gnomAD genome` = suppressWarnings(as.numeric(`MAF gnomAD genome`)) %>% replace_na(0),
         `MAF gnomAD exome` = suppressWarnings(as.numeric(`MAF gnomAD exome`)) %>% replace_na(0)) %>% 
  filter(`MAF gnomAD genome` < MAF,
         region %in% region_include,
         !(change %in% change_exclude))

saveRDS(qualvars, paste0(opts$out_prefix, "_qualvars_filtered_exonic.rds"))
