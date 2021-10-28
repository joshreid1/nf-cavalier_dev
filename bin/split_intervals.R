#!/usr/bin/env Rscript

stopifnot(require(docopt),
          require(rtracklayer),
          require(GenomicRanges),
          require(tidyverse))

doc <- "
Usage:
  split_intervals.R <ref_fai> <gaps_bed> <n> <pref>

Options:
  ref_fai       Reference fasta fai file
  gaps_bed      Bed file with Gap positions
  n             Maximum/Ideal number of regions to split
  pref          Output file prefix
"

opts <- docopt(doc)

n <- as.integer(opts$n)

read_tsv(opts$ref_fai, 
         col_names = c('chrom', 'len', 'offset', 'lb', 'lw'),
         col_types = c('cidii')) %>% 
  with(GRanges(chrom,IRanges(1, len)))  %>% 
  { suppressWarnings(GenomicRanges::setdiff(., import.bed(opts$gaps_bed))) } %>%
  GenomicRanges::reduce() %>% 
  as_tibble() %>% 
  mutate(tot = cumsum(as.numeric(width)),
         set = 1 + tot %/% (ceiling(last(tot) / n)),
         set = str_replace_all(
           format(as.integer(as.factor(set))), '\\s', '0')) %>%
  nest(data = -set) %>% 
  pwalk(function(set, data) {
    select(data, 1:3) %>% 
      write_tsv(str_c(opts$pref, '-', set, '.intervals.tsv'),
                col_names = FALSE)
  })
