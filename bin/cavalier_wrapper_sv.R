#!/usr/bin/env Rscript

stopifnot(require(tidyverse),
          require(magrittr),
          require(docopt),
          require(assertthat))
          # require(cavalier))

doc <- "
Usage:
  cavalier_wrapper.R <vcf> <ped> <sample_bam>... [options]

Options:
  vcf                         Input VEP vcf file.
  ped                         Pedigree file
  sample_bams                 Sample names and bam files in format name=/path/to/bam.
  --out=<f>                   Output directory [default: out].
  --genome=<f>                Reference genome for IGV snapshot [default: hg38].
  --gene-lists=<f>            Comma serparated list of gene list names.
  --min-impact=<f>            Minimum VEP impact [default: MODERATE].
  --maf-dom=<f>               Maximum MAF for dominant [default: 0.0001].
  --maf-de-novo=<f>           Maximum MAF for dominant [default: 0.0001].
  --maf-rec=<f>               Maximum MAF for recessive [default: 0.01].
  --maf-comp-het=<f>          Maximum MAF for compound het [default: 0.01].
  --min-depth=<f>             Minimum depth for variant [default: 5].
  --max-cohort-ac=<f>         Maximum allele count within cohort [default: Inf].
  --max-cohort-af=<f>         Maximum allele frequency within cohort [default: Inf].
"
opts <- ("AH029.subset.vcf.gz AH029.ped S35565_1=S35565_1.merged.bam \
    --out AH029 \
    --genome hg38 \
    --gene-lists PAA_3531_v1.52.tsv \
    --maf-dom 0.0001 \
    --maf-de-novo 0.0001 \
    --maf-rec 0.01 \
    --maf-comp-het 0.01 \
    --max-cohort-af 0.10 \
    --min-impact MODERATE") %>%
  str_split('\\s+', simplify = T) %>%
  str_trim() %>%
  docopt(doc, .)
# opts <- docopt(doc)
# print options
message('Using options:')
opts[names(opts) %>% 
       keep(~str_detect(., '^[:alpha:]'))] %>% 
  { class(.) <- c('list', 'docopt'); .} %>% 
  print()

# set and check options
assert_that(file.exists(opts$vcf),
            file.exists(opts$ped),
            opts$genome %in% c('hg19', 'hg38'))

set_cavalier_opt(ref_genome = opts$genome)
insecure()

sample_bams <- 
  setNames(str_extract(opts$sample_bam, '[^=]+$'),
           str_extract(opts$sample_bam, '^[^=]+'))

assert_that(all(file.exists(sample_bams)))

maf_dom <- as.numeric(opts$`--maf-dom`)
maf_rec <- as.numeric(opts$`--maf-rec`)
maf_comp_het <- as.numeric(opts$`--maf-comp-het`)
maf_de_novo <- as.numeric(opts$`--maf-de-novo`)
min_depth <- as.numeric(opts$`--min-depth`)
max_cohort_ac <- as.numeric(opts$`--max-cohort-ac`)
max_cohort_af <- as.numeric(opts$`--max-cohort-af`)

min_impact <- ordered(opts$min_impact,
                      c('MODIFIER', 'LOW', 'MODERATE', 'HIGH'))

lists <- c(str_split(opts$gene_lists, ',', simplify = T))
assert_that(all(file.exists(lists)))

list_df <-
  map_df(lists, function(l) {
    read_tsv(l, col_types = cols(.default = "c")) %T>% 
      { assert_that(all(c('list_id', 'list_name', 'gene') %in% names(.))) }
  }) %>% 
  mutate(list_name_url = case_when(
    str_starts(list_id, 'PAA:') ~ str_c('https://panelapp.agha.umccr.org/panels/', 
                                        str_extract(list_id, '(?<=PAA:)\\d+')),
    str_starts(list_id, 'PAE:') ~ str_c('https://panelapp.genomicsengland.co.uk/panels/', 
                                        str_extract(list_id, '(?<=PAE:)\\d+')),
    str_starts(list_id, 'HP:') ~ str_c('https://hpo.jax.org/app/browse/term/', 
                                       list_id)
  )) %>% 
  rename(symbol = gene) %>% 
  filter(!is.na(symbol)) %>% 
  mutate(symbol = hgnc_sym2sym(symbol),
         list_id_url = list_name_url) %>%
  distinct() %>%
  group_by(list_id, list_name, list_id_url, list_name_url, symbol) %>%
  (function(x) 
    `if`(n_groups(x) < nrow(x),
         summarise(x,
                   across(everything(), ~ str_c(sort(unique(na.omit(.))), collapse = '; ')),
                   .groups = 'drop'),
         ungroup(x))) %>%
  nest(data = -symbol)

    # load variants
vars <- load_vcf(opts$vcf,
                 caller = 'manta-jasmine',
                 annotater = 'VEP')

# slide_layout
layout <-
  `if`(length(sample_bams) == 1,
       slide_layout(c('var_info', 'igv', 'gtex'),
                    c('omim', 'custom_panel_data'),
                    heights = c(5,2)),
       bind_rows(
         slide_layout(c('var_info', 'pedigree', 'gtex'),
                      c('omim', 'custom_panel_data'),
                      heights = c(5,2)),
         slide_layout(c('igv'),
                      slide_num = 2L))
  )

# get candidates 
cand_vars <-
  vars %>% 
  (function(x) {
    if (opts$exclude_benign_missense) {
      x %>% 
        mutate(idx = seq_along(variant_id)) %>% 
        filter(consequence == 'missense_variant') %>% 
        select(idx, sift, polyphen, clin_sig) %>% 
        mutate(sift = sift > ordered('tolerated', levels(sift)),
               polyphen = polyphen > ordered('benign', levels(polyphen)),
               clin_sig = clin_sig > ordered('likely_benign ', levels(clin_sig)),
               all_na = is.na(sift) & is.na(polyphen) & is.na(clin_sig),
               any_pass = replace_na(sift | polyphen | clin_sig, FALSE) | all_na) %>% 
        filter(any_pass) %>% 
        pull(idx) %>% 
        { filter(x, consequence != 'missense_variant' | (seq_along(variant_id) %in% .)) }
    } else { x }
  }) %>%
  filter(AC <= max_cohort_ac,
         AF <= max_cohort_af,
         impact >= min_impact,
         gene %in% list_df$symbol) %>% 
  add_inheritance(ped_file = opts$ped,
                  af_column = 'af_gnomad',
                  af_compound_het = maf_comp_het,
                  af_dominant = maf_dom,
                  af_recessive = maf_rec,
                  min_depth = min_depth) %>%
  filter(!is.na(inheritance)) %>%
  left_join(select(list_df, gene = symbol, panel_data = data),
            by = 'gene') %>%
  mutate(title = str_c(opts$out, ' - ', gene),
         cohort_AC_AF = str_c(AC, ' (', round(AF, 2), ')'))

# create slides
create_slides(cand_vars,
              output = str_c(opts$out, '.pptx'),
              bam_files = sample_bams,
              ped_file = opts$ped,
              layout = layout,
              var_info = c(cavalier::get_var_info(),
                           cohort_AC_AF = 'cohort_AC_AF'))
