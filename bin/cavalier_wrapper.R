#!/usr/bin/env Rscript

stopifnot(require(tidyverse),
          require(magrittr),
          require(docopt),
          require(assertthat),
          require(cavalier))

doc <- "
Usage:
  cavalier_wrapper.R <vcf> <ped> <sample_bam>... [options]

Options:
  vcf                         Input VEP vcf file.
  ped                         Pedigree file
  sample_bams                 Sample names and bam files in format name=/path/to/bam.
  --out=<f>                   Output directory [default: out].
  --genome=<f>                Reference genome for IGV snapshot [default: hg19].
  --gene-lists=<f>            Comma serparated list of gene list names.
  --min-impact=<f>            Minimum VEP impact [default: MODERATE].
  --maf-dom=<f>               Maximum MAF for dominant [default: 0.0001].
  --maf-de-novo=<f>           Maximum MAF for dominant [default: 0.0001].
  --maf-rec=<f>               Maximum MAF for recessive [default: 0.01].
  --maf-comp-het=<f>          Maximum MAF for compound het [default: 0.01].
  --min-depth=<f>             Minimum depth for variant [default: 5].
  --max-cohort-ac=<f>         Maximum allele count within cohort [default: Inf].
  --max-cohort-af=<f>         Maximum allele frequency within cohort [default: Inf].
  --omim-genemap2=<f>         Path to OMIM_genemap2_file.
"
# opts <- ("AH017.subset.vcf.gz AH017.ped S35338_1=S35338_1.merged.bam S35339_1=S35339_1.merged.bam S35340_2=S35340_2.merged.bam S35341_2=S35341_2.merged.bam S36153_1=S36153_1.merged.bam \
#     --out AH017 \
#     --gene-lists custom_2_autoimmune_and_demyelinating.tsv \
#     --omim-genemap2 genemap2.txt \
#     --maf-dom 0.0001 \
#     --maf-rec 0.01 \
#     --maf-comp-het 0.01 \
#     --max-cohort-af 0.10") %>%
#  str_split('\\s+', simplify = T) %>%
#  str_trim() %>%
#   docopt(doc, .)

opts <- docopt(doc)
# print options
message('Using options:')
opts[names(opts) %>% 
       keep(~str_detect(., '^[:alpha:]'))] %>% 
  { class(.) <- c('list', 'docopt'); .} %>% 
  print()

set_cavalier_opt(singularity_img = '~/links/singularity_cache/jemunro-cavalier-dev.img',
                 singularity_cmd = '/stornext/System/data/apps/singularity/singularity-3.7.3/bin/singularity')

# set and check options
assert_that(file.exists(opts$vcf),
            file.exists(opts$ped))

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
    read_tsv(l, col_types = cols()) %T>% 
      { assert_that(all(c('list_id', 'list_name', 'gene') %in% names(.))) }
  }) %>% 
  rename(symbol = gene) %>% 
  mutate(symbol = hgnc_sym2sym(symbol)) %>% 
  nest(data = -symbol)
  
# load variants
vars <- load_vcf(opts$vcf, annot_source = 'VEP')

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

# get candidates and print slides
cand_vars <-
  vars %>% 
  filter(replace_na(sift != 'tolerated' |  polyphen != 'benign', TRUE),
         replace_na(!(is.na(sift) & polyphen == 'benign'), TRUE),
         replace_na(!(is.na(polyphen) & polyphen == 'tolerated'), TRUE),
         AC <= max_cohort_ac,
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
         cohort_AC_AF = str_c(AC, ' (', round(AF, 2), ')'),
         gene_ensembl = str_c(gene, ' (', ensembl_gene, ')')) %>% 
    create_slides(
      output = str_c(opts$out, '.pptx'),
      bam_files = sample_bams,
      ped_file = opts$ped,
      genemap2_file = opts$omim_genemap2, 
      layout = layout, 
      var_info = c(Gene = 'gene_ensembl',
                   cavalier:::get_var_info() %>% 
                     discard(~ . == 'gene'),
                   cohort_AC_AF = 'cohort_AC_AF'))
  