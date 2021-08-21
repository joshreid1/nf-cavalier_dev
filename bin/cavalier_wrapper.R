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
  --genome=<f>                Reference genome for IGV snapshot [default: hg19].
  --gene-lists=<f>            Comma serparated list of gene list names.
  --maf-dom=<f>               Maximum MAF for dominant [default: 0.0001].
  --maf-de-novo=<f>           Maximum MAF for dominant [default: 0.0001].
  --maf-rec=<f>               Maximum MAF for recessive [default: 0.01].
  --maf-comp-het=<f>          Maximum MAF for compound het [default: 0.01].
  --min-depth=<f>             Minimum depth for variant [default: 5].
  --max-cohort-ac=<f>         Maximum allele count within cohort [default: Inf].
  --max-cohort-af=<f>         Maximum allele frequency within cohort [default: Inf].
  --omim-genemap2=<f>         Path to OMIM_genemap2_file.
"
# opts <- ("AH050.subset.vcf.gz AH050.ped ST732_3=ST732_3.merged.bam ST733_3=ST733_3.merged.bam ST734_3=ST734_3.merged.bam ST934_3=ST934_3.merged.bam \
#     --out AH050 \
#     --gene-lists agha_202_genetic_epilepsy.tsv,agha_250_intellectual_disability_syndromic_and_non_syndromic.tsv,ge_285_intellectual_disability.tsv,ge_402_genetic_epilepsy_syndromes.tsv,ge_486_paediatric_disorders.tsv,ge_490_hypotonic_infant.tsv,ge_496_white_matter_disorders_childhood_onset.tsv \
#     --omim-genemap2 genemap2.txt \
#     --maf-dom 0.0001 \
#     --maf-rec 0.01 \
#     --maf-comp-het 0.01 \
#     --max-cohort-af 0.10") %>%
#  str_split('\\s+', simplify = T) %>%
#  str_trim() %>%
#   docopt(doc, .)

# opts <- docopt(doc)
# print options
message('Using options:')
opts[names(opts) %>% 
       keep(~str_detect(., '^[:alpha:]'))] %>% 
  { class(.) <- c('list', 'docopt'); .} %>% 
  print()

# set and check options
assert_that(file.exists(opts$vcf),
            file.exists(opts$ped))

sample_ids <- str_extract(opts$sample_bam, '^[^=]+')
sample_bams <- str_extract(opts$sample_bam, '[^=]+$')
assert_that(all(file.exists(sample_bams)))

maf_dom <- as.numeric(opts$`--maf-dom`)
maf_rec <- as.numeric(opts$`--maf-rec`)
maf_comp_het <- as.numeric(opts$`--maf-comp-het`)
maf_de_novo <- as.numeric(opts$`--maf-de-novo`)
min_depth <- as.numeric(opts$`--min-depth`)
max_cohort_ac <- as.numeric(opts$`--max-cohort-ac`)
max_cohort_af <- as.numeric(opts$`--max-cohort-af`)

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

cand_vars <-
  vars %>% 
  filter(!(sift == 'tolerated' &  polyphen == 'benign'),
         !(polyphen == 'benign' & is.na(sift)),
         !(sift == 'tolerated' & is.na(polyphen)),
         gene %in% list_df$symbol) %>%
  add_inheritance(ped_file = opts$ped) %>% 
  filter(!is.na(inheritance))



# create cavalier output if any variants remain
dir.create(opts$out, recursive = TRUE, showWarnings = FALSE)
if (nrow(candvars)) {
  output_cols <- c('Inheritance', "Variant", "Amino acid", "change", "Depth (R,A)", "Cohort AF",
                   "GnomAD AF", "SIFT", "Polyphen2", "Grantham", "RVIS", "GeVIR")
  candvars %>%
  create_igv_snapshots(sample_bam, "hg19", vcfs = opts$vcf) %>%
  #   create_igv_snapshots(sample_bam,  "hg19",
  #                        vcfs = opts$vcf,
  #                        # overwrite = TRUE,
  #                        singularity_img = '~/links/singularity_cache/jemunro-cavalier-dev.img',
  #                        singularity_bin = '/stornext/System/data/apps/singularity/singularity-3.7.3/bin/singularity') %>%
  #   mutate(panel_data = map(seq_along(gene), ~tibble())) %>% 
    mutate(Inheritance = `inheritance model`,
           Variant = str_c(chromosome, ':', position, ':', reference, '>', alternate),
           `Depth (R,A)` = `proband depth (R,A)`,
           `Amino acid` = str_replace(Amino_acids, '/', '>'),
           Polyphen2 = str_c(Polyphen2, ' (', Polyphen2_score, ')'),
           SIFT = str_c(SIFT, ' (', SIFT_score, ')'),
           title = str_c('Sample: ', sample_id, ', Gene: ', gene),
           `Cohort AF` = str_c(round(AF, 4), ' (', AC, ')')
    ) %>%
    rename(`GnomAD AF` = MAF_gnomAD) %>% 
    as.data.frame() %>%
    create_cavalier_output(opts$out, sampleID, output_cols,
                           hide_missing_igv = TRUE,
                           layout = "individual",
                           genemap2 = opts$`--omim-genemap2`,
                           GTEx_median_rpkm = opts$`--gtex-rpkm`,
                           title_col = 'title', 
                           add_data_col = 'panel_data',
                           output = c('html', 'ppt'))
}

