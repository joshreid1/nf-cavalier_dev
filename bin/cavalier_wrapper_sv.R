#!/usr/bin/env Rscript

stopifnot(require(tidyverse),
          require(magrittr),
          require(docopt),
          require(assertthat),
          require(cavalier))

doc <- "
Usage:
  cavalier_wrapper_sv.R <vcf> <ped> <sample_bam>... [options]

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
  --large-event=<f>           Minimum size for a large CNV event to be kept regardless of gene intersections [default: 1e6]
"
# opts <- ("AH042-M.subset.vcf.gz AH042-M.ped S35669_1=S35669_1.merged.bam \
#     --out AH042-M \
#     --genome hg38 \
#     --gene-lists PAA_137_v0.9939.tsv \
#     --maf-dom 0.0001 \
#     --maf-de-novo 0.0001 \
#     --maf-rec 0.01 \
#     --maf-comp-het 0.01 \
#     --max-cohort-af 0.10 \
#     --min-impact MODERATE") %>%
#   str_split('\\s+', simplify = T) %>%
#   str_trim() %>%
#   docopt(doc, .)
opts <- docopt(doc)

csq_keep <- 
  str_c(
    # "5_prime_UTR_variant",
        "coding_sequence_variant",
        # "intron_variant",
        sep = '|')

# print options
message('Using options:')
opts[names(opts) %>% 
       keep(~str_detect(., '^[:alpha:]'))] %>% 
  { class(.) <- c('list', 'docopt'); .} %>% 
  print()

# set and check options
invisible(
  assert_that(file.exists(opts$vcf),
              file.exists(opts$ped),
              opts$genome %in% c('hg19', 'hg38')))

set_cavalier_opt(ref_genome = opts$genome)
insecure()

sample_bams <- 
  setNames(str_extract(opts$sample_bam, '[^=]+$'),
           str_extract(opts$sample_bam, '^[^=]+'))

invisible(assert_that(all(file.exists(sample_bams))))

maf_dom <- as.numeric(opts$`--maf-dom`)
maf_rec <- as.numeric(opts$`--maf-rec`)
maf_comp_het <- as.numeric(opts$`--maf-comp-het`)
maf_de_novo <- as.numeric(opts$`--maf-de-novo`)
min_depth <- as.numeric(opts$`--min-depth`)
max_cohort_ac <- as.numeric(opts$`--max-cohort-ac`)
max_cohort_af <- as.numeric(opts$`--max-cohort-af`)
min_large_event <- as.numeric(opts$large_event)

min_impact <- ordered(opts$min_impact,
                      c('MODIFIER', 'LOW', 'MODERATE', 'HIGH'))

lists <- c(str_split(opts$gene_lists, ',', simplify = T))
invisible(assert_that(all(file.exists(lists))))

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
vars <-
  load_vcf(opts$vcf,
           caller = 'manta',
           annotater = 'VEP',
           SVO = TRUE) %>% 
  mutate(af_gnomad = pmax(svo_af, af_gnomad, na.rm = T))

# add large event annotation


# slide_layout
layout <-
  `if`(length(sample_bams) == 1,
       slide_layout(c('var_info', 'gtex'),
                    c('omim', 'custom_panel_data'),
                    heights = c(5,2)),
       slide_layout(c('var_info', 'pedigree', 'gtex'),
                    c('omim', 'custom_panel_data'),
                    heights = c(5,2))
  )

# get candidates 
cand_vars <-
  vars %>% 
  filter(AC <= max_cohort_ac,
         AF <= max_cohort_af) %>% 
  (function(x) {
    # variants that intersect gene lists
    gl_vars <-
      filter(x, 
           gene %in% list_df$symbol,
           impact > min_impact | str_detect(consequence, csq_keep))
    # large event CNVs
    le_vars <-
      anti_join(x, gl_vars, 'variant_id') %>% 
      filter(SVTYPE %in% c('DEL', 'DUP'),
             abs(SVLEN) >= min_large_event) %>% 
      select(variant_id, chrom, pos, ref, alt, AF, AC, AN, END, SVTYPE, SVLEN, genotype, af_gnomad) %>% 
      distinct() %>% 
      mutate(gene = 'Large CNV')
    # combine
    bind_rows(gl_vars, le_vars)
  }) %>% 
  add_inheritance(ped_file = opts$ped,
                  af_column = 'af_gnomad',
                  af_compound_het = maf_comp_het,
                  af_dominant = maf_dom,
                  af_recessive = maf_rec,
                  min_depth = NULL) %>%
  filter(!is.na(inheritance)) %>%
  left_join(select(list_df, gene = symbol, panel_data = data),
            by = 'gene') %>%
  mutate(title = str_c(opts$out, ' - ', gene),
         cohort_AC_AF = str_c(AC, ' (', round(AF, 2), ')')) %>% 
  group_by(variant_id) %>% 
  mutate(other_genes = map(gene, ~setdiff(gene, .)) %>% map_chr(str_c, collapse = ',')) %>% 
  ungroup()

# For SV, makes more sense to summarise across all affected genes

# create slides
create_slides(cand_vars,
              output = str_c(opts$out, '.pptx'),
              bam_files = sample_bams,
              ped_file = opts$ped,
              layout = layout,
              var_info = c(cavalier::get_var_info(sv=TRUE),
                           cohort_AC_AF = 'cohort_AC_AF'))

cand_vars %>% 
  mutate(family = opts$out) %>% 
  select(family, chrom, pos, SVTYPE, SVLEN, END) %>% 
  distinct() %>% 
  write_csv(str_c(opts$out, '.candidates.csv'))

vid <- unique(cand_vars$variant_id)
write_lines(length(vid), str_c(opts$out, '.num_cand'))
out_vcf <- str_c(opts$out, '.candidates.vcf.gz')
if (length(vid)) {
  gds <- cavalier:::vcf_to_gds(opts$vcf)
  SeqArray::seqSetFilter(gds, variant.id = unique(cand_vars$variant_id))
  SeqArray::seqGDS2VCF(gds, str_c(opts$out, '.candidates.vcf.gz'),
                       verbose = F)
} else {
  file.create(out_vcf)
}

