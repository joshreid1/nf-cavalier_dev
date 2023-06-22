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
  --exclude-benign-missense   Exclude missense variants that are predicted to be benign by sift, polyphen and ClinVar
  --sv                        Run in structural variant mode.
  --no-slides                 Don't create pptx output.
  --out=<f>                   Output file prefix [default: out].
  --family=<f>                Name of sample/family [default: Family].
  --genome=<f>                Reference genome for IGV snapshot [default: hg38].
  --caller=<f>                Name of variant caller [default: GATK].
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
  --include-sv-csv            Flag to retain all coding sequence variants (i.e. whole/partial exon)
  --sv-chr-exclude=<f>        Chromosomes to exclude from SV callsets [default: chrM,chrY,M,Y].
  --cache-dir=<f>             Cavalier cache dir [default: ~/.cavalier].
"
opts <- docopt(doc)
# print options
message('Using options:')
opts[names(opts) %>%
       keep(~str_detect(., '^[:alpha:]'))] %>%
  { class(.) <- c('list', 'docopt'); .} %>%
  print()

# set and check options
sample_bams <-
  setNames(str_extract(opts$sample_bam, '[^=]+$'),
           str_extract(opts$sample_bam, '^[^=]+'))

lists <- c(str_split(opts$gene_lists, ',', simplify = T))

invisible(
  assert_that(file.exists(opts$vcf),
              file.exists(opts$ped),
              all(file.exists(lists)),
              all(file.exists(sample_bams)),
              opts$genome %in% c('hg19', 'hg38')))

maf_dom <- as.numeric(opts$`--maf-dom`)
maf_rec <- as.numeric(opts$`--maf-rec`)
maf_comp_het <- as.numeric(opts$`--maf-comp-het`)
maf_de_novo <- as.numeric(opts$`--maf-de-novo`)
min_depth <- as.numeric(opts$`--min-depth`)
max_cohort_ac <- as.numeric(opts$`--max-cohort-ac`)
max_cohort_af <- as.numeric(opts$`--max-cohort-af`)
min_large_event <- as.numeric(opts$large_event)
min_impact <- ordered(opts$min_impact, c('MODIFIER', 'LOW', 'MODERATE', 'HIGH'))
sv_chr_exclude <-  c(str_split(opts$sv_chr_exclude, ',', simplify = T))

set_cavalier_opt(ref_genome = opts$genome)
set_cavalier_opt(cache_dir = opts$cache_dir)
# set_cavalier_opt(
#   singularity_img = '~/singularity_cache/bahlolab-cavalier-dev.img',
#   singularity_cmd = '/stornext/System/data/apps/singularity/singularity-3.7.3/bin/singularity')
insecure()

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


update_filter_stats <- function(data, set) {
  n <- n_distinct(data$variant_id)
  filter_stats <-
    get('filter_stats', envir = .GlobalEnv) %>% 
    add_row(set = set, n = n)
  assign('filter_stats', filter_stats, envir = .GlobalEnv)
  return(data)
}


if (!opts$sv) { # SNPS
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
  # load variants
  vars <-
    load_vcf(opts$vcf,
             samples = names(sample_bams),
             caller = opts$caller,
             annotater = 'VEP') %>%
    mutate(AN = pmax(0, AN - rowSums(mutate_all(genotype, ~str_count(., '[01]')))),
           AC = pmax(0, AC - rowSums(mutate_all(genotype, ~str_count(., '[1]')))),
           AF = if_else(AN > 0, AC / AN, 0))
  
  
  filter_stats <- tibble(set = 'all', n = n_distinct(vars$variant_id))
  
  # get candidates
  cand_vars <-
    vars %>%
    filter(is.na(af_gnomad) | af_gnomad <= maf_rec) %>% 
    update_filter_stats('gnomad_allele_freq') %>% 
    filter(AC <= max_cohort_ac,
           AF <= max_cohort_af) %>% 
    update_filter_stats('cohort_allele_freq') %>% 
    filter(impact >= min_impact) %>% 
    update_filter_stats('min_impact') %>% 
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
    update_filter_stats('benign_missense') %>%
    filter(gene %in% list_df$symbol) %>%
    update_filter_stats('gene_list') %>% 
    add_inheritance(ped_file = opts$ped,
                    af_column = 'af_gnomad',
                    af_compound_het = maf_comp_het,
                    af_dominant = maf_dom,
                    af_recessive = maf_rec,
                    af_de_novo = maf_dom,
                    min_depth = min_depth) %>%
    filter(!is.na(inheritance)) %>%
    update_filter_stats('inheritance') %>% 
    left_join(select(list_df, gene = symbol, panel_data = data),
              by = 'gene') %>%
    mutate(title = str_c(opts$family, ' - ', gene),
           cohort_AC_AF = str_c(AC, ' (', round(AF, 2), ')')) %>%
    arrange(gene, chrom, pos)
  
  
  # create slides
  if (!opts$no_slides) {
    create_slides(cand_vars,
                  output = str_c(opts$out, '.pptx'),
                  bam_files = sample_bams,
                  ped_file = opts$ped,
                  layout = layout,
                  var_info = c(cavalier::get_var_info(),
                               cohort_AC_AF = 'cohort_AC_AF'))
  } else {
    file.create(str_c(opts$out, '.pptx'))
  }
  
  cand_vars %>%
    mutate(set = 'SNP',
           family = opts$family,
           genotype = map_chr(seq_len(nrow(cand_vars)), function(i) {
             str_c(names(genotype), ':', genotype[i,], collapse = ';')
           }),
           baf = map_chr(seq_len(nrow(cand_vars)), function(i) {
             str_c(names(depth_ref), ':', depth_alt[i,], '/', depth_alt[i,] + depth_ref[i,], collapse = ';')
           })) %>% 
    select(set, family, gene, consequence, inheritance, id, hgvs_genomic, hgvs_coding, hgvs_protein,
           genotype, baf) %>%
    distinct() %>% 
    write_csv(str_c(opts$out, '.candidates.csv'))
  
  write_csv(filter_stats, str_c(opts$out, '.filter_stats.csv'))
  
} else { # SVS
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
  # load variants
  vars <-
    load_vcf(opts$vcf,
             samples = names(sample_bams),
             caller = 'manta',
             annotater = 'VEP',
             SVO = TRUE) %>%
    mutate(af_gnomad = pmax(svo_af, af_gnomad, na.rm = T)) %>%
    mutate(AN = pmax(0, AN - rowSums(mutate_all(genotype, ~str_count(., '[01]')))),
           AC = pmax(0, AC - rowSums(mutate_all(genotype, ~str_count(., '[1]')))),
           AF = if_else(AN > 0, AC / AN, 0)) %>% 
    filter(!chrom %in% sv_chr_exclude)
  
  filter_stats <- tibble(set = 'all', n = n_distinct(vars$variant_id))
  
  # get candidates
  cand_vars <-
    vars %>%
    filter(is.na(af_gnomad) | af_gnomad <= maf_rec) %>% 
    update_filter_stats('gnomad_allele_freq') %>% 
    filter(AC <= max_cohort_ac,
           AF <= max_cohort_af) %>%
    update_filter_stats('cohort_allele_freq')  %>% 
    (function(x) {
      # variants that intersect gene lists
      gl_vars <-
        filter(x,
               gene %in% list_df$symbol,
               impact >= min_impact |
                 (opts$include_sv_csv & str_detect(consequence, 'coding_sequence_variant')))
      # large event CNVs
      le_vars <-
        anti_join(x, gl_vars, 'variant_id') %>%
        filter(SVTYPE %in% c('DEL', 'DUP'),
               abs(SVLEN) >= min_large_event) %>%
        select(variant_id, chrom, pos, ref, alt, AF, AC, AN, END, SVTYPE, SVLEN, genotype, af_gnomad, id) %>%
        distinct() %>%
        mutate(gene = 'Large CNV')
      # combine
      bind_rows(gl_vars, le_vars)
    }) %>%
    update_filter_stats('min_impact') %>% 
    filter(gene %in% c('Large CNV', list_df$symbol)) %>% 
    update_filter_stats('gene_list') %>% 
    add_inheritance(ped_file = opts$ped,
                    af_column = 'af_gnomad',
                    af_compound_het = maf_comp_het,
                    af_dominant = maf_dom,
                    af_recessive = maf_rec,
                    min_depth = NULL) %>%
    filter(!is.na(inheritance)) %>%
    update_filter_stats('inheritance') %>% 
    annotate_gaps() %>%
    filter(is.na(gap_type)) %>%
    update_filter_stats('gaps') %>% 
    left_join(select(list_df, gene = symbol, panel_data = data),
              by = 'gene') %>%
    mutate(title = str_c(opts$family, ' - ', gene),
           cohort_AC_AF = str_c(AC, ' (', round(AF, 2), ')')) %>%
    group_by(variant_id) %>%
    mutate(other_genes = map(gene, ~setdiff(gene, .)) %>% map_chr(str_c, collapse = ',')) %>%
    ungroup() %>%
    arrange(gene, chrom, pos)
  
  # For SV, makes more sense to summarise across all affected genes
  
  # create slides
  if (!opts$no_slides) {
    create_slides(cand_vars,
                  output = str_c(opts$out, '.pptx'),
                  bam_files = sample_bams,
                  ped_file = opts$ped,
                  layout = layout,
                  var_info = c(cavalier::get_var_info(sv=TRUE),
                               cohort_AC_AF = 'cohort_AC_AF'))
  } else {
    file.create(str_c(opts$out, '.pptx'))
  }
  
  # write candidate info
  cand_vars %>%
    mutate(set = 'SV',
           family = opts$family) %>%
    select(set, family, gene, consequence, inheritance, id, SVTYPE, chrom, pos, END, SVLEN) %>%
    distinct() %>%
    write_csv(str_c(opts$out, '.candidates.csv'))
  
  write_csv(filter_stats, str_c(opts$out, '.filter_stats.csv'))
}

vid <- unique(cand_vars$variant_id)
write_lines(length(vid), str_c(opts$out, '.num_cand'))
out_vcf <- str_c(opts$out, '.candidates.vcf.gz')
if (length(vid)) {
  gds <- cavalier:::vcf_to_gds(opts$vcf)
  SeqArray::seqResetFilter(gds, verbose = F)
  SeqArray::seqSetFilter(gds, variant.id = vid)
  SeqArray::seqGDS2VCF(gds, out_vcf, verbose = F)
} else {
  file.create(out_vcf)
}