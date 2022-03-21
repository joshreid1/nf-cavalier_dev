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
  --out=<f>                   Output file prefix [default: out].
  --family=<f>                Name of sample/family [default: Family].
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
#   singularity_img = '~/links/singularity_cache/bahlolab-cavalier-dev.img',
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
             caller = 'GATK',
             annotater = 'VEP') %>% 
    mutate(AN = AN - rowSums(mutate_all(genotype, ~str_count(., '[01]'))),
           AC = AC - rowSums(mutate_all(genotype, ~str_count(., '[1]'))),
           AF = AC / AN)
  
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
    mutate(title = str_c(opts$family, ' - ', gene),
           cohort_AC_AF = str_c(AC, ' (', round(AF, 2), ')')) %>% 
    arrange(gene, chrom, pos)
  
  
  # create slides
  create_slides(cand_vars,
                output = str_c(opts$out, '.pptx'),
                bam_files = sample_bams,
                ped_file = opts$ped,
                layout = layout,
                var_info = c(cavalier::get_var_info(),
                             cohort_AC_AF = 'cohort_AC_AF'))
  
  # TBD: write candidate info
  cand_vars %>%
    mutate(set = 'SNP',
           family = opts$family) %>%
    select(set, family, gene, consequence, inheritance, id, hgvs_genomic, hgvs_protein) %>%
    distinct() %>%
    write_csv(str_c(opts$out, '.candidates.csv'))

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
    mutate(AN = AN - rowSums(mutate_all(genotype, ~str_count(., '[01]'))),
           AC = AC - rowSums(mutate_all(genotype, ~str_count(., '[1]'))),
           AF = AC / AN)
  
  # get candidates 
  cand_vars <-
    vars %>% 
    filter(!chrom %in% sv_chr_exclude) %>% 
    filter(AC <= max_cohort_ac,
           AF <= max_cohort_af) %>% 
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
    annotate_gaps() %>% 
    filter(is.na(gap_type)) %>% 
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
  create_slides(cand_vars,
                output = str_c(opts$out, '.pptx'),
                bam_files = sample_bams,
                ped_file = opts$ped,
                layout = layout,
                var_info = c(cavalier::get_var_info(sv=TRUE),
                             cohort_AC_AF = 'cohort_AC_AF'))
  
  # write candidate info
  cand_vars %>% 
    mutate(set = 'SV',
           family = opts$family) %>%
    select(set, family, gene, consequence, inheritance, id, SVTYPE, chrom, pos, END, SVLEN) %>%
    distinct() %>% 
    write_csv(str_c(opts$out, '.candidates.csv'))
  
}

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
