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
  --caller=<f>                Name of variant caller [default: GATK].
  --gene-lists=<f>            Comma serparated list of gene list names.
  --min-impact=<f>            Minimum VEP impact [default: MODERATE].
  --vcfanno_config=<f>        vcfanno config.
  --filter_by_annotation=<f>  Filter variants by column name.
  --min_value=<f>             Filter variants by min column value.
  --max_value=<f>             Filter variants by max column value.
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
  --cavalier-options=<f>      Additional Cavalier options in json format
"
opts <- docopt(doc)
set_options_from_json(opts$cavalier_options)
# print options
message('Using options:')
opts[names(opts) %>%
       keep(~str_detect(., '^[:alpha:]'))] %>%
  { class(.) <- c('list', 'docopt'); .} %>%
  print()

# for debugging
saveRDS(opts, '.opts.rds')
saveRDS(as.list(cavalier:::cavalier_opts), '.cavalier_opts.rds')

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
              get_cavalier_opt("ref_genome") %in% c('hg19', 'hg38')))

maf_dom <- as.numeric(opts$`--maf-dom`)
maf_rec <- as.numeric(opts$`--maf-rec`)
maf_comp_het <- as.numeric(opts$`--maf-comp-het`)
maf_de_novo <- as.numeric(opts$`--maf-de-novo`)
min_depth <- as.numeric(opts$`--min-depth`)
max_cohort_ac <- as.numeric(opts$`--max-cohort-ac`)
max_cohort_af <- as.numeric(opts$`--max-cohort-af`)
vcfanno_config <- as.character(opts$`--vcfanno_config`)
filter_by_annotation <- as.character(opts$`--filter_by_annotation`)
min_filter_value <- as.numeric(opts$`--min_value`)
max_filter_value <- as.numeric(opts$`--max_value`)
min_large_event <- as.numeric(opts$large_event)
min_impact <- ordered(opts$min_impact, c('MODIFIER', 'LOW', 'MODERATE', 'HIGH'))
sv_chr_exclude <-  c(str_split(opts$sv_chr_exclude, ',', simplify = T))

list_df <-
  map_df(lists, function(fn) {
    list_df <- read_tsv(fn, col_types = cols(.default = "c"))
    assert_that('ensembl_gene_id' %in% colnames(list_df))
    
    list_df %>% 
      select(
        ensembl_gene_id,
        any_of(c('list_id', 'list_name', 'list_version', 'inheritance')),
        any_of(starts_with('meta_'))
      ) %>% 
      rename_with(~ str_replace(., '^meta_', '')) %>% 
      filter(!is.na(ensembl_gene_id)) %>% 
      distinct()
  }) %>%
  bind_rows(tibble(list_id = character(), list_name = character())) %>% 
  mutate(list_name_url = case_when(
    str_starts(list_id, 'PAA:') ~ str_c('https://panelapp.agha.umccr.org/panels/',
                                        str_extract(list_id, '(?<=PAA:)\\d+')),
    str_starts(list_id, 'PAE:') ~ str_c('https://panelapp.genomicsengland.co.uk/panels/',
                                        str_extract(list_id, '(?<=PAE:)\\d+')),
    str_starts(list_id, 'HP:') ~ str_c('https://hpo.jax.org/app/browse/term/',
                                       list_id),
    str_starts(list_id, 'G4E:') ~ 'https://bahlolab.github.io/Genes4Epilepsy/'
  )) %>%
  mutate(list_id_url = list_name_url) %>%
  nest(panel_data = -ensembl_gene_id)

filter_stats <- tibble(set = character(), n = integer())
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
  
  # Extract list of annotations from vcfanno config file to add on the slides
  add_anno = c()
  if (vcfanno_config != 'NULL') {
    file_content <- readLines(vcfanno_config)
    names_line <- file_content[grep("^names", file_content)] # TODO: update this line to make it more robust
    matches <- gregexpr('"([^"]*)"', names_line)
    extracted_strings <- regmatches(names_line, matches)
    add_anno = unlist(lapply(extracted_strings, function(x) gsub('"', '', x)))
  }
  
  # load variants
  vars <-
    load_vcf(opts$vcf,
             samples = names(sample_bams),
             caller = opts$caller,
             annotater = 'VEP') %>%
    mutate(AN = pmax(0, AN - rowSums(mutate_all(genotype, ~str_count(., '[01]')))),
           AC = pmax(0, AC - rowSums(mutate_all(genotype, ~str_count(., '[1]')))),
           AF = if_else(AN > 0, AC / AN, 0))
  additional_filter <- function(data, filter_column, min_value, max_value) {
    if (filter_column != 'NULL' && filter_column %in% colnames(data)) {
      if (max_value == 'Inf') {
        max_value = Inf
      }
      if (min_value <= max_value) {
        data <- data %>%
        filter((!!sym(filter_column) <= max_value & !!sym(filter_column) >= min_value) | is.na(!!sym(filter_column)))
      }
    }
    if (!(filter_column %in% colnames(data))) {
      message(paste0('Missing annotation: ', filter_column))
      message(paste0('Available annotations: ', colnames(data)))
    }
    return(data)
  }
  # get candidates
  cand_vars <-
    vars %>%
    update_filter_stats('all') %>%
    filter(is.na(af_gnomad) | af_gnomad <= maf_rec) %>% 
    additional_filter(filter_column = filter_by_annotation, min = min_filter_value, max = max_filter_value) %>%
    update_filter_stats('gnomad_allele_freq') %>% 
    filter(AC <= max_cohort_ac,
           AF <= max_cohort_af) %>% 
    update_filter_stats('cohort_allele_freq') %>% 
    filter(impact >= min_impact) %>% 
    update_filter_stats('min_impact') %>% 
    # bening_missense 
    mutate(
      non_benign_sift =  sift > ordered('tolerated', levels(sift)),
      non_benign_polyphen = polyphen > ordered('benign', levels(polyphen)),
      non_benign_clin_sig = clin_sig > ordered('likely_benign ', levels(clin_sig)),
      non_benign_all_na = is.na(non_benign_sift) & is.na(non_benign_polyphen) & is.na(non_benign_clin_sig),
      benign_missense = 
        consequence == 'missense_variant' &
        !(replace_na(non_benign_sift | non_benign_polyphen | non_benign_clin_sig, FALSE) | non_benign_all_na)
    ) %>% 
    select(-starts_with('non_benign')) %>% 
    filter(!(opts$exclude_benign_missense & benign_missense)) %>% 
    update_filter_stats('benign_missense') %>%
    semi_join(list_df, by = 'ensembl_gene_id') %>% 
    update_filter_stats('gene_list') %>% 
    add_inheritance(ped_file = opts$ped,
                    af_column = 'af_gnomad',
                    gene_column = 'ensembl_gene_id',
                    uid_column = 'hgvs_genomic',
                    af_compound_het = maf_comp_het,
                    af_dominant = maf_dom,
                    af_recessive = maf_rec,
                    af_de_novo = maf_dom,
                    min_depth = min_depth)  %>% 
    filter(!is.na(inheritance)) %>%
    update_filter_stats('inheritance') %>% 
    left_join(list_df, by = 'ensembl_gene_id') %>%
    mutate(title = str_c(opts$family, ' - ', symbol),
           cohort_AC_AF = str_c(AC, ' (', round(AF, 2), ')')) %>%
    arrange(symbol, chrom, pos)
  
  
  # create slides
  if (!opts$no_slides) {
    create_slides(variants = cand_vars,
                  output = str_c(opts$out, '.pptx'),
                  bam_files = sample_bams,
                  ped_file = opts$ped,
                  layout = layout,
                  var_info = c(cavalier::get_var_info(),
                               cohort_AC_AF = 'cohort_AC_AF',
                               add_anno))
    
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
    select(set, family, symbol, consequence, inheritance, id, hgvs_genomic, hgvs_coding, hgvs_protein,
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

  # get candidates
  cand_vars <-
    vars %>%
    update_filter_stats('all') %>% 
    filter(is.na(af_gnomad) | af_gnomad <= maf_rec) %>% 
    update_filter_stats('gnomad_allele_freq') %>% 
    filter(AC <= max_cohort_ac,
           AF <= max_cohort_af) %>%
    update_filter_stats('cohort_allele_freq') %>% 
    (function(data) {
      bind_rows(
        # variants that intersect gene lists
        data %>% 
          semi_join(list_df, by = 'ensembl_gene_id'),
        # large event CNVs
        data %>% 
          filter(SVTYPE %in% c('DEL', 'DUP'),
                 abs(SVLEN) >= min_large_event) %>%
          select(variant_id, chrom, pos, ref, alt, AF, AC, AN, END, SVTYPE, SVLEN, genotype, af_gnomad, id, hgvs_genomic) %>%
          distinct() %>%
          mutate(symbol = 'Large CNV')
      )
    }) %>% 
    update_filter_stats('gene_list') %>% 
    filter(
      impact >= min_impact |
        (opts$include_sv_csv & str_detect(consequence, 'coding_sequence_variant') |
           symbol == 'Large CNV'
        )
    ) %>% 
    update_filter_stats('min_impact') %>% 
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
    left_join(list_df, by = 'ensembl_gene_id') %>%
    mutate(title = str_c(opts$family, ' - ', symbol),
           cohort_AC_AF = str_c(AC, ' (', round(AF, 2), ')')) %>% 
    group_by(variant_id) %>%
    mutate(other_genes = map(symbol, ~ setdiff(symbol, .)) %>% map_chr(str_c, collapse = ',')) %>%
    ungroup() %>%
    arrange(symbol, chrom, pos)
  
  # For SV, makes more sense to summarise across all affected genes
  
  # create slides
  if (!opts$no_slides) {
    create_slides(variants = cand_vars,
                  output = str_c(opts$out, '.pptx'),
                  bam_files = sample_bams,
                  ped_file = opts$ped,
                  layout = layout,
                  is_sv = TRUE,
                  var_info = c(cavalier::get_var_info(sv=TRUE),
                               cohort_AC_AF = 'cohort_AC_AF'))
  } else {
    file.create(str_c(opts$out, '.pptx'))
  }
  
  # write candidate info
  cand_vars %>%
    mutate(set = 'SV',
           family = opts$family) %>%
    select(set, family, symbol, consequence, inheritance, id, SVTYPE, chrom, pos, END, SVLEN) %>%
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