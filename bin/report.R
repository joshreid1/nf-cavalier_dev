#!/usr/bin/env Rscript

stopifnot(
  require(tidyverse),
  require(magrittr),
  require(docopt),
  require(assertthat),
  require(cavalier)
)

doc <- "
Usage:
  report.R <snv_vars> <ped> <sample_bams> <gene_lists> <conf> ... [options]

Options:
  snv_vars                    Input Variant TSV file.
  ped                         Pedigree file.
  sample_bams                 Sample names and bam files in format name=/path/to/bam.
  gene_lists                  Comma serparated list of gene list filenames.
  config                      Config Json file with additional params.
  --out=<f>                   Output file prefix [default: out].
  --family=<f>                Name of sample/family [default: Family].
  --cav-opts=<f>              Json file with additional options for cavalier package.
  --no-slides                 Don't create pptx output.
"
opts <- docopt(doc)

message('Using CLI options:')
opts[names(opts) %>%
       keep(~str_detect(., '^[:alpha:]'))] %>%
  { class(.) <- c('list', 'docopt'); .} %>%
  print()
  
# setwd('/vast/scratch/users/munro.j/nextflow/work/8f/114b239dfab0375bf1175e75af30fa')
# opts <- list(
#   snv_vars = 'test_snp.chr6_chr19.clean.vcfanno.flt.vep.family_Plu_PK86442.tsv.gz',
#   ped  = 'Plu_PK86442.ped',
#   sample_bams = 'T23547=T23547.bam,T23548=T23548.bam',
#   conf = 'report_conf.json',
#   gene_lists = 'G4E_ALL_v2025-03.tsv,list_entrez.mapped.tsv',
#   family = 'Plu_PK8644',
#   out = 'out',
#   cav_opts = 'cavalier_opts.json'
# )

if (!is.null(opts$cav_opts)) {
  set_options_from_json(opts$cav_opts)
}

sample_bams <-
  c(str_split(opts$sample_bams, pattern = ',', simplify = T)) %>% 
  (function(x) setNames(str_extract(x, '[^=]+$'), str_extract(x, '^[^=]+')))

invisible(
  stopifnot(
    file.exists(opts$snv_vars),
    file.exists(opts$ped),
    all(file.exists(sample_bams))
  )
)

########### READ CONF ###############
conf <- jsonlite::read_json(opts$conf)

conf$snv_freq_filters <- 
  map(conf$snv_freq_filters, function(x) {
    map(x, function(y) {
      map(y, function(z) {
        `if`(is.numeric(z), z, Inf)
      })
    })
  })

conf$snv_report <- unlist(conf$snv_report)

conf$snv_subsets <-
  map(conf$snv_subsets, function(x) {
    x$report <- unlist(x$report)
    if(is.null(x$report)) {
      x$report <- character()
    }
    x
  })

############## READ GENE LISTS ###############
list_df <-
  c(str_split(opts$gene_lists, ',', simplify = T)) %>% 
  map_df(function(fn) {
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


update_removed <- function(curr, prev, reason,
                           removed_df = tibble(),
                           keys = c('CHROM', 'POS', 'REF', 'ALT', 'Gene')
                           ) {
  bind_rows(
    removed_df,
    select(prev, all_of(keys)) %>% 
      anti_join(curr, by = keys) %>% 
      mutate(REASON = reason)
  ) %>% distinct()
}

############## READ SNV VARIANTS ###############
snv_vars <- 
  read_tsv(
    opts$snv_vars,
    na = '.',
    col_types = cols(.default = col_character())
  ) %>% 
  readr::type_convert(guess_integer = TRUE) %>% 
  mutate(
    genotype = 
      pick(starts_with('FMT_GT_')) %>% 
      rename_with(~str_remove(., 'FMT_GT_')),
    GENOTYPE = 
      genotype %>% 
      mutate_all(function(x) {
        case_when(
          # TODO - handle hemizygous variants???
          str_count(x, '[01]') == 2 & str_count(x, '[1]') == 2 ~ 'HOM',
          str_count(x, '[01]') == 2 & str_count(x, '[1]') == 1 ~ 'HET',
          str_count(x, '[01]') == 2 & str_count(x, '[1]') == 0 ~ 'REF'
        ) %>% factor(c('REF', 'HET', 'HOM'))
      }),
    DEPTH = 
      pick(starts_with('FMT_AD_')) %>% 
      rename_with(~str_remove(., 'FMT_AD_')) %>% 
      mutate_all(function(x) map_int( str_split(x, ","), ~ sum(as.integer(.x)) )),
  ) %>%
  # select(!starts_with(c('FMT_GT_', 'FMT_AD_'))) %>% 
  mutate(
    # remove family counts from cohort counts for filtering
    AN = pmax(0, AN - 2*rowSums(!is.na(GENOTYPE))),
    AC = pmax(0, AC - rowSums(GENOTYPE == 'HET', na.rm = TRUE) - 2*rowSums(GENOTYPE == 'HOM', na.rm = TRUE)),
    AF = if_else(AN > 0, AC / AN, 0),
  ) %>% 
  distinct() %>% 
  mutate(variant_id = str_c(HGVSg, ':', Gene)) %>% 
  mutate(family_id = opts$family, .before = 1)

###### FILTER BASED ON DEPTH ######################
snv_vars_flt <-
  filter(
    snv_vars,
    rowSums(DEPTH >= conf$snv_min_depth) == ncol(DEPTH)
  )
snv_removed <- update_removed(snv_vars_flt, snv_vars, '1_LOW_DEPTH')
snv_vars <- snv_vars_flt

######## INTERSECT WITH GENE LIST ##################
snv_vars_flt <-
  filter(
    snv_vars, 
    Gene %in% list_df$ensembl_gene_id
  )
snv_removed <- update_removed(snv_vars_flt, snv_vars, '2_GENE_LIST', snv_removed)
snv_vars <- snv_vars_flt

############### CUSTOM DATA MUTATIONS ############### 
if (length(names(conf$snv_mutate)) > 0) {
  for (col in names(conf$snv_mutate)) {
    expr <- rlang::parse_expr(conf$snv_mutate[[col]])
    message('mutating snvs with: ', col, ' = ', rlang::expr_text(expr))
    snv_vars <- mutate(snv_vars, !! col := !! expr)
  }
}

###### INHERITANCE AND FREQUENCY FILTERING ##########
# requires pop_AF/pop_AC/AF/AC
for(col in c('AC', 'AF', 'pop_AC', 'pop_AF')) {
  if (!col %in% colnames(snv_vars)) {
    warning('snv ', col, ' is not defined - no filtering on ', col, ' will be performed (tip: add with mutate)' )
    snv_vars[[col]] <- 0
  }
}
ped_df   <- read_ped(opts$ped)
aff  <- intersect(get_affected(ped_df)  , colnames(snv_vars$GENOTYPE))
una  <- intersect(get_unaffected(ped_df), colnames(snv_vars$GENOTYPE))

snv_vars_flt <-
  snv_vars %>% 
  mutate(inheritance = case_when(
    (rowSums(snv_vars$GENOTYPE[, aff, drop = F] == 'HET') == length(aff) &
      rowSums(snv_vars$GENOTYPE[, una, drop = F] == 'REF') == length(una)) ~ 'dominant',
    (rowSums(snv_vars$GENOTYPE[, aff, drop = F] == 'HOM') == length(aff) &
       rowSums(snv_vars$GENOTYPE[, una, drop = F] != 'HOM') == length(una)) ~ 'recessive',
  )) %>% 
  filter(!is.na(inheritance)) %>% 
  (function(x) bind_rows(
    x,
    filter(x, inheritance == 'dominant') %>% mutate(inheritance = 'compound')
  )) %>% 
  filter(
    or(# pop_AF
      inheritance %in% c('recessive', 'compound') & pop_AF < conf$snv_freq_filters$pop$AF$recessive,
      inheritance == 'dominant'                   & pop_AF < conf$snv_freq_filters$pop$AF$dominant
    ),
    or(# pop_AC
      inheritance %in% c('recessive', 'compound') & pop_AC < conf$snv_freq_filters$pop$AC$recessive,
      inheritance == 'dominant'                   & pop_AC < conf$snv_freq_filters$pop$AC$dominant
    ),
    or(# AF
      inheritance %in% c('recessive', 'compound') & AF     < conf$snv_freq_filters$cohort$AF$recessive,
      inheritance == 'dominant'                   & AF     < conf$snv_freq_filters$cohort$AF$dominant
    ),
    or(# AC
      inheritance %in% c('recessive', 'compound') & AC     < conf$snv_freq_filters$cohort$AC$recessive,
      inheritance == 'dominant'                   & AC     < conf$snv_freq_filters$cohort$AC$dominant
    )
  ) %>% 
  chop(inheritance) %>% 
  mutate(
    inheritance = map(inheritance, unique) %>% map_chr(str_c, collapse = '&'),
    inheritance = if_else(str_detect(inheritance, 'dominant'), 'dominant', inheritance)) 

snv_removed <- update_removed(snv_vars_flt, snv_vars, '3_FREQ_AND_INHERITANCE', snv_removed)
snv_vars <- snv_vars_flt

############### SUBSET FILTERING ############### 
snv_vars$subset <- NA_character_
for (col in names(conf$snv_subsets)) {
  # col <- 'missense'
  sub <- filter(snv_vars, is.na(subset))
  for (ex in conf$snv_subsets[[col]]$filter) {
    expr <- rlang::parse_expr(ex)
    message('filtering for "', col, '" snvs with:\n  ', rlang::expr_text(expr))
    sub <- filter(sub, !! expr)
  }
  snv_vars <- 
    mutate(snv_vars,
    subset = if_else(variant_id %in% sub$variant_id, col, subset)
  )
}
# snv_vars<- snv_vars %>% mutate(subset = if_else(SpliceAI_max > 0.1, 'splicing', subset))
snv_vars_flt <- filter(snv_vars, !is.na(subset))
snv_removed <- update_removed(snv_vars_flt, snv_vars, '4_SUBSET_FILTERS', snv_removed)
snv_vars <- snv_vars_flt

############### COMPOUND HET FILTER ############### 
snv_vars_flt <-
  snv_vars %>% 
  add_count(Gene, dom_comp = inheritance %in% c('dominant', 'compound'), name = 'n_compound') %>%
  filter(inheritance != 'compound' | n_compound > 1) %>% 
  select(-n_compound, -dom_comp)
snv_removed <- update_removed(snv_vars_flt, snv_vars, '5_COMPOUND_HET', snv_removed)
snv_vars <- snv_vars_flt

############### PASS ########################
snv_removed <- update_removed(filter(snv_vars, FALSE), snv_vars, 'PASS', snv_removed)

############### CREATE SLIDES ###############

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

slide_data <-
  snv_vars %>% 
  select(
    subset,
    title,
    CHROM, POS, REF, ALT, SYMBOL,
    genotype,
    all_of(
      union(
        conf$snv_report,
        map(conf$snv_subsets, 'report') %>% unname() %>% unlist()
      )
    )
  ) %>% 
  mutate(
    # add keys required for create_slides
    chrom = CHROM,
    pos = POS,
    ref = REF,
    alt = ALT,
    ensembl_gene_id = Gene,
    symbol = coalesce(
      hgnc_entrez2sym(hgnc_ensembl2entrez(Gene)),
      SYMBOL
    )
  ) %>% 
  left_join(
    list_df,
    by = c('Gene' = 'ensembl_gene_id')
  ) %>% 
  split.data.frame(.$subset)

if (!opts$no_slides) {
  slides <- cavalier::get_slide_template()
  
  for (subset in names(slide_data)) {
    message('Creating slides for "', subset, '" snvs')
    new_slides <- tempfile(str_c('.', subset, '.'), tmpdir = getwd(), fileext = '.pptx')
    create_slides(
      variants = slide_data[[subset]],
      slide_template = slides,
      output = new_slides,
      bam_files = sample_bams,
      ped_file = opts$ped,
      layout = layout, #[1:5,]
      var_info = c(conf$snv_report, conf$snv_subsets[[subset]]$report) %>% (\(x) x[!duplicated(x)]),
    )
    slides <- new_slides
  }
}

if (length(slide_data) == 0 | opts$no_slides) {
  # no variants/ no slides
  create_slides(variants = tibble(), layout = tibble(), output = str_c(opts$out, '.snv.pptx'))
} else {
  file.rename(slides, str_c(opts$out, '.snv.pptx'))
}

snv_vars %>% 
  select(where(~ !is.data.frame(.))) %>% 
  write_csv(str_c(opts$out, '.snv_candidates.csv'))

snv_removed %>% 
  count(REASON) %>% 
  write_csv(str_c(opts$out, '.snv_filter_stats.csv'))

snv_removed %>% 
  write_csv(str_c(opts$out, '.snv_filter_reason.csv.gz'))

############ FOR IGV  ##################
write_tsv(
  transmute(snv_vars,
            CHROM = CHROM, 
            START = POS - 1L,
            END = START + if_else(nchar(ALT) > nchar(REF), 2L, nchar(ALT)),
            title = str_c(SYMBOL, ' ', HGVSg)),
  str_c(opts$out, '.igv.bed.gz'),
  col_names = F
)



