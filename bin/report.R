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
  report.R <snv_vars> <ped> <sample_bams> <gene_lists> <config> ... [options]

Options:
  snv_vars                    Input Variant TSV file.
  ped                         Pedigree file.
  sample_bams                 Sample names and bam files in format name=/path/to/bam.
  gene_lists                  Comma serparated list of gene list filenames.
  config                      Config Json file with additional params.
  --out=<prefix>              Output file prefix [default: out].
  --cav-opts=<json>           Json file with additional options for cavalier package.
  --func-source=<rscript>     Custom functions source files, comma separated. 
  --no-slides                 Don't create pptx output.                 
"
opts <- docopt(doc)

message('Using CLI options:')
opts[names(opts) %>%
       keep(~str_detect(., '^[:alpha:]'))] %>%
  { class(.) <- c('list', 'docopt'); .} %>%
  print()

if (!is.null(opts$cav_opts)) {
  set_options_from_json(opts$cav_opts)
}

################# CHECK ARGS ######################
sample_bams <-
  c(str_split(opts$sample_bams, pattern = ',', simplify = T)) %>% 
  (function(x) setNames(str_extract(x, '[^=]+$'), str_extract(x, '^[^=]+')))

src_files  <- c(str_split(opts$func_source, pattern = ',', simplify = T))
gene_lists <- c(str_split(opts$gene_lists, ',', simplify = T))

stopifnot(
  file.exists(opts$snv_vars),
  file.exists(opts$ped),
  all(file.exists(sample_bams)),
  all(file.exists(src_files)),
  all(file.exists(gene_lists))
)

################ CHECK CONF ###################
conf <- jsonlite::read_json(opts$conf, simplifyVector=TRUE)
conf$snv$fields <- unlist(conf$snv$fields)
message('Parsed config as:\n', conf)

stopifnot(
  is.list(conf$snv),
  rlang::is_scalar_character(conf$snv$functions),
  rlang::is_named(conf$snv$fields),
  is.character(conf$snv$fields)
)

snv_functions <- c(str_split(conf$snv$functions, pattern = ',', simplify = T))

######### HELPER FUNCTIONS ###################
update_stats <- function(curr, prev, reason, removed_df = tibble(), keys = c('CHROM', 'POS', 'REF', 'ALT', 'Gene', 'SYMBOL')) {
  distinct(
    bind_rows(
      removed_df,
      select(prev, all_of(keys)) %>% 
        anti_join(curr, by = keys) %>% 
        mutate(REASON = reason)
    )
  )
}

############## READ PEDIGREE   ################
PEDIGREE   <- read_ped(opts$ped)

############## READ GENE LISTS ###############
GENE_LIST <-
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
  nest(panel_data = -ensembl_gene_id) %>% 
  mutate(hgnc_symbol = hgnc_ensembl2sym(ensembl_gene_id), .after = ensembl_gene_id)


######## START SUBPROCESS AND CHECK FUNCS #############
script_session <- callr::r_session$new()

session_funcs <-
  script_session$run(
    function(SRC_FILES, SCRIPT_FUNCS) {
      for (SRC in SRC_FILES){
        source(SRC)
      }
      lapply(
        setNames(SCRIPT_FUNCS, SCRIPT_FUNCS),
        function(x) tryCatch(get(x), error = function(e) NULL)
      )
    },
    args = list(src_files, snv_functions)
  ) %>% 
  keep(is.function) %>% 
  names()

if (! all(snv_functions %in% session_funcs)) {
  stop(
    "missing functions: ", 
    str_c(setdiff(snv_functions, session_funcs), collapse = ', '),
    '\ncheck script: ',  str_c(src_files, collapse = ', '))
}

############# APPLY SNV FUNCTIONS #################
ARGS <-
  list(
    FILEPATH  = opts$snv_vars,
    PEDIGREE  = PEDIGREE,
    GENE_LIST = GENE_LIST,
    CONFIG = conf$snv,
    VARIANTS = NULL
  )

snv_stats <- tibble()

for (func in snv_functions) {
  message("Running function: ", func)
  VARIANTS <-
    script_session$run(
      function(FUNCTION, ARGS) {
        f <- get(FUNCTION)
        do.call(f, ARGS)
      },
      args = list(func, ARGS)
    )
  message(
    '  - Retured ', format(nrow(VARIANTS), big.mark = ','), ' SNV variants with ', ncol(VARIANTS), ' columns'
  )
  if (!is.null(ARGS$VARIANTS)) {
    snv_stats <- update_stats(VARIANTS, ARGS$VARIANTS, func, snv_stats)
  }
  ARGS$VARIANTS <- VARIANTS
}

############### Update stats ########################
snv_stats <-
  update_stats(filter(VARIANTS, FALSE), VARIANTS, 'PASS', snv_stats) %>% 
  left_join(
     select(VARIANTS, CHROM, POS, REF, ALT, Gene, SYMBOL, TYPE),
     by = c('CHROM', 'POS', 'REF', 'ALT', 'Gene', 'SYMBOL')
  )

snv_variants <- VARIANTS

############### CREATE SLIDES ###############
# slide_layout
layout <-
  `if`(
    length(sample_bams) == 1,
    slide_layout(
      c('var_info', 'igv', 'gtex'),
      c('omim', 'custom_panel_data'),
      heights = c(23, 8),
      title_height = 0.09,
      pad = 0.015
    ),
    bind_rows(
      slide_layout(
        c('var_info', 'pedigree', 'gtex'),
        c('omim', 'custom_panel_data'),
        heights = c(23, 8),
        title_height = 0.09,
        pad = 0.015
      ),
      slide_layout(
        c('igv'),
        slide_num = 2L,
        title_height = 0.09,
        pad = 0.015
      )
    )
  )

slide_data <-
  snv_variants %>% 
  select(
    title,
    CHROM, POS, REF, ALT, SYMBOL,
    genotype = GT,
    # this will rename these as required
    all_of(conf$snv$fields),
    any_of(
      setNames(
        str_c(conf$snv$fields, '_url'),
        str_c(names(conf$snv$fields), '_url')
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
    GENE_LIST,
    by = 'ensembl_gene_id'
  )

if ((nrow(slide_data) > 0) & !opts$no_slides) {
  
  message('Creating slides for ', nrow(slide_data), ' variants')
  create_slides(
    variants = slide_data,
    output = str_c(opts$out, '.snv.pptx'),
    bam_files = sample_bams,
    ped_file = opts$ped,
    layout = layout,
    var_info = c(names(conf$snv$fields), keep(colnames(slide_data), str_ends, '_url')),
    # var_info = c(names(conf$snv$fields)),
  )
} else {
  if (opts$no_slides) {
    message('Slides disabled, creating placeholder file')
    file.create(str_c(opts$out, '.snv.pptx'))
  } else {
    message('No variants remain, creating placeholder slides')
    create_slides(variants = tibble(), layout = tibble(), output = str_c(opts$out, '.snv.pptx'))
  }
}

################ CREATE OUTPUTS ##################
snv_variants %>% 
  select(where(~ !is.data.frame(.))) %>% 
  write_csv(str_c(opts$out, '.snv_candidates.csv'))

snv_stats %>% 
  count(REASON, TYPE) %>% 
  write_csv(str_c(opts$out, '.snv_filter_stats.csv'))

snv_stats %>% 
  write_csv(str_c(opts$out, '.snv_filter_reason.csv.gz'))

write_tsv(
  transmute(snv_variants,
            CHROM = CHROM, 
            START = POS - 1L,
            END = START + if_else(nchar(ALT) > nchar(REF), 2L, nchar(ALT)),
            title = str_c(SYMBOL, ' - ', broad_id)),
  str_c(opts$out, '.igv.bed.gz'),
  col_names = F
)



