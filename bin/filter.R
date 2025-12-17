#!/usr/bin/env Rscript

stopifnot(
  require(tidyverse),
  require(magrittr),
  require(docopt),
  require(cavalier)
)

# # NB these are set by Nextflow and should be left commented out except for testing
# options(
#   # FMT
#   FILTER_SHORT_MIN_DP = 5,
#   FILTER_SHORT_MIN_GQ = 10,
#   # POPULATION SHORT AF/AC/HOM
#   FILTER_SHORT_POP_DOM_MAX_AF = 0.0001,
#   FILTER_SHORT_POP_REC_MAX_AF = 0.01,
#   FILTER_SHORT_POP_DOM_MAX_AC = 20,
#   FILTER_SHORT_POP_REC_MAX_AC = 100,
#   FILTER_SHORT_POP_DOM_MAX_HOM = 5,
#   FILTER_SHORT_POP_REC_MAX_HOM = 10,
#   # COHORT SHORT AF/AC
#   FILTER_SHORT_COH_DOM_MAX_AF = NULL,
#   FILTER_SHORT_COH_REC_MAX_AF = NULL,
#   FILTER_SHORT_COH_DOM_MAX_AC = NULL,
#   FILTER_SHORT_COH_REC_MAX_AC = NULL,
#   # POPULATION SV AF/AC/HOM
#   FILTER_STRUC_POP_DOM_MAX_AF = 0.0001,
#   FILTER_STRUC_POP_REC_MAX_AF = 0.01,
#   FILTER_STRUC_POP_DOM_MAX_AC = 20,
#   FILTER_STRUC_POP_REC_MAX_AC = 100,
#   FILTER_STRUC_POP_DOM_MAX_HOM = NULL,
#   FILTER_STRUC_POP_REC_MAX_HOM = NULL,
#   # COHORT SV AF/AC
#   FILTER_STRUC_COH_DOM_MAX_AF = 0.01,
#   FILTER_STRUC_COH_REC_MAX_AF = 0.01,
#   FILTER_STRUC_COH_DOM_MAX_AC = NULL,
#   FILTER_STRUC_COH_REC_MAX_AC = NULL,
#   # CLINVAR
#   FILTER_SHORT_CLINVAR_KEEP_PAT  = '(p|P)athogenic(?!ity)',
#   FILTER_SHORT_CLINVAR_DISC_PAT  = '(b|B)enign',
#   FILTER_SHORT_CLINVAR_ANYWHERE = TRUE,
#   # Clingen PP3 supporting  doi: 10.1016/j.ajhg.2022.10.013
#   FILTER_SHORT_LOF              = TRUE,
#   FILTER_SHORT_MISSENSE         = TRUE,
#   FILTER_SHORT_SPLICING         = TRUE,
#   FILTER_SHORT_MIN_CADD_PP      = 25.3,
#   # Clingen PP3 supporting https://doi.org/10.1016/j.ajhg.2023.06.002
#   FILTER_SHORT_MIN_SPLICEAI_PP  = 0.20,
#   FILTER_SHORT_VEP_MIN_IMPACT   = 'MODERATE',
#   FILTER_SHORT_VEP_CONSEQUENCES = NULL
# )

MAIN <- function(opts) {
  
  message('Using CLI options:')
  PRINT_STR(
    opts[str_detect(names(opts), '^[:alpha:]')] 
  )

  ############### SET OPTIONS ######################
  do.call(options, jsonlite::fromJSON(opts$options))
  message('Using filter options:')
  PRINT_STR(
    options() %>% {.[str_starts(names(.), 'FILTER_')]}
  )
  if (!is.null(opts$cav_opts)) {
    cavalier::set_options_from_json(opts$cav_opts)
  }
  
  gene_lists <- c(str_split(opts$gene_lists, ',', simplify = T))
  
  ################# CHECK ARGS ######################
  stopifnot(
    file.exists(opts$short_var),
    file.exists(opts$ped),
    all(file.exists(gene_lists))
  )
  
  ############## READ PEDIGREE   ################
  
  pedigree <- cavalier::read_ped(opts$ped)
  
  ##############  GENE LIST ###############
  gene_set <- unique(na.omit(read_lines(opts$gene_set)))
  message('Searching for variants in ', length(gene_set), ' genes')
  
  ########## SHORT VARS #################
  short_cand <- LOAD_AND_FILTER_SHORT(
    short_var = opts$short_var, 
    gene_set = gene_set,
    pedigree = pedigree
  )
  
  ############## TODO: STRUC VARS #################
  
  
  ########## COMPOUND FILTER(BOTH) ################
  
  FC <- FILTER_COMPOUND(SHORT_VAR = short_cand)
  
  ######### OUTPUTS ###############################

  if (!is.null(FC$SHORT)) {
    
    saveRDS(FC$SHORT, str_c(opts$output, '.short.filtered_variants.rds'))
    
    FC$SHORT %>% 
      select(where(~ !is.data.frame(.))) %>% 
      write_csv(
        file = str_c(opts$output, '.short.filtered_variants.csv')
      )
    
    write_csv(
      .GlobalEnv$.tracking.SHORT$filtered,
      file = str_c(opts$output, '.short.reason_filtered.csv.gz')
    )
    
    FC$SHORT %>% 
      transmute(
        CHROM = CHROM, 
        START = POS - 1L,
        END = START + if_else(nchar(ALT) > nchar(REF), 2L, nchar(ALT)),
        title = variant_id,
      ) %>% 
      distinct() %>% 
      write_tsv(
        str_c(opts$output, '.short.igv.bed.gz'),
        col_names = F
      )

    cat(nrow(FC$SHORT), file = str_c(opts$output, '.short.count'))
    
  } else {
    file.create(str_c(opts$output, '.empty.short.filtered_variants.rds'))
    file.create(str_c(opts$output, '.empty.short.filtered_variants.csv'))
    file.create(str_c(opts$output, '.empty.short.reason_filtered.csv.gz'))
    file.create(str_c(opts$output, '.empty.short.igv.bed.gz'))
    cat("0", file = str_c(opts$output, '.short.count'))
  }
  
  if (!is.null(FC$STRUC)) {
    
    saveRDS(FC$STRUC, str_c(opts$output, '.struc.filtered_variants.rds'))
    
    FC$STRUC %>% 
      select(where(~ !is.data.frame(.))) %>% 
      write_csv(
        file = str_c(opts$output, '.struc.filtered_variants.csv')
      )
    
    write_csv(
      .GlobalEnv$.tracking.STRUC$filtered,
      file = str_c(opts$output, '.struc.reason_filtered.csv.gz')
    )

    cat(nrow(FC$STRUC), file = str_c(opts$output, '.struc.count'))

  } else {
    file.create(str_c(opts$output, '.empty.struc.filtered_variants.rds'))
    file.create(str_c(opts$output, '.empty.struc.filtered_variants.csv'))
    file.create(str_c(opts$output, '.empty.struc.reason_filtered.csv.gz'))
    cat("0", file = str_c(opts$output, '.struc.count'))
  }
  
  message('FILTERING COMPLETE')
}

############ LOAD_AND_FILTER_SHORT  ###############################
LOAD_AND_FILTER_SHORT <- function(short_var, gene_set, pedigree) {
  
  short_cand <- 
    LOAD_SHORT(short_var) %>% 
    FILTER_SHORT_TYPE() %>% 
    FILTER_GENE(gene_set, set = "SHORT") %>%
    FILTER_INHERITANCE(pedigree, set = "SHORT") %>%
    FILTER_SHORT_FMT()
}

###################### LOAD_SHORT ###############################
LOAD_SHORT <- function(FILEPATH = NULL, ...) {
  # check libraries and args
  stopifnot(
    require(tidyverse),
    rlang::is_scalar_character(FILEPATH),
    file.exists(FILEPATH)
  )
  
  # expected types for columns from default annotation sources
  col_spec <- cols(
    CHROM = col_character(),
    POS = col_integer(),
    REF = col_character(),
    ALT = col_character(),
    AC = col_integer(),
    AF = col_double(),
    AN = col_integer(),
    CADD = col_double(),
    gnomad_AC = col_integer(),
    gnomad_AF = col_double(),
    gnomad_FAF95 = col_double(),
    gnomad_nhomalt = col_integer(),
    phyloP100 = col_double(),
    SYMBOL = col_character(),
    Gene = col_character(),
    VARIANT_CLASS = col_character(),
    Consequence = col_character(),
    IMPACT = col_character(),
    Feature_type = col_character(),
    Feature = col_character(),
    BIOTYPE = col_character(),
    EXON = col_character(),
    INTRON = col_character(),
    HGVSc = col_character(),
    HGVSp = col_character(),
    HGVSg = col_character(),
    Amino_acids = col_character(),
    HGNC_ID = col_character(),
    MANE = col_logical(),
    MANE_SELECT = col_logical(),
    MANE_PLUS_CLINICAL = col_logical(),
    CCDS = col_logical(),
    ENSP = col_character(),
    SIFT = col_character(),
    PolyPhen = col_character(),
    CLNSIG = col_character(),
    CLNGENE = col_character(),
    CLNVID = col_integer(),
    Existing_variation = col_character(),
    am_class = col_character(),
    am_pathogenicity = col_double(),
    REVEL = col_double(),
    SpliceAI_pred_DP_AG = col_integer(),
    SpliceAI_pred_DP_AL = col_integer(),
    SpliceAI_pred_DP_DG = col_integer(),
    SpliceAI_pred_DP_DL = col_integer(),
    SpliceAI_pred_DS_AG = col_double(),
    SpliceAI_pred_DS_AL = col_double(),
    SpliceAI_pred_DS_DG = col_double(),
    SpliceAI_pred_DS_DL = col_double(),
    SpliceAI_pred_SYMBOL = col_character(),
    `_5UTR_annotation` = col_character(),
    `_5UTR_consequence` = col_character(),
    Existing_InFrame_oORFs = col_integer(),
    Existing_OutOfFrame_oORFs = col_integer(),
    Existing_uORFs = col_integer(),
    .default = col_character()
  )
  
  VARIANTS_OUT <-
    # read tsv with custom col_spec
    read_tsv(
      file = FILEPATH,
      na = '.',
      col_types = col_spec
    ) %>% 
    # try to guess column types not in col_spec, but better to add to col_spec
    readr::type_convert(
      guess_integer = TRUE
    ) %>% 
    # convert FMT fields to nested data frames
    FMT_TO_DF() %>%
    # add/modify columns
    mutate(
      CLNSIG = if_else(CLNSIG_GENE_MATCH(CLNGENE, SYMBOL, Gene), CLNSIG, NA_character_),
      AN = AN - (GT %>% mutate_all(~ str_count(., '[01]')) %>% rowSums(na.rm = TRUE)),
      AC = AC - (GT %>% mutate_all(~ str_count(., '[1]')) %>% rowSums(na.rm = TRUE)),
      AF = replace_na(AC / AN, 0),
      # useful for filtering
      pop_AF = pmax(replace_na(gnomad_AF, 0), replace_na(gnomad_FAF95, 0)),
      pop_AC =  replace_na(gnomad_AC, 0),
      pop_hom = replace_na(gnomad_nhomalt, 0),
      SpliceAI_max = 
        replace_na(
          pmax(SpliceAI_pred_DS_AG, SpliceAI_pred_DS_AL, SpliceAI_pred_DS_DG, SpliceAI_pred_DS_DL),
          0
        ),
      SIFT_score     = str_extract(SIFT    , '(?<=\\()[0-9\\.]+(?=\\)$)') %>% as.numeric(),
      PolyPhen_score = str_extract(PolyPhen, '(?<=\\()[0-9\\.]+(?=\\)$)') %>% as.numeric(),
      # not meaningful for insertions
      phyloP100 = replace(phyloP100, VARIANT_CLASS == "insertion", NA),
      variant_id = str_c(CHROM, POS, REF, ALT, sep = '-')
    )
  
  return(
    VARIANTS_OUT %>% TRACK_REASON(init=TRUE)
  )
}

####################### FMT_TO_DF #########################
# helper to convert sample format fields to a dataframe with colnames
FMT_TO_DF <- function(input_data) {
  stopifnot(
    is.data.frame(input_data)
  )
  
  format_fields <- 
    colnames(input_data) %>% 
    keep(str_detect, '^FMT_[A-Za-z]+(?=_)') %>% 
    str_extract('(?<=^FMT_)[A-Za-z]+(?=_)') %>% 
    unique()
  
  output <- input_data
  
  for (FMT in format_fields) {
    prefix <- str_c('FMT_', FMT, '_')
    output <-
      output %>% 
      mutate(
        !!FMT :=
          pick(starts_with(prefix)) %>% 
          rename_with(~ str_remove(., prefix)) %>% 
          {suppressWarnings(
            readr::type_convert(., guess_integer = T, col_types = cols())
          )} 
      )
  }
  
  return(output)
}

#################### CLNSIG_GENE_MATCH #################### 
# helper to restrict ClinVar CLNSIG annotation to matching genes
CLNSIG_GENE_MATCH <- function(CLNGENE, SYMBOL, Gene) {
  # match by symbol
  match1 <- map2_lgl(
    SYMBOL,
    str_extract_all(CLNGENE, '[^|:]+(?=:)'),
    ~ .x %in% .y
  )
  # match by entrez id
  match2 <- map2_lgl(
    as.character(cavalier::hgnc_ensembl2entrez(Gene)), 
    str_extract_all(CLNGENE, '(?<=:)[^|:]+'),
    ~ .x %in% .y
  )
  replace_na(match1 | match2, FALSE)
}

#################### TRACK_REASON #################### 
# keep track of reason variants are removed for debugging
TRACK_REASON <- function(data, reason = NA_character_, init=FALSE, set = 'SHORT') {
  
  force(data)
  
  if (set == 'SHORT') {
    keys = c('CHROM', 'POS', 'REF', 'ALT', 'Gene', 'SYMBOL')
  } else {
    keys = c('CHROM', 'POS', 'REF', 'ALT', 'SVTYPE', 'SVLEN', 'END')
  }
  env_name <- str_c('.tracking.', set)
  
  if (init) {
    .GlobalEnv[[env_name]] <- new.env()
    env <- .GlobalEnv[[env_name]]
    env$filtered <- tibble()
    env$prev <- data
    env$step <- 0L
    message('Initial short Gene-Variants: n=', nrow(data))
  } else{
    env <- .GlobalEnv[[env_name]]
    env$step <- env$step + 1L
    
    new_flt <- 
      select(env$prev, all_of(keys)) %>% 
      anti_join(data, by = keys) %>% 
      mutate(REASON = reason, STEP = env$step) %>% 
      distinct()
    
    env$filtered <-
      bind_rows(
        env$filtered,
        new_flt
      )
    message(reason, ' removed ', nrow(new_flt), ' short Gene-Variants, ', nrow(data), ' remain')
    env$prev <- data
  }
  
  return(data)
}

###################### FILTER_SHORT_FMT ###############################
# filter shorts by DP and GQ format fields
FILTER_SHORT_FMT <- function(VARIANTS) {
  stopifnot(
    require(tidyverse),
    is.data.frame(VARIANTS)
  )
  
  VARIANTS_OUT <- VARIANTS 
  
  MIN_DP <- getOption('FILTER_SHORT_MIN_DP')
  if (!is.null(MIN_DP)) {
    stopifnot( 'DP' %in% colnames(VARIANTS) )
    VARIANTS_OUT <-
      VARIANTS_OUT %>% 
      filter(do.call(pmin, unname(DP)) >= MIN_DP)
  }
  
  MIN_GQ <- getOption('FILTER_SHORT_MIN_GQ')
  if (!is.null(MIN_GQ)) {
    stopifnot( 'GQ' %in% colnames(VARIANTS) )
    VARIANTS_OUT <-
      VARIANTS_OUT %>% 
      filter(do.call(pmin, unname(GQ)) >= MIN_GQ)
  }
  
  return(
    VARIANTS_OUT %>% TRACK_REASON('FILTER_SHORT_FMT')
  )
  
}

###################### FILTER_SHORT_GENE ###############################
# filter for pressence in certain genes or clinvar pathogenic
FILTER_GENE <- function(VARIANTS, GENE_SET, set = 'SHORT') {
  stopifnot(
    require(tidyverse),
    is.data.frame(VARIANTS),
    is.character(GENE_SET)
  )

  FILTER_SHORT_CLINVAR_KEEP_PAT <- '$.' # never TRUE
  FILTER_SHORT_CLINVAR_ANYWHERE <- getOption('FILTER_SHORT_CLINVAR_ANYWHERE', FALSE)
  if (FILTER_SHORT_CLINVAR_ANYWHERE) {
    FILTER_SHORT_CLINVAR_KEEP_PAT <- getOption('FILTER_SHORT_CLINVAR_KEEP_PAT', '$.')
  }

  VARIANTS_OUT <-
    VARIANTS %>%
    filter(
      Gene %in% GENE_SET |
        str_detect(CLNSIG, FILTER_SHORT_CLINVAR_KEEP_PAT)
    )
  
  return(
    VARIANTS_OUT %>% TRACK_REASON(str_c('FILTER_GENE(', set, ')'))
  )
}

###################### FILTER_SHORT_TYPE ###############################
# Filter short variants by type of interest (i.e. VEP consequence)
FILTER_SHORT_TYPE <- function(VARIANTS) {
  stopifnot(
    require(tidyverse),
    is.data.frame(VARIANTS)
  )
  
  FILTER_SHORT_LOF      <- getOption('FILTER_SHORT_LOF'     , TRUE)
  FILTER_SHORT_MISSENSE <- getOption('FILTER_SHORT_MISSENSE', TRUE)
  FILTER_SHORT_SPLICING <- getOption('FILTER_SHORT_SPLICING', TRUE)
  
  FILTER_SHORT_CLINVAR_KEEP_PAT <- getOption('FILTER_SHORT_CLINVAR_KEEP_PAT', '$.')
  FILTER_SHORT_CLINVAR_DISC_PAT <- getOption('FILTER_SHORT_CLINVAR_DISC_PAT', '$.')
  
  FILTER_SHORT_VEP_MIN_IMPACT <- 
    getOption('FILTER_SHORT_VEP_MIN_IMPACT', 'MODERATE')
  
  FILTER_SHORT_VEP_IMPACTS <-
    list(
      MODIFIER = c('MODIFIER', 'LOW', 'MODERATE', 'HIGH'),
      LOW = c('LOW', 'MODERATE', 'HIGH'),
      MODERATE = c('MODERATE', 'HIGH'),
      HIGH = c('HIGH')
    )[[FILTER_SHORT_VEP_MIN_IMPACT]]
  
  FILTER_SHORT_VEP_CONSEQUENCES <- 
    getOption('FILTER_SHORT_VEP_CONSEQUENCES', '') %>% 
    str_split(',', simplify = T) %>% 
    c()
  
  VARIANTS_OUT <-
    VARIANTS %>% 
    filter(!replace_na(str_detect(CLNSIG, FILTER_SHORT_CLINVAR_DISC_PAT), FALSE)) %>% 
    mutate(
      TYPE = case_when(
        ############## LOF #############  
        (
          FILTER_SHORT_LOF &
            IMPACT == 'HIGH'
        ) ~ 'LOF',
        ########### MISSENSE ########### 
        (
          FILTER_SHORT_MISSENSE &
            str_detect(Consequence, "missense")
        ) ~ 'MISSENSE',
        ############# SPLICING ########### 
        (
          FILTER_SHORT_SPLICING & 
            (
              SpliceAI_max >= getOption('FILTER_SHORT_MIN_SPLICEAI_PP', -Inf) |
                (
                  str_detect(Consequence, "splice") &
                    IMPACT == 'MODERATE'
                ) 
            )
        ) ~ 'SPLICING',
        ############# OTHER ############# 
        # Captures anything else
        (
          IMPACT %in% FILTER_SHORT_VEP_IMPACTS |
            Consequence %in% FILTER_SHORT_VEP_CONSEQUENCES |
            str_detect(CLNSIG, FILTER_SHORT_CLINVAR_KEEP_PAT) |
            CADD  > getOption('FILTER_SHORT_MIN_CADD_PP', -Inf)
        ) ~ 'OTHER'
      )
    ) %>% 
    filter(!is.na(TYPE))
  
  return(
    VARIANTS_OUT %>% TRACK_REASON('FILTER_SHORT_TYPE')
  )
}

###################### FILTER_INHERITANCE ###############################
FILTER_INHERITANCE <- function(VARIANTS, PEDIGREE, set = 'SHORT') {
  stopifnot(
    require(tidyverse),
    is.data.frame(VARIANTS),
    is.data.frame(PEDIGREE),
    is.character(set)
  )
  
  pop_dom_max_af  <- getOption(str_c('FILTER_', set, '_POP_DOM_MAX_AF') , Inf)
  pop_rec_max_af  <- getOption(str_c('FILTER_', set, '_POP_REC_MAX_AF') , Inf)
  pop_dom_max_ac  <- getOption(str_c('FILTER_', set, '_POP_DOM_MAX_AC') , Inf)
  pop_rec_max_ac  <- getOption(str_c('FILTER_', set, '_POP_REC_MAX_AC') , Inf)
  pop_dom_max_hom <- getOption(str_c('FILTER_', set, '_POP_DOM_MAX_HOM'), Inf)
  pop_rec_max_hom <- getOption(str_c('FILTER_', set, '_POP_REC_MAX_HOM'), Inf)
  
  coh_dom_max_af  <- getOption(str_c('FILTER_', set, '_COH_DOM_MAX_AF') , Inf)
  coh_rec_max_af  <- getOption(str_c('FILTER_', set, '_COH_REC_MAX_AF') , Inf)
  coh_dom_max_ac  <- getOption(str_c('FILTER_', set, '_COH_DOM_MAX_AC') , Inf)
  coh_rec_max_ac  <- getOption(str_c('FILTER_', set, '_COH_REC_MAX_AC') , Inf)

  FILTER_SHORT_CLINVAR_KEEP_PAT <- getOption('FILTER_SHORT_CLINVAR_KEEP_PAT', '$.')

  VARIANTS_OUT <-
    VARIANTS %>% 
    mutate(
      GENOTYPE = 
        GT %>% 
        mutate_all(function(x) {
          case_when(
            # TODO: HEMI variants? (should just appear as HOM)
            str_count(x, '[01]') == 2 & str_count(x, '[1]') == 2 ~ 'HOM',
            str_count(x, '[01]') == 2 & str_count(x, '[1]') == 1 ~ 'HET',
            str_count(x, '[01]') == 2 & str_count(x, '[1]') == 0 ~ 'REF'
          ) %>% factor(c('REF', 'HET', 'HOM'))
        })
    )
  
  # Affected
  AFF <- 
    filter(PEDIGREE, affected) %>% 
    pull(id) %>% 
    intersect(colnames(VARIANTS_OUT$GENOTYPE))
  # Unaffected
  UNA <- 
    filter(PEDIGREE, !affected) %>% 
    pull(id) %>% 
    intersect(colnames(VARIANTS_OUT$GENOTYPE))
  
  VARIANTS_OUT <-
    VARIANTS_OUT %>% 
    mutate(
      inheritance = case_when(
        (
          rowSums(GENOTYPE[, AFF, drop = F] == 'HET') == length(AFF) &
            rowSums(GENOTYPE[, UNA, drop = F] == 'REF') == length(UNA)
        ) ~ 'dominant',
        (
          rowSums(GENOTYPE[, AFF, drop = F] == 'HOM') == length(AFF) &
            rowSums(GENOTYPE[, UNA, drop = F] != 'HOM') == length(UNA)
        ) ~ 'recessive',
      )) %>% 
    select(-GENOTYPE) %>% 
    mutate(
      inheritance = case_when(
         # Retain ClinVar pathogenic regardless of allele frequencies
         
        ( 
          inheritance == 'dominant' &
            pop_AF  < pop_dom_max_af &
            pop_AC  < pop_dom_max_ac &
            pop_hom < pop_dom_max_hom &
            AF      < coh_dom_max_af &
            AC      < coh_dom_max_ac 
        ) ~ 'dominant',
        ( 
          inheritance == 'recessive' &
            pop_AF  < pop_rec_max_af &
            pop_AC  < pop_rec_max_ac &
            pop_hom < pop_rec_max_hom &
            AF      < coh_rec_max_af &
            AC      < coh_rec_max_ac 
        ) ~ 'recessive',
        ( 
          inheritance == 'dominant' &
            pop_AF  < pop_rec_max_af &
            pop_AC  < pop_rec_max_ac &
            pop_hom < pop_rec_max_hom &
            AF      < coh_rec_max_af &
            AC      < coh_rec_max_ac 
        ) ~ 'compound',
        (
          inheritance == 'dominant' &
            str_detect(CLNSIG, FILTER_SHORT_CLINVAR_KEEP_PAT)
        ) ~ 'dominant',
        (
          inheritance == 'recessive' &
            str_detect(CLNSIG, FILTER_SHORT_CLINVAR_KEEP_PAT)
        ) ~ 'recessive',
      )
    ) %>% 
    filter(!is.na(inheritance))
  
  return(
    VARIANTS_OUT %>% TRACK_REASON(str_c('FILTER_INHERITANCE(', set, ')'))
  )
}

###################### FILTER_COMPOUND ###############################
FILTER_COMPOUND <- function(SHORT_VAR = NULL, STRUC_VAR = NULL) {
  stopifnot(
    require(tidyverse),
    is.data.frame(SHORT_VAR) | is.null(SHORT_VAR),
    is.data.frame(STRUC_VAR) | is.null(STRUC_VAR)
  )
  
  GENES <- character()
  
  if (!is.null(SHORT_VAR)) {
    GENES <- c(GENES, SHORT_VAR$Gene)
  }
  
  if (!is.null(STRUC_VAR)) {
    GENES <- c(GENES, STRUC_VAR$Gene)
  }
  
  # genes with multiple hits
  GENES <-
    tibble(GENES) %>% 
    count(GENES, name = 'n_hit') %>% 
    filter(n_hit > 1) %>% 
    pull()

  

  if (!is.null(SHORT_VAR)) {
    FILTER_SHORT_CLINVAR_KEEP_PAT <- getOption('FILTER_SHORT_CLINVAR_KEEP_PAT', '$.')

    SHORT_OUT <-
      SHORT_VAR %>% 
      filter(
        inheritance != 'compound' | 
          Gene %in% GENES | 
          str_detect(CLNSIG, FILTER_SHORT_CLINVAR_KEEP_PAT)
      )
    TRACK_REASON(SHORT_OUT, str_c('FILTER_COMPOUND(SHORT)'), set = 'SHORT')
  } else {
    SHORT_OUT <- NULL
  }
    
  if (!is.null(STRUC_VAR)) {
    STRUC_OUT <-
      STRUC_VAR %>% 
      filter(
        inheritance != 'compound' | Gene %in% GENES
      )
    TRACK_REASON(SHORT_OUT, str_c('FILTER_COMPOUND(STRUC)'), set = 'STRUC')
  } else {
    STRUC_OUT <- NULL
  }
    
  return(
    list(SHORT = SHORT_OUT, STRUC = STRUC_OUT)
  )
  
}

PRINT_STR <- function(x) {
  capture.output(str(x)) %>% 
    str_c(collapse = '\n') %>% 
    message()
}

doc <- "
Usage:
  filter.R <ped> <gene_set> <options> ... [options]

Options:
  ped                         Pedigree file.
  gene_set                    File container ensembl gene ids of intererst.
  options                     Filter options Json file.
  --short-var=<TSV>           Short Variants TSV input.
  --struc-var=<TSV>           Structural Variants TSV input.
  --output=<prefix>           Output file prefix [default: output].
  --cav-opts=<json>           Json file with additional options for cavalier package.
"

# run main function
invisible(MAIN(docopt(doc)))