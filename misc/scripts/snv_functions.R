
# All functions have access to the following arguments:
# - FILEPATH : path to input TSV from split-vep of VCF
# - VARIANTS : data frame with variants, one row per gene/allele/variant (except first function, e.g. SNV_LOAD)
# - CONFIG   : customisable configuration options
# - GENE_LIST: a data frame with columns ensembl_gene_id and hgnc_symbol for gene list filtering
# - PEDIGREE : a pedigree data frame, columns: famid, id, dadid, momid, sex, affected
#
# NB: functions must use named arguments, not positional, and include "..." at end


###################### SNV_LOAD ###############################
SNV_LOAD <- function(FILEPATH = NULL, PEDIGREE = NULL, ...) {
  # check libraries and args
  stopifnot(
    require(tidyverse),
    rlang::is_scalar_character(FILEPATH),
    file.exists(FILEPATH),
    is.data.frame(PEDIGREE)
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
    FMT_AD_T23547 = col_character(),
    FMT_AD_T23548 = col_character(),
    FMT_GT_T23547 = col_character(),
    FMT_GT_T23548 = col_character(),
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
    CLIN_SIG = col_character(),
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
    # add/modify columns
    mutate(
      # useful for filtering
      pop_AF = pmax(replace_na(gnomad_AF, 0), replace_na(gnomad_FAF95, 0)),
      pop_AC =  replace_na(gnomad_AC, 0),
      SpliceAI_max = 
        replace_na(
          pmax(SpliceAI_pred_DS_AG, SpliceAI_pred_DS_AL, SpliceAI_pred_DS_DG, SpliceAI_pred_DS_DL),
          0
        ),
      SIFT_score     = str_extract(SIFT    , '(?<=\\()[0-9\\.]+(?=\\)$)') %>% as.numeric(),
      PolyPhen_score = str_extract(PolyPhen, '(?<=\\()[0-9\\.]+(?=\\)$)') %>% as.numeric(),
      # not meaningful for insertions
      phyloP100 = replace(phyloP100, VARIANT_CLASS == "insertion", NA),
      # title for slides
      title = str_c(PEDIGREE$famid[1], SYMBOL, HGVSg, sep = " - "),
      # reporting summary columns
      gnomAD = str_c("AF=", signif(replace_na(gnomad_AF, 0), 2), "; AC=", replace_na(gnomad_AC, 0), "; Hom=",  replace_na(gnomad_nhomalt, 0)),
      cohort = str_c("AF=", signif(replace_na(AF, 0), 2), "; AC=", replace_na(AC, 0)),
      SpliceAI = str_c("AG=", SpliceAI_pred_DS_AG, "; DG=", SpliceAI_pred_DS_DG, "; AL=", SpliceAI_pred_DS_AL, "; DL=", SpliceAI_pred_DS_DL),
      AlphaMissense = str_c(am_class, "(", am_pathogenicity, ")"),
      # add whitespace for text wrapping
      HGVSc = str_replace(HGVSc, ":", ": "),
      HGVSp = str_replace(HGVSp, ":", ": "),
    )
  
  return(VARIANTS_OUT)
}

###################### SNV_FILTER_DEPTH ###############################
SNV_FILTER_DEPTH <- function(VARIANTS = NULL, CONFIG = NULL, ...) {
  stopifnot(
    require(tidyverse),
    is.data.frame(VARIANTS),
    is.list(CONFIG),
    rlang::is_scalar_integerish(CONFIG$min_depth),
    any(str_starts(colnames(VARIANTS), 'FMT_AD_'))
  )
  
  VARIANTS_OUT <-
    VARIANTS %>% 
    mutate(
      min_depth = 
        pick(starts_with('FMT_AD_')) %>% 
        mutate_all(function(x) map_int( str_split(x, ","), ~ sum(as.integer(.x)) )) %>% 
        (function(x) do.call(pmin, unname(x)))) %>% 
    filter(min_depth >= CONFIG$min_depth) %>% 
    select(-min_depth)
  
  return(VARIANTS_OUT)
}

###################### SNV_FILTER_GENES ###############################
SNV_FILTER_GENES <- function(VARIANTS = NULL, GENE_LIST = NULL, ...) {
  stopifnot(
    require(tidyverse),
    is.data.frame(VARIANTS),
    is.data.frame(GENE_LIST)
  )
  
  VARIANTS_OUT <-
    VARIANTS %>%
    filter(
      Gene %in% na.omit(GENE_LIST$ensembl_gene_id) |
        SYMBOL %in% na.omit(GENE_LIST$hgnc_symbol)
    )
  
  return(VARIANTS_OUT)
}

###################### SNV_FILTER_TYPE ###############################
SNV_FILTER_TYPE <- function(VARIANTS = NULL, ...) {
  stopifnot(
    require(tidyverse),
    is.data.frame(VARIANTS)
  )
  
  VARIANTS_OUT <-
    VARIANTS %>% 
    mutate(
      TYPE = case_when(
        ############## LOF #############  
        (
          IMPACT == 'HIGH'
        ) ~ 'LOF',
        ########### MISSENSE ########### 
        (
          str_detect(Consequence, "missense") &
            (
                CADD  > 20                           |
                REVEL > 0.5                          |
                str_detect(CLIN_SIG , "pathogenic" ) |
                str_detect(SIFT     , "deleterious") |
                !str_detect(PolyPhen, "benign"     ) |
                !str_detect(am_class, "benign"     )
            )
        ) ~ 'MISSENSE',
        ############# SPLICING ########### 
        (
          SpliceAI_max > 0.20 |
          (
            str_detect(Consequence, "splice") &
            IMPACT == 'MODERATE'
          ) 
        ) ~ 'SPLICING',
        ############# OTHER ############# 
        (
          (
            IMPACT == 'MODERATE' & 
            !str_detect(Consequence, "missense")
          )                                  |
          str_detect(CLIN_SIG, "pathogenic") |
          CADD      > 25                     |
          phyloP100 > 8
        ) ~ 'OTHER'
      )
    ) %>% 
    filter(!is.na(TYPE))
  
  return(VARIANTS_OUT)
}

###################### SNV_FILTER_INHERITANCE ###############################
SNV_FILTER_INHERITANCE <- function(VARIANTS = NULL, CONFIG = NULL, PEDIGREE = NULL, ...) {
  stopifnot(
    require(tidyverse),
    is.data.frame(VARIANTS),
    is.list(CONFIG),
    is.list(CONFIG$freq_thresholds),
    is.data.frame(PEDIGREE)
  )
  
  val_or_inf <- function(expr) {
    expr_sub <- substitute(expr)
    x <-
      tryCatch(
        eval(expr_sub, envir = parent.frame()),
        error = function(e) { NULL }
      )
    if (is.numeric(x)) {
      return(x)
    }
    return(Inf)
  }
  
  VARIANTS_OUT <-
    VARIANTS %>% 
    mutate(
      ## currently required by create_slides - should probably lift this requirement
      genotype = 
        pick(starts_with('FMT_GT_')) %>% 
        rename_with(~str_remove(., 'FMT_GT_')),
      GENOTYPE = 
        genotype %>% 
        mutate_all(function(x) {
          case_when(
            # TODO - hemizygous variants
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
      Inheritance = case_when(
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
    filter(!is.na(Inheritance)) %>% 
    (function(x) bind_rows(
      x,
      # "dominant" could also be compound het
      filter(x, Inheritance == 'dominant') %>% mutate(Inheritance = 'compound')
    )) %>% 
    filter(
      ########### pop_AF ##########
      (
        (Inheritance %in% c('recessive', 'compound') & pop_AF < val_or_inf(CONFIG$freq_thresholds$pop$AF$recessive)) |
        (Inheritance == 'dominant'                   & pop_AF < val_or_inf(CONFIG$freq_thresholds$pop$AF$dominant))
      ),
      ########### pop_AC ##########
      (
        (Inheritance %in% c('recessive', 'compound') & pop_AC < val_or_inf(CONFIG$freq_thresholds$pop$AC$recessive)) |
        (Inheritance == 'dominant'                   & pop_AC < val_or_inf(CONFIG$freq_thresholds$pop$AC$dominant))
      ),
      ########### AF ############
      (
        (Inheritance %in% c('recessive', 'compound') & AF     < val_or_inf(CONFIG$freq_thresholds$cohort$AF$recessive)) |
        (Inheritance == 'dominant'                   & AF     < val_or_inf(CONFIG$freq_thresholds$cohort$AF$dominant))
      ),
      ########### AC ############
      (
        (Inheritance %in% c('recessive', 'compound') & AC     < val_or_inf(CONFIG$freq_thresholds$cohort$AC$recessive)) |
        (Inheritance == 'dominant'                   & AC     < val_or_inf(CONFIG$freq_thresholds$cohort$AC$dominant))
      )
    ) %>% 
    chop(Inheritance) %>% 
    mutate(
      Inheritance = map(Inheritance, unique) %>% map_chr(str_c, collapse = '&'),
      Inheritance = if_else(str_detect(Inheritance, 'dominant'), 'dominant', Inheritance)
    ) %>% 
    add_count(Gene, dom_comp = Inheritance %in% c('dominant', 'compound'), name = 'n_compound') %>%
    filter(Inheritance != 'compound' | n_compound > 1) %>% 
    select(-n_compound, -dom_comp)
  
  return(VARIANTS_OUT)
}

