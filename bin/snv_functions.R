
# All functions have access to the following arguments:
# - FILEPATH : path to input TSV from split-vep of VCF
# - VARIANTS : data frame with variants, one row per gene/allele/variant (except first function, e.g. SNV_LOAD)
# - CONFIG   : customisable configuration options
# - GENE_LIST: a data frame with columns ensembl_gene_id and hgnc_symbol for gene list filtering
# - PEDIGREE : a pedigree data frame, columns: famid, id, dadid, momid, sex, affected
#
# NB: functions must use named arguments, not positional, and include "..." at end

####################### FMT_TO_MATRIX #########################
# helper to convert sample format fields to a dataframe with colnames
FMT_TO_DF <- function(input_data) {
  stopifnot(
    require(tidyverse),
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
          readr::type_convert(guess_integer = T, col_types = cols())
      ) %>% 
      select(-starts_with(prefix))
  }
  
  return(output)
}

###################### SNV_LOAD ###############################
SNV_LOAD <- function(FILEPATH = NULL, ...) {
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
    CLIN_SIG = col_character(),
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
      AN = AN - (GT %>% mutate_all(~ str_count(., '[01]')) %>% rowSums(na.rm = TRUE)),
      AC = AC - (GT %>% mutate_all(~ str_count(., '[1]')) %>% rowSums(na.rm = TRUE)),
      AF = AC / AN,
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
      phyloP100 = replace(phyloP100, VARIANT_CLASS == "insertion", NA)      
    )
  
  return(VARIANTS_OUT)
}

###################### SNV_FILTER_FMT ###############################
SNV_FILTER_FMT <- function(VARIANTS = NULL, CONFIG = NULL, ...) {
  stopifnot(
    require(tidyverse),
    is.data.frame(VARIANTS),
    is.list(CONFIG)
  )
  
  VARIANTS_OUT <- VARIANTS 
  
  if ('min_DP' %in% names(CONFIG)) {
    stopifnot( 'DP' %in% colnames(VARIANTS) )
    VARIANTS_OUT <-
      VARIANTS_OUT %>% 
      filter(do.call(pmin, unname(DP)) >= CONFIG$min_DP)
  }
  
  if ('min_GQ' %in% names(CONFIG)) {
    stopifnot( 'GQ' %in% colnames(VARIANTS) )
    VARIANTS_OUT <-
      VARIANTS_OUT %>% 
      filter(do.call(pmin, unname(GQ)) >= CONFIG$min_GQ)
  }
  
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
                str_detect(CLIN_SIG , "pathogenic") |
                CADD             > 22.7             | # Clingen BP4 supporting  doi: 10.1016/j.ajhg.2022.10.013
                REVEL            > 0.29             | # Clingen BP4 supporting  doi: 10.1016/j.ajhg.2022.10.013
                SIFT_score       < 0.327            | # Clingen BP4 supporting  doi: 10.1016/j.ajhg.2022.10.013
                PolyPhen_score   > 0.113            | # Clingen BP4 supporting  doi: 10.1016/j.ajhg.2022.10.013
                am_pathogenicity > 0.169              # Clingen BP4 supporting  doi: https://doi.org/10.1016/j.gim.2025.101402
            )
        ) ~ 'MISSENSE',
        ############# SPLICING ########### 
        (
          SpliceAI_max >= 0.20 | # Clingen PP3 supporting https://doi.org/10.1016/j.ajhg.2023.06.002
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
          CADD      > 25.3   # Clingen PP3 supporting  doi: 10.1016/j.ajhg.2022.10.013
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
      GENOTYPE = 
        GT %>% 
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
      ########### pop_hom ##########
      (
        (Inheritance %in% c('recessive', 'compound') & pop_hom < val_or_inf(CONFIG$freq_thresholds$pop$hom$recessive)) |
        (Inheritance == 'dominant'                   & pop_hom < val_or_inf(CONFIG$freq_thresholds$pop$hom$dominant))
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
    select(-n_compound, -dom_comp) %>% 
    group_by(Gene) %>% 
    mutate(Inheritance = case_when(
      Inheritance == 'dominant' & n() > 1 ~ 'dominant_or_compound',
      TRUE ~ Inheritance
    )) %>% 
    ungroup()
  
  return(VARIANTS_OUT)
}

SNV_REPORT <- function(VARIANTS = NULL, PEDIGREE = NULL, ...) {
  # check libraries and args
  stopifnot(
    require(tidyverse),
    is.data.frame(PEDIGREE),
    is.data.frame(VARIANTS)
  )
  
  VARIANTS_OUT <-
    VARIANTS  %>% 
    # add/modify columns
    mutate(
      # reporting summary columns
      # title for slides
      broad_id = str_c(CHROM, POS, REF, ALT, sep = '-'),
      title = str_c(PEDIGREE$famid[1], SYMBOL, broad_id, sep = " - "),
      gnomAD = str_c("AF=", signif(replace_na(gnomad_AF, 0), 2), "; AC=", replace_na(gnomad_AC, 0), "; Hom=",  replace_na(gnomad_nhomalt, 0)),
      Cohort = str_c("AF=", signif(replace_na(AF, 0), 2), "; AC=", replace_na(AC, 0)),
      SpliceAI = str_c("AG=", SpliceAI_pred_DS_AG, "; DG=", SpliceAI_pred_DS_DG, "; AL=", SpliceAI_pred_DS_AL, "; DL=", SpliceAI_pred_DS_DL),
      AlphaMissense = str_c(am_class, "(", am_pathogenicity, ")"),
      dbSNP = str_extract(Existing_variation, 'rs[0-9]+'),
      # add GRCh38 for mutatylzer compatibility
      HGVSg = str_replace(HGVSg, "^([^:]+):(.*)$", "GRCh38(\\1):\\2"),
      # prefer coding & protein, or coding, or genomic. Drop IDs to better fit
      HGVS = coalesce(
        str_c(str_remove(HGVSc, '^.+:(?=[cp])'), '; ', str_remove(HGVSp, '^.+:(?=[cp])')),
        str_remove(HGVSc, '^.+:(?=[cp])'),
        HGVSg
      ),
      # # add URLS to slides
      gnomAD_url = str_c(
        'https://gnomad.broadinstitute.org/variant/',
        URLencode(broad_id),
        '?dataset=gnomad_r4'
      ),
      HGVS_url = str_c(
        'https://mutalyzer.nl/normalizer/',
        URLencode(HGVSg)
      ),
      SpliceAI_url = str_c(
        'https://spliceailookup.broadinstitute.org/#hg=38&variant=',
        URLencode(broad_id)
      ),
      Gene_url = str_c(
        'https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=',
        URLencode(Gene)
      ),
      CLIN_SIG_url = str_c(
        "https://www.ncbi.nlm.nih.gov/clinvar/?term=",
        dbSNP
      )
    )
  
  return(VARIANTS_OUT)
}

