#!/usr/bin/env Rscript

stopifnot(
  require(tidyverse),
  require(magrittr),
  require(docopt),
  require(cavalier)
)

MAIN <- function(opts) {
  
  # opts <- list(
  #   ped = '/vast/scratch/users/munro.j/nextflow/work/27/ca82f64091092582cd2037898088bd/Plu_PK86442.ped',
  #   gene_lists = '/vast/scratch/users/munro.j/nextflow/work/2d/ef2a22c6f185fcd1d1a75c8caf1a0d/G4E_ALL_v2025-09.tsv',
  #   slide_info = '/vast/scratch/users/munro.j/nextflow/work/a3/3bca35cd4a32408c99dddeb26ba512/slide_info.json',
  #   short_var = '/vast/scratch/users/munro.j/nextflow/work/27/ca82f64091092582cd2037898088bd/Plu_PK86442.short.filtered_variants.rds',
  #   igv = '/vast/scratch/users/munro.j/nextflow/work/27/ca82f64091092582cd2037898088bd/SID_T23547.VID_chr6-31668630-T-TA.png,/vast/scratch/users/munro.j/nextflow/work/27/ca82f64091092582cd2037898088bd/SID_T23548.VID_chr6-31668630-T-TA.png',
  #   cav_opts = '/vast/scratch/users/munro.j/nextflow/work/2d/ef2a22c6f185fcd1d1a75c8caf1a0d/cavalier_options.cache.json',
  #   output  = 'test'
  # )

  message('Using CLI options:')
  PRINT_STR(
    opts[str_detect(names(opts), '^[:alpha:]')] 
  )

  ############### SET OPTIONS ######################
  slide_info <- jsonlite::fromJSON(opts$slide_info)
  message('Using slide info:')
  PRINT_STR(slide_info)
  
  if (!is.null(opts$cav_opts)) {
    cavalier::set_options_from_json(opts$cav_opts)
  }
  
  gene_lists <- c(str_split(opts$gene_lists, ',', simplify = T))
  igv_pngs <- c(str_split(opts$igv, ',', simplify = T))
  
  ################# CHECK ARGS ######################
  stopifnot(
    file.exists(opts$short_var),
    file.exists(opts$ped),
    all(file.exists(gene_lists)),
    all(file.exists(igv_pngs)),
    all(c('SHORT', 'STRUC') %in% names(slide_info)),
    'DEFAULT' %in% names(slide_info$SHORT),
    'DEFAULT' %in% names(slide_info$STRUC)
  )
  
  ############## READ PEDIGREE   ################
  
  pedigree <- cavalier::read_ped(opts$ped)
  
  ############## CHECK GENE LIST ###############
  
  gene_list <- 
    c(str_split(opts$gene_lists, ',', simplify = T)) %>% 
    map_df(read_tsv, col_types = cols(.default = col_character()))
  
  stopifnot(
    'ensembl_gene_id' %in% colnames(gene_list),
    nrow(filter(gene_list, !is.na(ensembl_gene_id))) > 0
  )
  
  if (nrow(filter(gene_list, is.na(ensembl_gene_id)))) {
    warning('Excluded ', nrow(filter(gene_list, is.na(ensembl_gene_id))), ' entries from Gene list with missing ensembl_gene_id')
  }

  short_cand <- readRDS(opts$short_var) 
  
  ########### IGV Plots ########################
  IGV <-
    tibble(
      png = igv_pngs,
      id = str_extract(png, '(?<=SID_).+(?=\\.VID_)'),
      variant_id = str_extract(png, '(?<=VID_).+(?=\\.png$)'),
    ) %>% 
    inner_join(
      short_cand %>% 
        select(variant_id, GT) %>% 
        unnest(GT) %>%
        pivot_longer(-variant_id, names_to = 'id', values_to = 'GT'),
      by = c('variant_id', 'id')
    ) %>% 
    transmute(
      id = str_c(id, ': ', GT),
      png,
      variant_id
    ) %>% 
    nest(IGV = -variant_id) %>% 
    mutate(IGV = map(IGV, cavalier::plot_png_facets, crop_left = 70, crop_right = 120))
  
  ########### PEDIGREE PLOTS ########################
  if (nrow(pedigree) > 2) {
    PEDIGREE <-
      short_cand %>% 
      select(variant_id, GT) %>% 
      unnest(GT) %>%
      pivot_longer(-variant_id, names_to = 'id', values_to = 'gt') %>% 
      right_join(.,
                 expand_grid(pedigree, variant_id = unique(.$variant_id),
                             by = 'id')
      ) %>% 
      mutate(gt = replace_na(gt, '???'),
             label = str_c(id, gt, sep = '\n')) %>% 
      nest(PEDIGREE = -variant_id) %>% 
      mutate(PEDIGREE = map(PEDIGREE, plot_ped))
  } else {
    PEDIGREE <-
      short_cand %>% 
      select(variant_id) %>% 
      distinct() %>% 
      mutate(PEDIGREE = map(variant_id, ~ NULL))
  }
  
  ########### VAR_INFO TABLE ########################
  SHORT_FIELDS_ALL <- 
    reduce(slide_info$SHORT, c) %>% 
    (function(x) x[unique(names(x))]) %>% 
    unlist()
  
  SHORT_FIELDS <-
    short_cand %>% 
    select(TYPE) %>% 
    distinct() %>% 
    mutate(FIELDS = map(TYPE, function(x) {
      fields <- names(slide_info$SHORT$DEFAULT)
      custom <- names(slide_info$SHORT[[x]])
      union(fields, custom)
    })) %>% 
    with(setNames(FIELDS, TYPE))
  
  VAR_INFO <-
    short_cand %>% 
    # add/modify columns
    mutate(
      # need to maintain Gene
      Gene2 = Gene,
      # reporting summary columns
      broad_id = str_c(CHROM, POS, REF, ALT, sep = '-'),
      title = str_c(opts$output, SYMBOL, broad_id, sep = " - "),
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
      CLNSIG_url = str_c(
        "https://www.ncbi.nlm.nih.gov/clinvar/variation/",
        CLNVID
      ),
      GTEX = map2(SYMBOL, Gene, cavalier::plot_gtex_expression),
    ) %>% 
    select(
      variant_id,
      Gene2,
      TYPE,
      all_of(SHORT_FIELDS_ALL),
      any_of(
        setNames(
          str_c(SHORT_FIELDS_ALL, '_url'),
          str_c(names(SHORT_FIELDS_ALL), '_url')
        )
      )
    ) %>% 
    nest(VAR_INFO = -c(variant_id, Gene2, TYPE)) %>% 
    rename(Gene = Gene2) %>% 
    mutate(VAR_INFO = map2(VAR_INFO, TYPE, function(VI, TYPE_) {
      select(VI, all_of(SHORT_FIELDS[[TYPE_]]), any_of(str_c(SHORT_FIELDS[[TYPE_]], '_url')))
    }))
  
  ########### GTEX PLOTS ########################
  GTEX <- 
    short_cand %>% 
    select(SYMBOL, Gene) %>% 
    distinct() %>% 
    mutate(GTEX = map2(SYMBOL, Gene, cavalier::plot_gtex_expression)) %>% 
    select(Gene, GTEX) %>% 
    group_by(Gene) %>% 
    slice(1)
  
  ########### OMIM TABLE ########################
  OMIM <- 
    short_cand %>% 
    select(variant_id, Gene) %>% 
    distinct() %>% 
    inner_join(
      get_omim_disease_map() %>% 
        select(Gene = ensembl_id, disease_id, inheritance, disease_name),
      by = 'Gene'
    ) %>% 
    distinct() %>% 
    transmute(
      Gene,
      OMIM = str_c(disease_name, ' [', replace_na(inheritance, '-'), ']'),
      OMIM_url = str_c('https://omim.org/entry/', str_extract(disease_id, '\\d+$'))
    ) %>% 
    nest(OMIM = -Gene)
  
  ########### GENE LIST TABLE ########################
  LISTS <-
    gene_list %>% 
    mutate(list_name_url = case_when(
      str_starts(list_id, 'PAA:') ~ str_c('https://panelapp-aus.org/panels/',
                                          str_extract(list_id, '(?<=PAA:)\\d+')),
      str_starts(list_id, 'PAE:') ~ str_c('https://panelapp.genomicsengland.co.uk/panels/',
                                          str_extract(list_id, '(?<=PAE:)\\d+')),
      str_starts(list_id, 'HP:') ~ str_c('https://hpo.jax.org/app/browse/term/',
                                         list_id),
      str_starts(list_id, 'G4E:') ~ 'https://bahlolab.github.io/Genes4Epilepsy/'
    )) %>%
    select(Gene = ensembl_gene_id, list_name, list_name_url, list_version, any_of ('inheritance')) %>% 
    semi_join(short_cand, by = 'Gene') %>% 
    nest(LISTS = -Gene)
    
      
  ########### SLIDE_DATA ########################
  SLIDE_DATA <-
    short_cand %>% 
    select(Gene, SYMBOL, variant_id) %>% 
    distinct() %>% 
    mutate(TITLE = str_c(opts$output, SYMBOL, variant_id, sep = ' - ')) %>% 
    select(TITLE, Gene, variant_id) %>% 
    left_join(VAR_INFO, by = c('Gene', 'variant_id')) %>% 
    left_join(IGV, by = 'variant_id') %>% 
    left_join(PEDIGREE, by = 'variant_id') %>% 
    left_join(OMIM, by = 'Gene') %>% 
    left_join(LISTS, by = 'Gene') %>% 
    left_join(GTEX, by = 'Gene') %>% 
    arrange(TITLE)
  
  ############# SLIDE_LAYOUT #####################  
  if (ncol(short_cand$GT) == 1) {
    layout <-
      cavalier::slide_layout(
        c('VAR_INFO', 'IGV', 'GTEX'),
        c('OMIM', 'LISTS'),
        heights = c(23, 8),
        title_height = 0.09,
        pad = 0.015,
        transpose = 'VAR_INFO'
      )
  } else {
    layout <-
      bind_rows(
        cavalier::slide_layout(
          c('VAR_INFO', 'PEDIGREE', 'GTEX'),
          c('OMIM', 'LISTS'),
          heights = c(23, 8),
          title_height = 0.09,
          pad = 0.015,
          transpose = 'VAR_INFO'
        ),
        cavalier::slide_layout(
          c('IGV'),
          slide_num = 2L,
          title_height = 0.09,
          pad = 0.015
        )
      )
  }
  
  slides <-
    cavalier::create_slides(
      slide_layout = layout,
      slide_data = SLIDE_DATA,
      output = str_c(opts$output, '.pptx')
    )
}


PRINT_STR <- function(x) {
  capture.output(str(x)) %>% 
    str_c(collapse = '\n') %>% 
    message()
}

doc <- "
Usage:
  make_slides.R <ped> <gene_lists> <options> ... [options]

Options:
  ped                         Pedigree file.
  gene_lists                  Comma serparated list of gene list filenames.
  options                     Slide options Json file.
  --short-var=<RDS>           Short Variants RDS input.
  --igv=<PNG>                 IGV screenshot PNGs, comma separated.
  --struc-var=<RDS>           Structural Variants RDS input.
  --output=<prefix>           Output file prefix [default: output].
  --slide-info=<json>         Json file specifying fields to include in VAR_INFO.
  --cav-opts=<json>           Json file with additional options for cavalier package.
"

# run main function
invisible(MAIN(docopt(doc)))
