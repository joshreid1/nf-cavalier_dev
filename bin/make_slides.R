#!/usr/bin/env Rscript

stopifnot(
  require(tidyverse),
  require(magrittr),
  require(docopt),
  require(cavalier)
)

MAIN <- function(opts) {
  
  setwd('/vast/scratch/users/munro.j/nextflow/work/22/f7f951e306d2e1450232ced94176d4')
  
  opts <- list(
    ped = 'Haf_PK271.ped',
    gene_lists =  'PAA_202.a9ffc3a8a68443fe561c4f3f820ed628.tsv',
    short_var = 'Haf_PK271.short.filtered_variants.rds',
    struc_var = 'Haf_PK271.struc.filtered_variants.rds',
    igv       = 'SID_T21432.VID_chr6-170283291-G-A.png,SID_T21437.VID_chr6-170283291-G-A.png',
    svpv      = 'SVPV_chr19.13323531.DEL.95499b5cc7.png,SVPV_chr19.13323546.DEL.95499b5cc7.png,SVPV_chr19.13324001.DEL.95499b5cc7.png',
    samplot   = 'samplot_Haf_PK271_chr19-13323531-DEL-21293.png,samplot_Haf_PK271_chr19-13323546-DEL-21281.png,samplot_Haf_PK271_chr19-13324001-DEL-21000.png',
    slide_info = 'slide_info.json', 
    cav_opts  = 'cavalier_options.e044d87757d209315972f2caa5d05f3c.json',
    output    = 'test'
  )

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
  svpv_pngs <- c(str_split(opts$svpv, ',', simplify = T))
  samplots  <- c(str_split(opts$samplot, ',', simplify = T))
  
  ################# CHECK ARGS ######################
  stopifnot(
    opts$short_var == 'NONE' | file.exists(opts$short_var),
    opts$struc_var == 'NONE' | file.exists(opts$struc_var),
    file.exists(opts$ped),
    all(file.exists(gene_lists)),
    all(igv_pngs  == 'NONE') | all(file.exists(igv_pngs)),
    all(svpv_pngs == 'NONE') | all(file.exists(svpv_pngs)),
    all(samplots  == 'NONE') | all(file.exists(samplots)),
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
  
  ########### LOAD VARIANTS ######################
  empty_cand <- tibble(SYMBOL = character(), Gene = character(), variant_id = character(), GT = list())
  
  short_cand <- empty_cand
  if (opts$short_var != 'NONE') {
    short_cand <- readRDS(opts$short_var) 
  }
  
  struc_cand <- empty_cand
  if (opts$struc_var != 'NONE') {
    struc_cand <- readRDS(opts$struc_var) 
  }
  
  ########### SHORT+STRUC SHARED FEATURES ###################
  ########### GTEX PLOTS ########################
  GTEX <- 
    bind_rows(
      short_cand %>% select(SYMBOL, Gene),
      struc_cand %>% select(SYMBOL, Gene),
    ) %>% 
    distinct() %>% 
    mutate(GTEX = map2(SYMBOL, Gene, cavalier::plot_gtex_expression)) %>% 
    select(Gene, GTEX) %>% 
    group_by(Gene) %>% 
    slice(1)
  
  ########### OMIM TABLE ########################
  OMIM <- 
    bind_rows(
      short_cand %>% select(Gene),
      struc_cand %>% select(Gene),
    ) %>% 
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
    semi_join( 
      bind_rows(
        short_cand %>% select(Gene),
        struc_cand %>% select(Gene),
      ),
      by = 'Gene'
    ) %>% 
    nest(LISTS = -Gene)
  
  ########### PEDIGREE PLOTS ########################
  if (nrow(pedigree) > 2) {
    PEDIGREE <-
      bind_rows(
        struc_cand %>% select(variant_id, GT),
        short_cand %>% select(variant_id, GT)
      ) %>% 
      distinct() %>% 
      unnest(GT) %>%
      pivot_longer(-variant_id, names_to = 'id', values_to = 'gt') %>% 
      right_join(.,
                 expand_grid(pedigree, variant_id = unique(.$variant_id),
                             by = 'id')
      ) %>% 
      mutate(gt = replace_na(gt, '???'),
             label = str_c(id, gt, sep = '\n')) %>% 
      nest(PEDIGREE = -variant_id) %>% 
      mutate(PEDIGREE = map(PEDIGREE, distinct)) %>%
      mutate(PEDIGREE = map(PEDIGREE, plot_ped))
  } else {
    PEDIGREE <-
      bind_rows(
        struc_cand %>% select(variant_id, GT),
        short_cand %>% select(variant_id, GT)
      ) %>% 
      distinct() %>% 
      mutate(PEDIGREE = map(variant_id, ~ NULL))
  }
  
  ########### CREATE SHORT SLIDES ##############
  if (opts$short_var != 'NONE') { 
    
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
    
    ########### SHORT VAR INFO TABLE ########################
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
    
    SHORT_VAR_INFO <-
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
    
    
    ########### SHORT SLIDE DATA ########################
    SHORT_SLIDE_DATA <-
      short_cand %>% 
      select(Gene, SYMBOL, variant_id) %>% 
      distinct() %>% 
      mutate(TITLE = str_c(opts$output, SYMBOL, variant_id, sep = ' - ')) %>% 
      select(TITLE, Gene, variant_id) %>% 
      left_join(SHORT_VAR_INFO, by = c('Gene', 'variant_id')) %>% 
      left_join(IGV, by = 'variant_id') %>% 
      left_join(PEDIGREE, by = 'variant_id') %>% 
      left_join(OMIM, by = 'Gene') %>% 
      left_join(LISTS, by = 'Gene') %>% 
      left_join(GTEX, by = 'Gene') %>% 
      arrange(TITLE)
    
    ############# SHORT SLIDE LAYOUT #####################  
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
        slide_data = SHORT_SLIDE_DATA,
        output = str_c(opts$output, '.pptx')
      )
  }
  
  if (opts$struc_var != 'NONE') { 
    ########### SVPV Plots ########################
    SVPV <-
      tibble( png = svpv_pngs) %>% 
      tidyr::extract(
        png,
        into = c("CHROM", "POS", "SVTYPE"),
        regex = "^SVPV_([^\\.]+)\\.([0-9]+)\\.([^\\.]+)",
        remove = FALSE
      ) %>% 
      mutate(POS = as.numeric(POS)) %>% 
      inner_join(
        struc_cand %>% 
          select(CHROM, POS, SVTYPE, variant_id) %>% 
          distinct(),
      ) %>% 
      transmute(
        id = 'x',
        png,
        variant_id
      ) %>% 
      nest(SVPV = -variant_id) %>% 
      mutate(SVPV = 
               SVPV %>% 
               map(cavalier::plot_png_facets, crop_left = 0, crop_right = 0) %>% 
               map(~ . + ggplot2::theme(strip.text = element_blank())) # kludge
      )
    
    ########### SAMPLOTS ########################
    SAMPLOT <-
      tibble(png = samplots) %>% 
      tidyr::extract(
        png,
        into = c("CHROM", "POS", "SVTYPE", "SVLEN"),
        regex = "_([^\\._]+)-([0-9]+)-([^\\.]+)-([_0-9]+)\\.png",
        remove = FALSE
      ) %>% 
      mutate(POS = as.numeric(POS), SVLEN = as.numeric(SVLEN)) %>% 
      inner_join(
        struc_cand %>% 
          select(CHROM, POS, SVTYPE, SVLEN, variant_id) %>% 
          mutate(SVLEN = abs(SVLEN)) %>% 
          distinct(),
      ) %>% 
      transmute(
        id = 'x',
        png,
        variant_id
      ) %>% 
      nest(SAMPLOT = -variant_id) %>% 
      mutate(SAMPLOT = 
               SAMPLOT %>% 
               map(cavalier::plot_png_facets, crop_left = 0, crop_right = 0) %>% 
               map(~ . + ggplot2::theme(strip.text = element_blank())) # kludge
      )
    
    ########### STRUC VAR INFO TABLE ########################
    
    STRUC_FIELDS_ALL <- 
      reduce(slide_info$STRUC, c) %>% 
      (function(x) x[unique(names(x))]) %>% 
      unlist()
    
    # STRUC_FIELDS <-
    #   struc_cand %>% 
    #   select(TYPE) %>% 
    #   distinct() %>% 
    #   mutate(FIELDS = map(TYPE, function(x) {
    #     fields <- names(slide_info$SHORT$DEFAULT)
    #     custom <- names(slide_info$SHORT[[x]])
    #     union(fields, custom)
    #   })) %>% 
    #   with(setNames(FIELDS, TYPE))
    
    STRUC_VAR_INFO <-
      struc_cand %>% 
      # add/modify columns
      mutate(
        # need to maintain Gene
        Gene2 = Gene,
        # reporting summary columns
        title = str_c(opts$output, SYMBOL, variant_id, sep = " - "),
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
  }
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
  --struc-var=<RDS>           Structural Variants RDS input.
  --igv=<PNG>                 IGV screenshot PNGs, comma separated.
  --svpv=<PNG>                SVPV PNGs, comma separated.
  --samplot=<PNG>             samplot PNGs, comma separated.
  --output=<prefix>           Output file prefix [default: output].
  --slide-info=<json>         Json file specifying fields to include in VAR_INFO.
  --cav-opts=<json>           Json file with additional options for cavalier package.
"

# run main function
# invisible(MAIN(docopt(doc)))
