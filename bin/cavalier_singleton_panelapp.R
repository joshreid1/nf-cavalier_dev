#!/usr/bin/env Rscript

stopifnot(require(tidyverse),
          require(docopt),
          require(cavalier))

doc <- "
Usage:
  cavalier_singleton_panelapp.R <vcf> <out> <sample_bam> [options]

Options:
  vcf                         Input VEP vcf file.
  out                         Output directory.
  sample_bam                  Sample name and bam file in format name=/path/to/bam.
  --genome=<f>                Reference genome for IGV snapshot [default: hg19].
  --gene-lists=<f>            Comma serparated list of gene list names.
  --maf-dom=<f>               Maximum MAF for dominant [default: 0.0001].
  --maf-rec=<f>               Maximum MAF for recessive [default: 0.01].
  --maf-comp-het=<f>          Maximum MAF for compound het [default: 0.01].
  --min-depth=<f>             Minimum depth for variant [default: 5].
  --max-cohort-ac=<f>         Maximum allele count within cohort [default: Inf].
  --max-cohort-af=<f>         Maximum allele frequency within cohort [default: Inf].
  --gtex-rpkm=<f>             Path to GTEx_median_rpkm_file.
  --omim-genemap2=<f>         Path to OMIM_genemap2_file.
"
# args <- ("S35872_1.subset.vcf.gz S35872_1 S35872_1=S35872_1.merged.bam \
#     --gene-lists AGHA-0152,AGHA-3279 \
#     --maf-dom 0.0001 \
#     --maf-rec 0.01 \
#     --maf-comp-het 0.01 \
#     --gtex-rpkm /stornext/Bioinf/data/lab_bahlo/public_datasets/GTEx/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz \
#     --omim-genemap2 /stornext/Bioinf/data/lab_bahlo/ref_db/human/OMIM/OMIM_2020-04-29/genemap2.txt \
#     --max-cohort-af 0.10") %>%
#  str_split('\\s+', simplify = T) %>%
#  str_trim()
# opts <- docopt(doc, args)

opts <- docopt(doc)
# print options
message('Using options:')
opts[names(opts) %>% 
       keep(~str_detect(., '^[:alpha:]'))] %>% 
  { class(.) <- c('list', 'docopt'); .} %>% 
  print()

# set options
sample_id <- str_extract(opts$sample_bam, '^[^=]+')
sample_bam <- str_extract(opts$sample_bam, '[^=]+$')
sampleID <- list("proband") %>% setNames(sample_id)
inheritance_MAF <- list("individual dominant"  = as.numeric(opts$`--maf-dom`),
                        "individual recessive" = as.numeric(opts$`--maf-rec`),
                        "individual comp het"  = as.numeric(opts$`--maf-comp-het`))
min_depth <- as.numeric(opts$`--min-depth`)
max_cohort_ac <- as.numeric(opts$`--max-cohort-ac`)
max_cohort_af <- as.numeric(opts$`--max-cohort-af`)
primary_panels <- c(str_split(opts$`--gene-lists`, ',', simplify = TRUE))

panelapp_tbl <- 
  read_rds('~/packages/panelapp/panel_app_table.rds') %>% 
  mutate(gene = {
    at <- which(gene %in% cavalier::HGNC_alias$alias)
    replace(gene, at, cavalier::HGNC_alias$symbol[match(gene[at], cavalier::HGNC_alias$symbol)])
    })

# find similar panels
min_sim <- 0.50
panelapp_sim <- read_rds('~/packages/panelapp/panel_app_sim.rds')
secondary_panels <-
  panelapp_sim %>% 
  filter(id %in% primary_panels,
         !sim_id %in% primary_panels,
         sim > min_sim,
         !is_subset) %>% 
  mutate(is_subset = {
    xid <- sim_id
    filter(panelapp_sim, id %in% xid, sim_id %in% xid) %>% 
      filter(is_subset) %>% 
      pull(sim_id) %>% 
      { xid %in% .}
  }) %>% 
  filter(!is_subset) %>% 
  mutate(similar_to = str_c(id, ' (', format(sim, digits = 2), ')')) %>% 
  select(panel_id = sim_id) %>% 
  distinct() %>% 
  pull(panel_id)

panels_tbl <- 
  panelapp_tbl %>% 
  filter(panel_id %in% c(primary_panels, secondary_panels)) %>% 
  nest(panel_data=(-gene))

# load variants
vars <- 
  load_vep_vcf(opts$vcf, sampleID) %>% 
  mutate(Polyphen2 = if_else(Polyphen2 == 'unknown', NA_character_, Polyphen2)) %>% 
  mutate(is_tolerated = (SIFT == 'tolerated' & Polyphen2 == 'benign') |
           (SIFT == 'tolerated' & is.na(Polyphen2)) |
           (is.na(SIFT) & Polyphen2 == 'benign'))

# filter variants
filtvars <-
  vars %>% 
  filter(IMPACT %in% c('MODERATE', 'HIGH'),
         (!is_tolerated) | is.na(is_tolerated),
         !is.na(gene),
         AC <= max_cohort_ac,
         AF <= max_cohort_af,
         ) %>% 
  as.data.frame() %>% 
  filter_variants(sampleID, inheritance_MAF, MAF_column="MAF_gnomAD", min_depth = min_depth) %>% 
  as_tibble()

# intersect with gene lists
candvars <-
  inner_join(filtvars, panels_tbl, by = "gene") %>% 
  arrange(gene, position)

# create cavalier output if any variants remain
dir.create(opts$out, recursive = TRUE, showWarnings = FALSE)
if (nrow(candvars)) {
  output_cols <- c('Inheritance', "Variant", "Amino acid", "change", "Depth (R,A)", "Cohort AF",
                   "GnomAD AF", "SIFT", "Polyphen2", "Grantham", "RVIS", "GeVIR")
  candvars %>%
  create_igv_snapshots(sample_bam, "hg19", vcfs = opts$vcf) %>%
  #   create_igv_snapshots(sample_bam,  "hg19",
  #                        vcfs = opts$vcf,
  #                        # overwrite = TRUE,
  #                        singularity_img = '~/links/singularity_cache/jemunro-cavalier-dev.img',
  #                        singularity_bin = '/stornext/System/data/apps/singularity/singularity-3.7.3/bin/singularity') %>%
  #   mutate(panel_data = map(seq_along(gene), ~tibble())) %>% 
    mutate(Inheritance = `inheritance model`,
           Variant = str_c(chromosome, ':', position, ':', reference, '>', alternate),
           `Depth (R,A)` = `proband depth (R,A)`,
           `Amino acid` = str_replace(Amino_acids, '/', '>'),
           Polyphen2 = str_c(Polyphen2, ' (', Polyphen2_score, ')'),
           SIFT = str_c(SIFT, ' (', SIFT_score, ')'),
           title = str_c('Sample: ', sample_id, ', Gene: ', gene),
           `Cohort AF` = str_c(round(AF, 4), ' (', AC, ')')
    ) %>%
    rename(`GnomAD AF` = MAF_gnomAD) %>% 
    as.data.frame() %>%
    create_cavalier_output(opts$out, sampleID, output_cols,
                           hide_missing_igv = TRUE,
                           layout = "individual",
                           genemap2 = opts$`--omim-genemap2`,
                           GTEx_median_rpkm = opts$`--gtex-rpkm`,
                           title_col = 'title', 
                           add_data_col = 'panel_data',
                           output = c('html', 'ppt'))
}

