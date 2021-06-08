#!/usr/bin/env Rscript

stopifnot(require(tidyverse),
          require(docopt),
          require(cavalier))

doc <- "
Usage:
  cavalier_singleton.R <vcf> <out> <sample_bam> [options]

Options:
  vcf                         Input VEP vcf file.
  out                         Output directory.
  sample_bam                  Sample name and bam file in format name=/path/to/bam.
  --genome=<f>                Reference genome for IGV snapshot [default: hg19].
  --gene-lists=<f>            Comma serparated list of gene list files.
  --remainder                 When using gene lists, output variants not on any list as well.
  --maf-dom=<f>               Maximum MAF for dominant [default: 0.0001].
  --maf-rec=<f>               Maximum MAF for recessive [default: 0.01].
  --maf-comp-het=<f>          Maximum MAF for compound het [default: 0.01].
  --gtex-rpkm=<f>             Path to GTEx_median_rpkm_file.
  --omim-genemap2=<f>         Path to OMIM_genemap2_file.
"
# args <- ("S35167_3.subset.vcf.gz S35167_3 S35167_3=S35167_3.merged.bam \
#     --gene-lists differences_of_sex_development.txt,primary_ovarian_insufficiency_premature_ovarian_failure.txt \
#     --maf-dom 0.0001 \
#     --maf-rec 0.01 \
#     --maf-comp-het 0.01 \
#     --gtex-rpkm /stornext/Bioinf/data/lab_bahlo/public_datasets/GTEx/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz \
#     --omim-genemap2 /stornext/Bioinf/data/lab_bahlo/ref_db/human/OMIM/OMIM_2020-04-29/genemap2.txt") %>%
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

sample_id <- str_extract(opts$sample_bam, '^[^=]+')
sample_bam <- str_extract(opts$sample_bam, '[^=]+$')
sampleID <- list("proband") %>% setNames(sample_id)
inheritance_MAF <- list("individual dominant"  = as.numeric(opts$maf_dom),
                        "individual recessive" = as.numeric(opts$maf_rec),
                        "individual comp het"  = as.numeric(opts$maf_comp_het))
dir.create(opts$out)
gene_lists <- 'AGHA-0289'

vars <- 
  load_vep_vcf(opts$vcf, sampleID) %>%
  mutate(Polyphen2 = if_else(Polyphen2 == 'unknown', NA_character_, Polyphen2)) %>% 
  mutate(is_tolerated = (SIFT == 'tolerated' & Polyphen2 == 'benign') |
           (SIFT == 'tolerated' & is.na(Polyphen2)) |
           (is.na(SIFT) & Polyphen2 == 'benign'))

filtvars <-
  vars %>% 
  filter(!is_tolerated,
         !is.na(gene)) %>% 
  as.data.frame() %>% 
  filter_variants(sampleID, inheritance_MAF, MAF_column="MAF_gnomAD") %>% 
  as_tibble()

# add gene lists to variants
if (!is.null(opts$gene_lists)) {
  gene_lists <-
    tibble(fn = c(str_split(opts$gene_lists, ',', simplify = TRUE))) %>% 
    mutate(gene_list = basename(fn) %>% str_remove('\\.txt'),
           gene = map(fn, scan, what = character())) %>% 
    unnest(gene) %>% 
    select(gene_list, gene)
  
  if (opts$remainder) {
    filtvars <- 
      left_join(filtvars, gene_lists, by = "gene") %>% 
      mutate(gene_list  = replace_na(gene_list, 'remainder'))
  } else {
    filtvars <- inner_join(filtvars, gene_lists, by = "gene")
  }
} else {
  filtvars <- mutate(filtvars, gene_list = 'all')
}

# create cavalier output if any variants remain
if (nrow(filtvars)) {
  output_cols <- c("gene_list", 'inheritance model', "variant", "amino_acid", "change",
                   "gene", "MAF_gnomAD", "SIFT", "Polyphen2", "Grantham", "RVIS")
  create_igv_snapshots(filtvars, sample_bam, "hg19", 'igv') %>%
    mutate(sample_id = sample_id,
           variant = str_c(chromosome, ':', position, ':', reference, '>', alternate),
           amino_acid = str_replace(Amino_acids, '/', '>'),
           Polyphen2 = str_c(Polyphen2, ' (', Polyphen2_score, ')'),
           SIFT = str_c(SIFT, ' (', SIFT_score, ')'),
           title = str_c('sample: ', sample_id),
    ) %>%
    as.data.frame() %>%
    split.data.frame(.$gene_list) %>%
    walk2(names(.), ., function(list_name, list_vars) {
      create_cavalier_output(list_vars, file.path(opts$out, list_name), sampleID, output_cols,
                             hide_missing_igv = TRUE, layout = "individual",
                             genemap2 = opts$omim_genemap2,
                             GTEx_median_rpkm = opts$gtex_rpkm,
                             title_col = 'title'
      )
    })
}


