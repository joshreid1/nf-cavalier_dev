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

opts <- docopt(doc)
# opts <- 
#   docopt(doc, c('/stornext/HPCScratch/home/munro.j/runs/udp/cav/output/vcf_family_subset/S33843_1.subset.vcf.gz',
#                 'out', 
#                 'S33843_1=/bam/file',
#                 '--gene-lists', '~/analyses/udp/gene_lists/ataxia_superpanel.txt',
#                 '--gtex-rpkm', '/stornext/Bioinf/data/lab_bahlo/public_datasets/GTEx/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz',
#                 '--omim-genemap2', '/stornext/Bioinf/data/lab_bahlo/ref_db/human/OMIM/OMIM_2019-05-04/genemap2.txt'))
message('Using options:')
print(opts)

sample_id = str_extract(opts$sample_bam, '^[^=]+')
sample_bam = str_extract(opts$sample_bam, '[^=]+$')
sampleID <- list("proband") %>% setNames(sample_id)
inheritance_MAF <- list("individual dominant"  = as.numeric(opts$maf_dom),
                        "individual recessive" = as.numeric(opts$maf_rec),
                        "individual comp het"  = as.numeric(opts$maf_comp_het))
output_dir <- opts$out

vars <- 
  load_vep_vcf_2(opts$vcf, sampleID) %>% 
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

if (opts$gene_list) {
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

filtvars <- create_igv_snapshots(filtvars, sample_bam, "hg19", 'igv')

output_cols <- c("chromosome", "position", "reference", "alternate", "gene", "region", "change", "MAF_gnomAD", "SIFT", "Polyphen2", "Grantham", "RVIS")

filtvars %>% 
  split.data.frame(.$gene_list) %>% 
  walk2(names(.), ., function(list_name, list_vars) {
    create_cavalier_output(list_vars, file.path(opts$out, list_name), sampleID, output_cols,
                           hide_missing_igv = TRUE, layout = "individual", 
                           genemap2 = opts$omim_genemap2, 
                           GTEx_median_rpkm = opts$gtex_rpkm) 
  })

  