# nf-cavalier

Nextflow Pipeline for singleton and family based candidate variant prioritisation based on gene lists using the Cavalier R package. This pipeline is a work in progress.

## Installation
* Clone this repositoty

## Usage
* Create and navigate to run working directory
* Create configuration file in run directory named `nextflow.config`:
  ```Nextflow
    params {
      // output directory
      outdir = 'output'
      
      // inputs
      snp_vcf = 'my_cohort.SNPs.vcf.gz'
      sv_vcf = 'my_cohort.SVs.vcf.gz'
      ped = 'families.ped'
      bams = 'bams.tsv'
      lists = 'my_list_file.tsv,HP:0001250'
      
      // filtering
      maf_dom = 0.0001
      maf_rec = 0.01
      maf_comp_het = 0.01
      maf_de_novo = 0.0001
      max_cohort_af = 1.0
      max_cohort_ac = 'Inf'
      min_impact = 'LOW'
      exclude_benign_missense = false
      include_sv_csv = true
  
      // reference config
      ref_hg38 = true
      ref_fasta = '/PATH/TO/GRCh38.fasta'
      pop_sv = '/PATH/TO/gnomad-sv.vcf.gz'
      ref_gene = '/PATH/TO/RefSeqGene.hg38.UCSC.txt'
      vep_cache = '/PATH/T0/vep-cache'
      vep_cache_ver = '104'
    }
    ```
* Note: Bahlo Lab members should use this [config](https://github.com/bahlolab/nextflow-config/blob/master/nf-cavalier/milton.config) as a starting point.

* **Params**  
  * `outdir` - Output directory
  * `vcf` - (Optional) Input VCF file with SNP variant calls for all samples.
  * `sv_vcf` - (Optional) Input VCF file with SV variant calls for all samples.
  * `ped` - A [Ped format file](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format) describing familial relationships, with 1/2 coding for unaffected/affected phenotypes (missing phenotype not supported).
  * `bams` - TSV file with first column containing individual ID, second column containing path to indexed BAM file (no header row/  column names).
  * `lists` - Comma separated list of gene lists to use for filtering. This may be a local TSV file, e.g. 'my_list_file.tsv' or a web based gene list, e.g. [PAA:289](https://panelapp.agha.umccr.org/panels/289/).
    * **Local TSV file**  
      * Path to a TSV file with mandatory named column. The file should have at least one of the following column names:
     `ensembl_gene_id`, `hgnc_id`, `entrez_id` or `symbol`. Optional column names are `list_id`, `list_name`, `list_version` and `inheritance`. Note that all gene IDs are ultimately converted to Ensembl Gene IDs using HGNC to match VEP annotation.
     
         e.g.

          ```
          list_id	list_name	list_version	symbol	inheritance
          PAA:202	Genetic Epilepsy	1.26	AARS1	AR
          PAA:202	Genetic Epilepsy	1.26	ABAT	AR
          PAA:202	Genetic Epilepsy	1.26	ABCA2	AR
          ```
          
          or:

            ```
            list_id	list_name	list_version	ensembl_gene_id	inheritance
            PAA:202	Genetic Epilepsy	1.26	ENSG00000090861	AR
            PAA:202	Genetic Epilepsy	1.26	ENSG00000183044	AR
            PAA:202	Genetic Epilepsy	1.26	ENSG00000107331	AR
            ```  
    * **Web List** - Cavalier will automatically retrieve the latest version of these web lists
      * **PanelApp**: PanelApp Australia or PanelApp Genomics England lists may be specified with "PAA:" or "PAE:" prefix respectively. e.g. [PAA:289](https://panelapp.agha.umccr.org/panels/289/)
      * **HPO**: Human phenotype ontology terms may be specified with the "HP:" prefex, e.g. [HP:0001250](https://hpo.jax.org/browse/term/HP:0001250)
      * **Genes4Epilepsy**: [Genes4Epilepsy](https://github.com/bahlolab/Genes4Epilepsy) lists may be specified with the "G4E:" prexif, e.g. "G4E:All" for All Epilepsy genes, or "G4E:Focal" for Focal epilepsy genes only
      * **HGNC**: Gene subsets by locus group can be extracted from HGNC, for example "HGNC:protein-coding" will give a list of all protein coding genes
  * `exclude_benign_missense` - exclude missense variants that are predicted/annotated as benign by all of Sift, 
  Polyphen and ClinVar. Missing annotations are ignored.
  * `include_sv_csv` - Include "coding_sequence_variant" SVs regardless of VEP Impact.

* First run:  
`nextflow run /PATH/TO/nf-cavalier`
* Resume run:  
`nextflow run /PATH/TO/nf-cavalier -resume`
* Recommended to run workflow in a `screen` session
  
