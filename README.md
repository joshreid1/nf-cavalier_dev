# nf-cavalier

Nextflow Pipeline for singleton and family based candidate variant prioritisation based on gene lists using the Cavalier R package. This pipeline is a work in progress.

## Installation
* Clone this repositoty

## Usage
* Create and navigate to run working directory
* Create configuration file in run directory named `nextflow.config`:
  ```Nextflow
    params {
      // inputs
      id = 'my_cohort'
      vcf = 'my_cohort.SNPs.vcf.gz'
      sv_vcf = 'my_cohort.SVs.vcf.gz'
      ped = 'families.ped'
      bams = 'bams.tsv'
      lists = 'gene_lists.tsv'
      
      // filtering
      maf_dom = 0.0001
      maf_rec = 0.01
      maf_comp_het = 0.01
      maf_de_novo = 0.0001
      max_cohort_af = 1.0
      max_cohort_ac = 'Inf'
      min_impact = 'MODERATE'
      exclude_benign_missense = true
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
* **Params**  
  * `id` - Unique name for run. Used to name output files.
  * `vcf` - (Optional) Input VCF file with SNP variant calls for all samples.
  * `sv_vcf` - (Optional) Input VCF file with SV variant calls for all samples.
  * `ped` - A [Ped format file](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format) describing familial relationships, with 1/2 coding for unaffected/affected phenotypes (missing phenotype not supported).
  * `bams` - TSV file with first column containing individual ID, second column containing path to indexed BAM file (no header row/  column names).
  * `lists` - TSV file with first column containing family ID (matching first column of PED file), second column containing the path to a gene list file (no header row/  column names) or an gene list identifier (e.g. [PAA:289](https://panelapp.agha.umccr.org/panels/289/), [PAE:20](https://panelapp.genomicsengland.co.uk/panels/20/) or [HP:0001251](https://hpo.jax.org/app/browse/term/HP:0001251)). e.g.:
    ```
    FAM-01  /path/to/genes.tsv
    FAM-02  PAA:289
    FAM-02  HP:0001251
    ```
    * Each gene list file should itself be a TSV file with mandatory named columns "list_id", "list_name" and "gene". The "gene" column should contain HGNC symbol for each gene. Additional metadata columns may also be included such as "version", e.g.:
        ```
        list_id list_name       version gene    inheritance     status
        AGHA#289        Ataxia_Superpanel       0.437   AAAS    BIALLELIC       GREEN
        AGHA#289        Ataxia_Superpanel       0.437   ABCB7   X-LINKED        GREEN
        AGHA#289        Ataxia_Superpanel       0.437   ABCD1   X-LINKED        GREEN
        ```
    * `exclude_benign_missense` - exclude missense variants that are predicted/annotated as benign by all of Sift, 
    Polyphen and ClinVar. Missing annotations are ignored.
    * `include_sv_csv` - Include "coding_sequence_variant" SVs regardless of VEP Impact.
* First run:  
`nextflow run /PATH/TO/nf-cavalier`
* Resume run:  
`nextflow run /PATH/TO/nf-cavalier -resume`
* Recommended to run workflow in a `screen` session
  
