# nf-cavalier

Nextflow Pipeline for singleton and family based candidate variant prioritisation based on gene lists using the Cavalier R package. This pipeline is a work in progress.

## Usage

* **Example `nextflow.config`**
  ```Nextflow
    params {
      // inputs
      id = 'UDP-hg38'
      vcf = 'UDP-hg38.pass.vcf.gz'
      ped = 'families.ped'
      bams = 'bams.tsv'
      lists = 'gene_lists.tsv'
      
      // filtering
      maf_dom = 0.0001
      maf_rec = 0.01
      maf_comp_het = 0.01
      maf_de_novo = 0.0001
      max_cohort_af = 0.10
      min_impact = 'MODERATE'
      
      // reference config
      ref_hg38 = true
      ref_fasta = '/stornext/Bioinf/data/lab_bahlo/ref_db/human/hg38/GATK/fasta_no_alt/hg38.no_alt.fasta'
      omim_genemap2 = '/stornext/Bioinf/data/lab_bahlo/ref_db/human/OMIM/OMIM_2021-08-17/genemap2.txt'
      vep_cache = '/stornext/Bioinf/data/lab_bahlo/ref_db/vep-cache'
      vep_cache_ver = '104'
    }
    ```
* **Params**  
  * `id` - Unique name for run. Used to name output files.
  * `vcf` - Input VCF file with variant calls for all samples
  * `ped` - A [Ped format file](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format) describing familial relationships, with 1/2 encoding for unaffected/affected phenotypes (missing phenotype no supported)
  * `bams` - TSV file with first column containing individual ID, second column containing path to indexed BAM file (no header row/  column names)
  * `lists` - TSV file with first column containing family ID (matching first column of PED file), second column containing the path to a gene list file (no header row/  column names). Each gene list file should itself be a TSV file with mandatory named columns "list_id", "list_name" and "gene". The "gene" column should contain HGNC symbol for each gene. Additional metadata columns may also be included such as "version", e.g.:
    ```
    list_id list_name       version gene    inheritance     status
    AGHA#289        Ataxia_Superpanel       0.437   AAAS    BIALLELIC       GREEN
    AGHA#289        Ataxia_Superpanel       0.437   ABCB7   X-LINKED        GREEN
    AGHA#289        Ataxia_Superpanel       0.437   ABCD1   X-LINKED        GREEN
    ```
  
