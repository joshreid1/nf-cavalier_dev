#!/usr/bin/env nextflow
/*
TODO:
    - Add download URLs for reference files as default to simplify configuration
 */
//input params
params.outdir = 'output'

params.fill_tags = true
params.remove_fields = 'INFO/CSQ'
params.sv_vcf = ''
params.ped = ''
params.bams = ''
params.lists = ''

// SNP args
params.snp_vcf = ''
params.snp_caller = 'GATK'
params.snp_n_shards = 200
params.snp_format_keep = 'FORMAT/GT,FORMAT/AD'
params.snp_info_keep = 'INFO/AC,INFO/AF,INFO/AN'
// optional
params.snp_vcfanno = [
    [   vcf: '/vast/projects/bahlo_cache/annotation/gnomAD/joint_sites_4.1.vcf.gz', 
        index: 'csi', 
        fields: [gnomad_AF: 'AF', gnomad_AC: 'AC', gnomad_FAF95: 'fafmax_faf95_max', gnomad_nhomalt: 'nhomalt']
    ],
    [   tsv: '/vast/projects/munro_data/ref-data/whole_genome_SNVs.tsv.gz',
        index: 'tbi',
        fields: [CADD: 6]
    ],
    [   tsv: '/vast/projects/bahlo_cache/annotation/CADD/1.7/gnomad.genomes.r4.0.indel.tsv.gz',
        index: 'tbi',
        fields: [CADD: 6]
    ],
    [   tsv: '/vast/projects/bahlo_cache/annotation/phyloP/hg38.phyloP100way.bed.gz',
        index: 'tbi',
        fields: [phyloP100: 4],
        op: 'mean'
    ]
]

params.snp_vcfanno_filter = 'gnomad_AF<0.01 || gnomad_AF="."'
params.vep_spliceai_snv   = '/vast/projects/munro_data/ref-data/spliceAI/1.3/spliceai_scores.raw.snv.hg38.vcf.gz'
params.vep_spliceai_indel = '/vast/projects/munro_data/ref-data/spliceAI/1.3/spliceai_scores.raw.indel.hg38.vcf.gz'
params.vep_alphamissense  = '/vast/projects/bahlo_cache/annotation/alphamissense/AlphaMissense_hg38.tsv.gz'
params.vep_revel          = '/vast/projects/bahlo_cache/annotation/REVEL/revel_1.3.hg38.vep.tsv.gz'
// params.vep_dbnsfp         = '/vast/projects/bahlo_cache/dbNSFP/dbNSFP5.1a_grch38.gz'
// params.vep_dbnsfp_fileds  = 'REVEL_score,AlphaMissense_score,AlphaMissense_pred,MetaSVM_score,MetaSVM_pred'


params.sv_n_shards  = 20

// filter params
params.maf_dom = 0.0001
params.maf_de_novo = 0.0001
params.maf_rec = 0.01
params.maf_comp_het = 0.01
params.max_cohort_af = 1.0
params.max_cohort_ac = 'Inf'
params.min_impact = 'MODERATE'
params.exclude_benign_missense = true
params.include_sv_csv = true
params.variants_override = null // disable filtering and report specific set of variants by id in TSV file

// ref params
params.ref_fasta = ''
// params.ref_hg38 = true
params.vep_cache = ''
params.vep_cache_ver = '114'
// params.vep_assembly = params.ref_hg38 ? 'GRCh38' : 'GRCh37'
params.pop_sv = ''
params.ref_gene = ''

// map additional options for cavalier R pacakge
params.cache_dir = 'cavalier_cache'
params.database_mode = 'fallback' // one of 'latest', 'fallback' or 'offline'
params.cavalier_options = [:]

// exec params
params.sv_types = 'DEL,DUP,INS,INV,BND'
params.sv_type_match = [DEL: ['DEL'], DUP: ['CNV', 'DUP']]
params.no_slides = false // skip making slides in cavalier


include { validate_params } from './functions/validate'

validate_params()

include { CAVALIER } from './workflows/cavalier'

workflow {
    CAVALIER()
}
