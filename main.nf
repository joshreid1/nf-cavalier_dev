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
// ref params
params.ref_fasta = ''
// params.ref_hg38 = true
params.vep_cache = ''
params.vep_cache_ver = '114'
// params.vep_assembly = params.ref_hg38 ? 'GRCh38' : 'GRCh37'
params.ref_gene = ''

// SNV args
params.snv_vcf = ''
params.snv_caller = 'GATK'
params.snv_n_shards = 200
params.snv_info = ['AC', 'AF', 'AN']
params.snv_format = ['GT', 'AD'] // TODO: ensure 'GT' is always present
// optional
params.snv_vcfanno = [
    [   vcf: '/vast/projects/bahlo_cache/annotation/gnomAD/joint_sites_4.1.vcf.gz', 
        csi: true, 
        fields: [gnomad_AF: 'AF', gnomad_AC: 'AC', gnomad_FAF95: 'fafmax_faf95_max', gnomad_nhomalt: 'nhomalt']
    ],
    [   tsv: '/vast/projects/bahlo_cache/annotation/CADD/1.7/whole_genome_SNVs.tsv.gz',
        fields: [CADD: 6]
    ],
    [   tsv: '/vast/projects/bahlo_cache/annotation/CADD/1.7/gnomad.genomes.r4.0.indel.tsv.gz',
        fields: [CADD: 6]
    ],
    [   tsv: '/vast/projects/bahlo_cache/annotation/phyloP/hg38.phyloP100way.bed.gz',
        fields: [phyloP100: 4],
        op: 'max'
    ]
]
params.snv_vcfanno_filter = 'gnomad_AF<0.01 || gnomad_AF="."'

params.vep_spliceai_snv   = '/vast/projects/bahlo_cache/annotation/spliceAI/1.3/spliceai_scores.raw.snv.hg38.vcf.gz'
params.vep_spliceai_indel = '/vast/projects/bahlo_cache/annotation/spliceAI/1.3/spliceai_scores.raw.indel.hg38.vcf.gz'
params.vep_alphamissense  = '/vast/projects/bahlo_cache/annotation/alphamissense/AlphaMissense_hg38.tsv.gz'
params.vep_revel          = '/vast/projects/bahlo_cache/annotation/REVEL/revel_1.3.hg38.vep.tsv.gz'
params.vep_utr_annotator  = '/vast/projects/bahlo_cache/annotation/UTRannotator/uORF_5UTR_GRCh38_PUBLIC.txt'

// Can append multiple, by default we always source "$projectDir/misc/scripts/snv_functions.R"
params.report_func_source  = null
/* functions run in order listed
    - first option has file path as first arg
    - remainder have a dataframe as first arg
    - second arg for all in config (see params.snv_config)
*/
// must pe defined in params.report_func_source or "$projectDir/misc/scripts/snv_functions.R"
params.snv_report_functions = 'SNV_LOAD,SNV_FILTER_DEPTH,SNV_FILTER_GENES,SNV_FILTER_TYPE,SNV_FILTER_INHERITANCE'
// all params.snv_report_* fields are collected and passed to snv_report_functions
params.snv_report_min_depth = 5
params.snv_report_freq_thresholds = [
    pop   : [AF: [dominant: 0.0001, recessive: 0.01 ], AC: [dominant: null , recessive: null ]],
    cohort: [AF: [dominant: null  , recessive: null],  AC: [dominant: null , recessive: null ]]
]
// fields to report for all variants
params.snv_report_fields =  [
    Gene         : 'Gene',
    Inheritance  : 'Inheritance', 
    Consequence  : 'Consequence',
    'HGVS coding': 'HGVSc',
    'HGVS prot.' : 'HGVSp',
    ClinVar      : 'CLIN_SIG', 
    gnomAD       : 'gnomAD', 
    Cohort       : 'cohort',
    CADD         : 'CADD',
    PhyloP100    : 'phyloP100',
    SIFT         : 'SIFT',
    PolyPhen     : 'PolyPhen',
    AlphaMissense: 'AlphaMissense',
    REVEL        : 'REVEL'
]
// additional options for cavalier R pacakge
params.cavalier_options = [
    database_mode: 'fallback',
    cache_dir: 'cavalier_cache'
]
// SV params
params.pop_sv = ''
params.sv_n_shards  = 20
params.include_sv_csv = true
params.sv_types = 'DEL,DUP,INS,INV,BND'
params.sv_type_match = [DEL: ['DEL'], DUP: ['CNV', 'DUP']]
params.no_slides = true // skip making slides in cavalier

include { validate_params } from './functions/validate'

validate_params()

include { CAVALIER } from './workflows/cavalier'

workflow {
    CAVALIER()
}
