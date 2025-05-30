#!/usr/bin/env nextflow

// main input params
params.outdir = 'output'
params.ped = ''
params.bams = ''
params.lists = ''
params.snv_vcf = ''
params.snv_vcf_annotated = null // skip annotation by providing pre-annotated VCF (output from another run)
params.sv_vcf = '' // not implemented

// ref fasta
params.ref_fasta = '/vast/projects/bahlo_cache/ref_genome/hg38_GATK/Homo_sapiens_assembly38.fasta'

// options for cavalier R pacakge
params.cavalier_options = [
    database_mode: 'fallback', // will try to get latest db/vers, but fallback to cache if unavailable (set to 'offline' no check)
    cache_dir: '.cavalier_cache', // will write any downloaded files to this dir
    read_only_cache_dir: '/vast/projects/bahlo_cache/cavalier_cache' // will use cache files from
]

/* =================== SNV/Indel ARGS =================== */
// split input VCF into shards for parallel processing
params.snv_n_shards = 200
// apply filter to input variants
params.snv_vcf_filter = "PASS,."
// INFO fields to keep from VCF
params.snv_info = ['AC', 'AF', 'AN']
// FORMAT fields to keep from VCF
params.snv_format = ['GT', 'GQ', 'DP'] // TODO: ensure 'GT' is always present
// fill AC, AF, and AN from VCF - recommended to set true if VCF does not have these set
params.snv_fill_tags = false
// vcfanno settings
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
// filter to apply after VCF anno, drastically reduces VEP runtime if we drop common variants
params.snv_vcfanno_filter = 'gnomad_AF<0.01 || gnomad_AF="."' // (keep if AF < 0.01 or AF is missing)
// VEP settings
params.vep_cache     = '/vast/projects/bahlo_cache/vep_cache'
params.vep_cache_ver = '114'
// check number of variants output by VEP equal to number input
params.vep_check     = true
// VEP plugins (setting to null will disable plugin)
params.vep_spliceai_snv   = '/vast/projects/bahlo_cache/annotation/spliceAI/1.3/spliceai_scores.raw.snv.hg38.vcf.gz'
params.vep_spliceai_indel = '/vast/projects/bahlo_cache/annotation/spliceAI/1.3/spliceai_scores.raw.indel.hg38.vcf.gz'
params.vep_alphamissense  = '/vast/projects/bahlo_cache/annotation/alphamissense/AlphaMissense_hg38.tsv.gz'
params.vep_revel          = '/vast/projects/bahlo_cache/annotation/REVEL/revel_1.3.hg38.vep.tsv.gz'
params.vep_utr_annotator  = '/vast/projects/bahlo_cache/annotation/UTRannotator/uORF_5UTR_GRCh38_PUBLIC.txt'

// all params.snv_report_* fields are collected and passed to snv_report_functions
// Can append multiple, by default we always source "$projectDir/misc/scripts/snv_functions.R"
params.report_func_source  = null
/* functions run in order listed
    - first option has file path as first arg
    - remainder have a dataframe as first arg
    - second arg for all in config (see params.snv_config)
*/
// functions must be defined in "./bin/snv_functions.R" or params.report_func_source
params.snv_report_functions = 'SNV_LOAD,SNV_FILTER_FMT,SNV_FILTER_GENES,SNV_FILTER_TYPE,SNV_FILTER_INHERITANCE,SNV_REPORT'
params.snv_report_min_DP = 5  // at least depth of 5 in all samples
params.snv_report_min_GQ = 10 // at least 90% confidence of sample geneotypes
// used by SNV_FILTER_INHERITANCE function, strictly less than
params.snv_report_freq_thresholds = [
    // adjust these based on expected frequencies of diease-of-interest
    pop: [
        AF : [dominant: 0.0001, recessive: 0.01 ], 
        AC : [dominant: 5     , recessive: 50   ],
        hom: [dominant: 2     , recessive: 10   ]
    ],
    cohort: [
        AF : [dominant: 0.01  , recessive: 0.05  ],
        AC : [dominant: 2     , recessive: 5     ]
    ]
]
// fields to report for all variants
params.snv_report_fields =  [
    Gene         : 'Gene',
    Inheritance  : 'Inheritance', 
    Consequence  : 'Consequence',
    HGVS         : 'HGVS',
    ClinVar      : 'CLIN_SIG',
    'gnomAD v4.1': 'gnomAD', 
    Cohort       : 'Cohort',
    PhyloP100    : 'phyloP100',
    CADD         : 'CADD',
    REVEL        : 'REVEL',
    AlphaMissense: 'AlphaMissense',
    SIFT         : 'SIFT',
    PolyPhen     : 'PolyPhen',
    SpliceAI     : 'SpliceAI'
]

/* =================== SV ARGS =================== */
// ( not currently implemented )
params.pop_sv = ''
params.ref_gene = ''
params.sv_n_shards  = 20
params.include_sv_csv = true
params.sv_types = 'DEL,DUP,INS,INV,BND'
params.sv_type_match = [DEL: ['DEL'], DUP: ['CNV', 'DUP']]
params.no_slides = false // skip making slides in cavalier

include { validate_params } from './functions/validate'

validate_params()

include { CAVALIER } from './workflows/cavalier'

workflow {
    CAVALIER()
}
