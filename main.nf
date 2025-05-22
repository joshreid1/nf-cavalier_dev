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
    [   tsv: '/vast/projects/munro_data/ref-data/whole_genome_SNVs.tsv.gz',
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

params.vep_spliceai_snv   = '/vast/projects/munro_data/ref-data/spliceAI/1.3/spliceai_scores.raw.snv.hg38.vcf.gz'
params.vep_spliceai_indel = '/vast/projects/munro_data/ref-data/spliceAI/1.3/spliceai_scores.raw.indel.hg38.vcf.gz'
params.vep_alphamissense  = '/vast/projects/bahlo_cache/annotation/alphamissense/AlphaMissense_hg38.tsv.gz'
params.vep_revel          = '/vast/projects/bahlo_cache/annotation/REVEL/revel_1.3.hg38.vep.tsv.gz'
params.vep_utr_annotator  = null
// params.vep_utr_annotator  = '/vast/projects/bahlo_cache/annotation/UTRannotator/uORF_5UTR_GRCh38_PUBLIC.txt'


// params.vep_dbnsfp         = '/vast/projects/bahlo_cache/dbNSFP/dbNSFP5.1a_grch38.gz'
// params.vep_dbnsfp_fileds  = 'REVEL_score,AlphaMissense_score,AlphaMissense_pred,MetaSVM_score,MetaSVM_pred'

// NB: filter_AF = pmax(replace_na(gnomad_FAF95, 0), replace_na(gnomad_AF, 0))

// filter params

// default missing values to 'Inf'
params.freq_filters = [
    snv: [
        pop   : [AF: [dominant: 0.0001, recessive: 0.01 ], AC: [dominant: 'Inf' , recessive: 'Inf']],
        cohort: [AF: [dominant: 'Inf' , recessive: 'Inf'], AC: [dominant: 'Inf' , recessive: 'Inf']]
    ],
    sv : [
        pop   : [AF: [dominant: 0.0001, recessive: 0.01 ], AC: [dominant: 'Inf' , recessive: 'Inf']],
        cohort: [AF: [dominant: 'Inf' , recessive: 'Inf'], AC: [dominant: 'Inf' , recessive: 'Inf']]
    ]
]

params.snv_mutate = [
    pop_AF: 'pmax(replace_na(gnomad_AF, 0), replace_na(gnomad_FAF95, 0))',
    pop_AC: 'replace_na(gnomad_AC, 0)',
    gnomAD: 'str_c("AF=", gnomad_AF, "; AC=", gnomad_AC, "; Hom=",  gnomad_nhomalt)',
    SpliceAI_max: 'pmax(SpliceAI_pred_DS_AG, SpliceAI_pred_DS_AL, SpliceAI_pred_DS_DG, SpliceAI_pred_DS_DL)',
    SpliceAI: 'str_c("AG=", SpliceAI_pred_DS_AG, "; DG=", SpliceAI_pred_DS_DG, "; AL=", SpliceAI_pred_DS_AL, "; DL=", SpliceAI_pred_DS_DL)'  
]

params.snv_report = ['gnomAD', 'Consequence', 'CADD']
params.snv_subsets = [
    splicing: [
        or: [ 
            'Consequence %in% c("splice_acceptor_variant", "splice_donor_variant")',
            'SpliceAI_max_score > 0.5'
        ],
        report: ['PhyloP100', 'SpliceAI_*']
    ],
    missense: [
        or: ['Consequence == "missense_variant"'],
        report: ['REVEL', 'am_*', 'SIFT', 'PolyPhen']
    ],
    lof: [
        or: ['IMPACT == "HIGH"'],
        report: ['PhyloP100', 'SpliceAI_*']
    ],
    generic: [
        or: [
            'IMPACT %in% c("MODERATE", "HIGH")',
            'CADD >= 20',
        ],
        report: ['PhyloP100', 'SpliceAI_*']
    ]
]

//TODO: How flexible should this be?
//  - whats the minimum level of code to support the most analyses?
//  - subsets is pretty complicated, maybe these should be predefined
//  - easier to assume run with specified vcfanno and specified vep, than be flexible
//  - maintain flexibility for devs but not for users


// TODO:
// - framework for interpretting various predictor scores described here
// - https://spliceailookup.broadinstitute.org/#variant=NM_001089.3(ABCA3)%3Ac.875A%3ET%20(p.Glu292Val)&hg=38&bc=basic&distance=500&mask=0&ra=0


// params.snv_filter = [
//     // entries combined with and, unless inside an 'or' block
//     or: [
//         'IMPACT >= MODERATE',
//         'MAX_CLIN_SIG >= "pathogenic"',
//         'CADD >= 20',
//         'am_pathogenicity >= 0.564',
//         'REVEL >= 0.7',
//         'SpliceAI_max_score >= 0.5'
//     ]
// ]

params.sv_n_shards  = 20





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
