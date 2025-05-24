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
// params.vep_utr_annotator  = null


// TODO: default missing values to 'Inf'
params.snv_freq_filters = [
    pop   : [AF: [dominant: 0.0001, recessive: 0.01 ], AC: [dominant: null , recessive: null ]],
    cohort: [AF: [dominant: null,   recessive: null],  AC: [dominant: null , recessive: null ]]
]

params.snv_min_depth = 5

params.snv_mutate = [
    title:  'str_c(family_id, SYMBOL, HGVSg, sep = " | ")',
    pop_AF: 'pmax(replace_na(gnomad_AF, 0), replace_na(gnomad_FAF95, 0))',
    pop_AC: 'replace_na(gnomad_AC, 0)',
    gnomAD: 'str_c("AF=", signif(replace_na(gnomad_AF, 0), 2), "; AC=", replace_na(gnomad_AC, 0), "; Hom=",  replace_na(gnomad_nhomalt, 0))',
    cohort: 'str_c("AF=", signif(replace_na(AF, 0), 2), "; AC=", replace_na(AC, 0))',
    phyloP100: 'replace(phyloP100, VARIANT_CLASS == "insertion", NA)',
    SpliceAI_max: 'replace_na(pmax(SpliceAI_pred_DS_AG, SpliceAI_pred_DS_AL, SpliceAI_pred_DS_DG, SpliceAI_pred_DS_DL), 0)',
    SpliceAI: 'str_c("AG=", SpliceAI_pred_DS_AG, "; DG=", SpliceAI_pred_DS_DG, "; AL=", SpliceAI_pred_DS_AL, "; DL=", SpliceAI_pred_DS_DL)',
    AlphaMissense: 'str_c(am_class, "(", am_pathogenicity, ")")',
    HGVSc: 'str_replace(HGVSc, ":", ":\\n")',
    HGVSp: 'str_replace(HGVSp, ":", ":\\n")'
]

params.snv_report = [
    Gene         : 'Gene',
    Inheritance  : 'inheritance', 
    Consequence  : 'Consequence',
    'HGVS coding': 'HGVSc',
    'HGVS prot.' : 'HGVSp',
    'ClinVar'    : 'CLIN_SIG', 
    'gnomAD'     : 'gnomAD', 
    'Cohort'     : 'cohort',
    'CADD'       : 'CADD',
    PhyloP100    : 'phyloP100'
]

params.snv_subsets = [
    lof: [
        filter:
        """
        IMPACT == 'HIGH'
        """,
        report: []  
    ],
    missense: [
        filter: [
            """
            str_detect(Consequence, "missense")
            """,
            """
            CADD  > 25                          |
            REVEL > 0.5                         |
            str_detect(CLIN_SIG, "pathogenic" ) |
            str_detect(SIFT    , "deleterious") |
            str_detect(PolyPhen, "damaging"   ) |
            str_detect(am_class, "pathogenic" )
            """
        ],
        report: [
            SIFT:'SIFT', 
            PolyPhen: 'PolyPhen', 
            AlphaMissense: 'AlphaMissense',
            REVEL: 'REVEL'
        ]
    ],
    splicing: [
        filter: [
            """
            ( 
                str_detect(Consequence, "splice") &
                IMPACT == 'MODERATE' 
            ) |
            SpliceAI_max > 0.33
            """
        ],
        report: [SpliceAI: 'SpliceAI']
    ],
    generic: [
        filter: [
            """
            (
                IMPACT == 'MODERATE' & 
                !str_detect(Consequence, "missense")
            )                                  |
            str_detect(CLIN_SIG, "pathogenic") |
            CADD      > 25                     |
            phyloP100 > 6
            """
        ],
        report: [SpliceAI: 'SpliceAI']
    ]
]

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

// additional options for cavalier R pacakge
params.cavalier_options = [
    database_mode: 'fallback',
    cache_dir: 'cavalier_cache'
]

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
