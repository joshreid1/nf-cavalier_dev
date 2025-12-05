#!/usr/bin/env nextflow

// main input params
params.outdir = 'output'
params.ped = ''
params.bams = ''
params.lists = ''
params.short_vcf = ''
params.short_vcf_annotated = null // skip annotation by providing pre-annotated VCF (output from another run)
params.sv_vcf = '' // not implemented
params.sv_vcf_annotated = null // not implemented
params.annotate_only = false // will not filter and report variants

// ref fasta
params.ref_fasta = '/vast/projects/bahlo_cache/ref_genome/hg38_GATK/Homo_sapiens_assembly38.fasta'

// used by cavalier R package - can point to a unified cache else will create in working directory
params.cavalier_cache_dir = 'cavalier_cache'

// additional options for cavalier R pacakge
params.cavalier_options = [ // will try to get latest db/vers, but fallback to cache if unavailable (set to 'offline' no check)
    read_only_cache_dir: '/vast/projects/bahlo_cache/cavalier_cache' // will use cache files from but won't write to
]

/* =================== SHORT/Indel ARGS =================== */
// split input VCF into shards for parallel processing
params.short_n_shards = 200
// apply filter to input variants
params.short_vcf_filter = "PASS,."
// INFO fields to keep from VCF
params.short_info = ['AC', 'AF', 'AN']
// FORMAT fields to keep from VCF
params.short_format = ['GT', 'GQ', 'DP'] // TODO: ensure 'GT' is always present
// fill AC, AF, and AN from VCF - recommended to set true if VCF does not have these set
params.short_fill_tags = false
// vcfanno settings
params.vcfanno_binary = 'https://github.com/brentp/vcfanno/releases/download/v0.3.7/vcfanno_linux64'
params.short_vcfanno = [
    [   vcf: '/vast/projects/bahlo_cache/annotation/gnomAD/joint_sites_4.1.vcf.gz', 
        csi: true, 
        fields: [gnomad_AF: 'AF', gnomad_AC: 'AC', gnomad_FAF95: 'fafmax_faf95_max', gnomad_nhomalt: 'nhomalt']
    ],
    [   tsv: '/vast/projects/bahlo_cache/annotation/CADD/1.7/whole_genome_SNVs.tsv.gz',
        fields: [CADD: 6]
    ],
    [   tsv: '/vast/projects/bahlo_cache/annotation/CADD/1.7/gnomad.genomes.r4.0.indel.tsv.gz',
        csi: true,
        fields: [CADD: 6]
    ],
    [   tsv: '/vast/projects/bahlo_cache/annotation/CADD/1.7/gnomad.joint.r4.1.indel.tsv.bgz',
        csi: true,
        fields: [CADD: 6]
    ],
    [   tsv: '/vast/projects/bahlo_cache/annotation/phyloP/hg38.phyloP100way.bed.gz',
        fields: [phyloP100: 4],
        op: 'max'
    ],
    [   vcf: '/vast/projects/bahlo_cache/annotation/ClinVar/clinvar_20251123.vcf.gz',
        csi: true,
        fields: [CLNSIG: 'CLNSIG', CLNGENE: 'GENEINFO',  CLNVID: 'ID']
    ]
]

// filter to apply after VCF anno, drastically reduces VEP runtime if we drop common variants
params.short_vcfanno_filter = 'gnomad_AF<0.01 || gnomad_AF="."' // (keep if AF < 0.01 or AF is missing)
// VEP settings
params.vep_cache     = '/vast/projects/bahlo_cache/vep_cache'
params.vep_cache_ver = '115'
// check number of variants output by VEP equal to number input
params.vep_check     = true
// VEP plugins (setting to null will disable plugin)
params.vep_spliceai_snv   = '/vast/projects/bahlo_cache/annotation/spliceAI/1.3/spliceai_scores.raw.snv.hg38.vcf.gz'
params.vep_spliceai_indel = '/vast/projects/bahlo_cache/annotation/spliceAI/1.3/spliceai_scores.raw.indel.hg38.vcf.gz'
params.vep_alphamissense  = '/vast/projects/bahlo_cache/annotation/alphamissense/AlphaMissense_hg38.tsv.gz'
params.vep_revel          = '/vast/projects/bahlo_cache/annotation/REVEL/revel_1.3.hg38.vep.tsv.gz'
params.vep_utr_annotator  = '/vast/projects/bahlo_cache/annotation/UTRannotator/uORF_5UTR_GRCh38_PUBLIC.txt'

// filters, can be disabled by setting = null
params.FILTER_SHORT_MIN_DP = 5
params.FILTER_SHORT_MIN_GQ = 10
// short variant population (gnomAD) filters
params.FILTER_SHORT_POP_DOM_MAX_AF = 0.0001
params.FILTER_SHORT_POP_REC_MAX_AF = 0.01
params.FILTER_SHORT_POP_DOM_MAX_AC = 20
params.FILTER_SHORT_POP_REC_MAX_AC = 100
params.FILTER_SHORT_POP_DOM_MAX_HOM = 5
params.FILTER_SHORT_POP_REC_MAX_HOM = 10
// short variant cohort filters
params.FILTER_SHORT_COH_DOM_MAX_AF = null
params.FILTER_SHORT_COH_REC_MAX_AF = null
params.FILTER_SHORT_COH_DOM_MAX_AC = null
params.FILTER_SHORT_COH_REC_MAX_AC = null
// ClinVar significance filters (regex)
params.FILTER_SHORT_CLINVAR_KEEP_PAT  = '(p|P)athogenic(?!ity)'
params.FILTER_SHORT_CLINVAR_DISC_PAT  = '(b|B)enign'
// misc filters
params.FILTER_SHORT_CLINVAR_ANYWHERE = true
params.FILTER_SHORT_LOF              = true
params.FILTER_SHORT_MISSENCE         = true
params.FILTER_SHORT_SPLICING         = true
params.FILTER_SHORT_MIN_CADD_PP      = 25.3 // Clingen PP3 supporting  doi: 10.1016/j.ajhg.2022.10.013
params.FILTER_SHORT_MIN_SPLICEAI_PP  = 0.20 //Clingen PP3 supporting https://doi.org/10.1016/j.ajhg.2023.06.002
params.FILTER_SHORT_VEP_MIN_IMPACT   = 'MODERATE'
params.FILTER_SHORT_VEP_CONSEQUENCES = null

// structural variant population (gnomAD) filters
params.FILTER_STRUC_POP_DOM_MAX_AF = 0.0001
params.FILTER_STRUC_POP_REC_MAX_AF = 0.01
params.FILTER_STRUC_POP_DOM_MAX_AC = 20
params.FILTER_STRUC_POP_REC_MAX_AC = 100
params.FILTER_STRUC_POP_DOM_MAX_HOM = null
params.FILTER_STRUC_POP_REC_MAX_HOM = null
// strucutural  variant cohort filters
params.FILTER_STRUC_COH_DOM_MAX_AF = 0.01
params.FILTER_STRUC_COH_REC_MAX_AF = 0.01
params.FILTER_STRUC_COH_DOM_MAX_AC = null
params.FILTER_STRUC_COH_REC_MAX_AC = null

// fields to report for all variants
params.short_report_fields =  [
    Gene         : 'Gene',
    Inheritance  : 'Inheritance', 
    Consequence  : 'Consequence',
    HGVS         : 'HGVS',
    ClinVar      : 'CLNSIG',
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

include { SETUP    } from './workflows/setup'
include { ANNOTATE } from './workflows/annotate'
include { CAVALIER } from './workflows/cavalier'


workflow {

    validate_params()

    SETUP()

    ANNOTATE()

    if (!params.annotate_only) {
        CAVALIER(
            SETUP.out.lists,
            SETUP.out.cavalier_opts,
            ANNOTATE.out.short_vcf,
            ANNOTATE.out.struc_vcf
        )
    }

}
