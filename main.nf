#!/usr/bin/env nextflow
nextflow.enable.dsl=2
/*
TODO:
    - Add download URLs for reference files (similar to PLASTER) as default to simplify configuration
 */
//imput params
params.id = ''
params.outdir = 'output'
params.snp_vcf = ''
params.snp_caller = 'GATK'
params.fill_tags = false
params.remove_fields = 'INFO/CSQ'
params.sv_vcf = ''
params.ped = ''
params.bams = ''
params.lists = ''
params.chunk_size = 200000
params.sv_chunk_size = 10000

// filter params
params.max_af = 0.1
params.maf_dom = 0.0001
params.maf_de_novo = 0.0001
params.maf_rec = 0.01
params.maf_comp_het = 0.01
params.max_cohort_af = 1.0
params.max_cohort_ac = 'Inf'
params.min_impact = 'MODERATE'
params.exclude_benign_missense = true
params.include_sv_csv = true

// ref params
params.ref_fasta = ''
params.ref_hg38 = true
params.vep_cache = ''
params.vep_cache_ver = ''
params.vep_assembly = params.ref_hg38 ? 'GRCh38' : 'GRCh37'
params.pop_sv = ''
params.ref_gene = ''
params.cavalier_cache_dir = "$workflow.homeDir/.cavalier"

// exec params
params.sv_types = ['DEL', 'DUP', 'INS', 'INV', 'BND']
params.sv_type_match = [DEL: ['DEL'], DUP: ['CNV', 'DUP']]
params.no_slides = false // skip making slides in cavalier

include { vcf_channel; families_channel } from './nf/functions'
include { GetSamples } from './nf/GetSamples'
include { CheckInputs } from './nf/CheckInputs'
include { CleanAndChunk } from './nf/CleanAndChunk'
include { Annotate } from './nf/Annotate'
include { Report } from './nf/Report'

workflow {

    vcfs = vcf_channel()

    vcf_samples = GetSamples(vcfs)

    ann_vcf = vcf_samples |
        CheckInputs |
        combine(vcfs, by: 0) |
        CleanAndChunk |
        Annotate |
        combine(families_channel(vcf_samples), by: 0) |
        Report
}
