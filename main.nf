#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.vcf_input = ''
params.id = ''
params.n_split = 100
params.ped = 'samples.ped'

include { vcf_split } from './proc/vcf_split'
include { vep } from './proc/vep'
include { vep_filter } from './proc/vep_filter'
include { vcf_merge } from './proc/vcf_merge'
include { vcf_sample_subset } from './proc/vcf_sample_subset'
include { vcfanno } from './proc/vcfanno'
include { annovar } from './proc/annovar'
include { cavalier_prep } from './proc/cavalier_prep'
include { cavalier_merge } from './proc/cavalier_merge'


workflow {
    data = Channel.fromList([
        [params.id,
         file(params.vcf_input, checkIfExists:true),
         file(params.vcf_input + '.tbi', checkIfExists:true)]
    ])

    samples = Channel.fromPath(params.ped).splitCsv(sep: '\t').map { it[1] }

    data |
        vcf_split |
        flatten |
        map { [it.name.replaceFirst(params.id + '-', '').replaceFirst('.vcf.gz', ''), it] } |
        ( vcfanno & vep )

    vcfanno.out |
        annovar |
        cavalier_prep |
        toSortedList |
        map { [params.id, it]} |
        cavalier_merge

    vep.out |
        vep_filter |
        toSortedList |
        transpose |
        toList |
        map { [params.id, it[1], it[2]] } |
        vcf_merge |
        map { it[1] } |
        combine(samples) |
        vcf_sample_subset
}