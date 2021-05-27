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
include { vcf_flatten_multi } from './proc/vcf_flatten_multi'
include { vcf_family_subset } from './proc/vcf_family_subset'
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
    // list of family members, starting with proband
    // for now only tested on singletons
    families = Channel.fromPath(params.ped).splitCsv(sep: '\t').map { [[it[1]]] }

    data |
        vcf_split |
        flatten |
        map { [it.name.replaceFirst(params.id + '-', '').replaceFirst('.vcf.gz', ''), it] } |
        vcf_flatten_multi |
        vep |
        vep_filter |
        toSortedList |
        transpose |
        toList |
        map { [params.id, it[1], it[2]] } |
        vcf_merge |
        map { it[1] } |
        combine(families) |
        vcf_family_subset

//    vcfanno |
//        annovar |
//        cavalier_prep |
//        toSortedList |
//        map { [params.id, it]} |
//        cavalier_merge

}