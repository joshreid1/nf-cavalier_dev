#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.vcf_input = ''
params.id = ''
params.n_split = 100

include { vcfsplit } from './proc/vcfsplit'
include { vcfanno } from './proc/vcfanno'
include { annovar } from './proc/annovar'
include { cavalier_prep } from './proc/cavalier_prep'
include { cavalier_merge } from './proc/cavalier_merge'


workflow {
    data = Channel.fromList([[params.id, file(params.vcf_input, checkIfExists:true)]])
        .flatMap { (1..params.n_split).collect { i -> [i] + it } }

    data | vcfsplit | vcfanno | annovar | cavalier_prep

    cavalier_prep.out.toSortedList().map { [params.id, it]} | cavalier_merge
}