#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.vcf_input = ''
params.id = ''
params.n_split = 100
params.sample_manifest = 'sample_manifest.tsv'

include { vcf_split } from './tasks/vcf_split'
include { vcf_flatten_multi } from './tasks/vcf_flatten_multi'
include { vep } from './tasks/vep'
include { vep_filter } from './tasks/vep_filter'
include { vcf_merge } from './tasks/vcf_merge'
include { vcf_family_subset } from './tasks/vcf_family_subset'
include { cavalier_singleton } from './tasks/cavalier_singleton'

//include { vcfanno } from './tasks/vcfanno'
//include { annovar } from './tasks/annovar'
//include { cavalier_prep } from './tasks/cavalier_prep'

header = ['sample', 'bam', 'lists']
sample_manifest = file(params.sample_manifest, checkIfExists: true, type: 'file')
    .toFile().readLines().with { lines ->
        assert header == lines[0].split('\t')
        lines.drop(1).collect { [header, it.split('\t')].transpose().collectEntries { k, v -> [(k): v] } }
    }

workflow {
    data = Channel.fromList([
        [params.id,
         file(params.vcf_input, checkIfExists:true),
         file(params.vcf_input + '.tbi', checkIfExists:true)]
    ])
    // list of family members, starting with proband
    // ultimately will extract from pedigree, for now only working on singletons
//    families = Channel.fromPath(params.ped).splitCsv(sep: '\t').map { [[it[1]]] }
    families = Channel.from(sample_manifest).map { [[it.sample]] }
    samples = Channel.from(sample_manifest)
        .map {
            [ it.sample,
              file(it.bam, checkIfExists:true, type: 'file'),
              it.lists.split(',').collect {file(it, checkIfExists:true, type: 'file') }
            ]
        }

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
        vcf_family_subset |
        join(samples) |
        cavalier_singleton

//    vcfanno |
//        annovar |
//        cavalier_prep |
//        toSortedList |
//        map { [params.id, it]} |
//        cavalier_merge

}