#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/* TODO:
    inputs:
        - vcf
        - bam_manifest [sample, bam]
            - csv file
        - pedigree [family, sample, pid, mid, sex, phenotype]
            - tsv file, no column names
        - gene_lists [family, list_name, list_path]
            - csv file, one list per line
    gene_lists:
        - csv file, single required column 'gene'
        - additional metadata columns will be reported by cavalier
    checks:
        - check pedigree (don't necessarily expected all family members to have sample)
            - check at least one affected sample per family in VCF
        - check samples in vcf, pedigree and bam_manifest
        - check families in gene lists
        - report similar to peddy
    pedigree:
        - function read_ped(), return list of map
        - function group_ped(), return list of [family, [affected...], [unaffected...]]
        - collectFile() to split into sub pedigrees for input to cavalier
        - multipe models can be used for a given family by using different family_ids for each
    bam_manifest:
        - combine with pedigree to add in family
 */

params.vcf_input = ''
params.id = ''
params.n_split = 100
params.sample_manifest = 'sample_manifest2.tsv'

include { read_tsv } from './functions.nf'

include { vcf_split } from './tasks/vcf_split'
include { vcf_flatten_multi } from './tasks/vcf_flatten_multi'
include { vep } from './tasks/vep'
include { vep_filter } from './tasks/vep_filter'
include { vcf_merge } from './tasks/vcf_merge'
include { vcf_family_subset } from './tasks/vcf_family_subset'
include { cavalier_singleton } from './tasks/cavalier_singleton'

workflow {
    data = Channel.fromList([
        [params.id,
         file(params.vcf_input, checkIfExists:true),
         file(params.vcf_input + '.tbi', checkIfExists:true)]
    ])
    sample_manifest = read_tsv(
        file(params.sample_manifest, checkIfExists: true, type: 'file'),
        ['sample', 'bam', 'lists']
    )
    // list of family members, starting with proband
    // ultimately will extract from pedigree, for now only working on singletons
    families = Channel.from(sample_manifest).map { [[it.sample]] }
    samples = Channel.from(sample_manifest)
        .map {
            [ it.sample,
              file(it.bam, checkIfExists:true, type: 'file'),
              file(it.bam + '.bai', checkIfExists:true, type: 'file'),
              it.lists
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
}