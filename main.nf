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
   tweaks:
        - pre-filter VCF by select only samples in manifest/pedigree and only alternate alleles
 */

params.id = ''
params.vcf = ''
params.ped = ''
params.bams = ''
params.lists = ''
params.n_split = 100
params.vep_cache = ''
params.vep_cache_ver = ''
params.vep_assembly = ''
params.max_af = 0.10
params.vep_impact = ['MODERATE', 'HIGH']

include { path; read_tsv; get_families } from './nf/functions'

include { vcf_sample_list } from './nf/vcf_sample_list'
include { vcf_split } from './nf/vcf_split'
include { vcf_flatten_multi } from './nf/vcf_flatten_multi'
include { vep } from './nf/vep'
include { vep_filter } from './nf/vep_filter'
include { vcf_merge } from './nf/vcf_merge'
include { vcf_family_subset } from './nf/vcf_family_subset'
include { cavalier_singleton } from './nf/cavalier_singleton'

vcf = path(params.vcf)
tbi = path(params.vcf + '.tbi')
ped = read_tsv(path(params.ped), ['fid', 'iid', 'pid', 'mid', 'sex', 'phe'])
bams = read_tsv(path(params.bams), ['iid', 'bam'])
lists = read_tsv(path(params.lists), ['fid', 'list'])

workflow {

//    data = Channel.fromList([
//        [params.id,
//         file(params.vcf_input, checkIfExists:true),
//         file(params.vcf_input + '.tbi', checkIfExists:true)]
//    ])
//    sample_manifest = read_tsv(
//        file(params.sample_manifest, checkIfExists: true, type: 'file'),
//        ['sample', 'bam', 'lists']
//    )
//    // list of family members, starting with proband
//    // ultimately will extract from pedigree, for now only working on singletons
//    families = Channel.from(sample_manifest).map { [[it.sample]] }
//    samples = Channel.from(sample_manifest)
//        .map {
//            [ it.sample,
//              file(it.bam, checkIfExists:true, type: 'file'),
//              file(it.bam + '.bai', checkIfExists:true, type: 'file'),
//              it.lists
//            ]
//        }
    families = vcf_sample_list(vcf) |
        map { [it.toFile().readLines() as ArrayList] } |
        combine(get_families(ped)) |
        map { sm, fam, af, un -> [fam, af.intersect(sm), un.intersect(sm)] } |
        // Note: families silently dropped here if no affected members in VCF
        filter { it[1].size() > 0 }

    Channel.value([vcf, tbi]) |
        vcf_split |
        flatten |
        map { [it.name.replaceFirst(params.id + '-', '').replaceFirst('.vcf.gz', ''), it] } |
        vcf_flatten_multi |
        vep |
        vep_filter |
        toSortedList |
        transpose |
        toList |
        map { it[1..2] } |
        vcf_merge |
        map { it[0] } |
        combine(families) |
        vcf_family_subset |
        view
//        join(samples) |
//        cavalier_singleton
}