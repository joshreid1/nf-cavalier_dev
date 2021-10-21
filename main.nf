#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//imput params
params.id = ''
params.vcf = ''
params.ped = ''
params.bams = ''
params.lists = ''
// filter params
params.max_af = 0.10
params.vep_impact = ['LOW', 'MODERATE', 'HIGH']
params.maf_dom = 0.0001
params.maf_de_novo = 0.0001
params.maf_rec = 0.01
params.maf_comp_het = 0.01
params.max_cohort_af = 1.0
params.min_impact = 'MODERATE'
// ref params
params.ref_fasta = ''
params.ref_hg38 = true
params.vep_cache = ''
params.vep_cache_ver = ''
params.vep_assembly = params.ref_hg38 ? 'GRCh38' : 'GRCh37'
// exec params
params.n_split = 100

include { path; read_tsv; get_families; date_ymd } from './nf/functions'

include { vcf_sample_list } from './nf/vcf_sample_list'
include { check_latest_version } from './nf/check_latest_version'
include { get_latest_version } from './nf/get_latest_version'
include { vcf_split } from './nf/vcf_split'
include { vcf_flatten_multi } from './nf/vcf_flatten_multi'
include { vep } from './nf/vep'
include { vep_filter } from './nf/vep_filter'
include { vcf_merge } from './nf/vcf_merge'
include { vcf_family_subset } from './nf/vcf_family_subset'
include { cavalier } from './nf/cavalier'

vcf = path(params.vcf)
tbi = path(params.vcf + '.tbi')
ped = read_tsv(path(params.ped), ['fid', 'iid', 'pid', 'mid', 'sex', 'phe'])
bams = read_tsv(path(params.bams), ['iid', 'bam'])
lists = read_tsv(path(params.lists), ['fid', 'list'])
ref_fasta = path(params.ref_fasta)
ref_fai = path(params.ref_fasta + '.fai')
vep_cache = path(params.vep_cache)

workflow {

    families = vcf_sample_list(vcf) |
        map { [it.toFile().readLines() as ArrayList] } |
        combine(get_families(ped)) |
        map { sm, fam, af, un -> [fam, af.intersect(sm), un.intersect(sm)] } |
        // Note: families silently dropped here if no affected members in VCF
        filter { it[1].size() > 0 }

    ped_channel = Channel.from(ped) |
        map { it.values() as ArrayList } |
        collectFile(newLine:true) {
            [ "${it[0]}.ped", it.join('\t')]
        } |
        map { [it.name.replaceAll('.ped', ''), it] }

    if (lists.any { it.list  ==~ '^(HP|PA[AE]):.+'}) {
        lists = Channel.from(lists) |
            map { it.values() as ArrayList } |
            branch {
                web: it[1]  ==~ '^(HP|PA[AE]):.+'
                file: true
            }
        list_channel =
            lists.web |
                map { it[1] } |
                unique |
                collectFile(name: 'list_ids.txt', newLine: true) |
                combine([date_ymd()]) |
                check_latest_version |
                splitCsv(sep: '\t', skip: 1, strip: true) |
                get_latest_version |
                combine(lists.web.map {it[[1,0]]}, by:0) |
                map { [it[2], it[1]] } |
                mix(lists.file.map { [it[0], path(it[1])] }) |
                groupTuple(by: 0)
    } else {
        list_channel = Channel.from(lists) |
            map { it.values() as ArrayList } |
            map { [it[0], path(it[1])] } |
            groupTuple(by: 0)
    }

    bam_channel = Channel.from(bams) |
        map { [it.iid, path(it.bam), path(it.bam + '.bai')] } |
        combine(ped.collect { [it.iid, it.fid] }, by: 0) |
        map { it[[3 ,0,1, 2]] } |
        groupTuple(by: 0)

    Channel.value([vcf, tbi]) |
        vcf_split |
        flatten |
        map { [it.name.replaceFirst(params.id + '-', '').replaceFirst('.vcf.gz', ''), it] } |
        vcf_flatten_multi |
        combine([[ref_fasta, ref_fai, vep_cache]]) |
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
        map { it[0..1] } |
        combine(ped_channel, by:0) |
        combine(list_channel, by:0) |
        combine(bam_channel, by:0) |
        cavalier
}