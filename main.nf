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
params.maf_dom = 0.0001
params.maf_de_novo = 0.0001
params.maf_rec = 0.01
params.maf_comp_het = 0.01
params.max_cohort_af = 1.0
params.min_impact = 'MODERATE'
params.exclude_benign_missense = true
// ref params
params.ref_fasta = ''
params.ref_hg38 = true
params.vep_cache = ''
params.vep_cache_ver = ''
params.vep_assembly = params.ref_hg38 ? 'GRCh38' : 'GRCh37'
// exec params
params.chunk_size = 200000

include { path; read_tsv; get_families; date_ymd } from './nf/functions'
include { check_inputs } from './nf/check_inputs'
include { update_list_versions } from './nf/update_list_versions'
include { pull_latest_list } from './nf/pull_latest_list'
include { vcf_sample_set } from './nf/vcf_sample_set'
include { vcf_count } from './nf/vcf_count'
include { split_intervals } from './nf/split_intervals'
include { vcf_split_norm } from './nf/vcf_split_norm'

include { vep } from './nf/vep'
include { vcf_merge } from './nf/vcf_merge'
include { vcf_family_subset } from './nf/vcf_family_subset'
include { cavalier } from './nf/cavalier'

vcf = path(params.vcf)
tbi = path(params.vcf + '.tbi')
ped = read_tsv(path(params.ped), ['fid', 'iid', 'pid', 'mid', 'sex', 'phe'])
bams = read_tsv(path(params.bams), ['iid', 'bam'])
lists = read_tsv(path(params.lists), ['fid', 'list'])
ref_fa = path(params.ref_fasta)
ref_fai = path(params.ref_fasta + '.fai')
gaps = params.ref_hg38 ?
    path("${workflow.projectDir}/data/hg38.gaps.bed.gz") :
    path("${workflow.projectDir}/data/hg19.gaps.bed.gz")
vep_cache = path(params.vep_cache)

workflow {

    vcf_samples = vcf_sample_set(vcf) | map { it.toFile().readLines() as ArrayList }

    check_inputs(ped, bams, lists, vcf_samples)

    families = vcf_samples |
        flatMap { sam ->
            get_families(ped).collect { fam, af, un ->
                [fam, af.intersect(sam), un.intersect(sam)] }
            .findAll { it[1].size() > 0 }
        }

    ped_channel = Channel.from(ped) |
        map { it.values() as ArrayList } |
        collectFile(newLine:true) {
            [ "${it[0]}.ped", it.join('\t')]
        } |
        map { [it.name.replaceAll('.ped', ''), it] }

    bam_channel = Channel.from(bams) |
        map { [it.iid, path(it.bam), path(it.bam + '.bai')] } |
        combine(ped.collect { [it.iid, it.fid] }, by: 0) |
        map { it[[3,0,1,2]] } |
        groupTuple(by: 0)

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
                update_list_versions |
                splitCsv(sep: '\t', skip: 1, strip: true) |
                pull_latest_list |
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

    Channel.value([vcf, tbi]) |
        vcf_count |
        map { Math.ceil((it.toFile().text as int) / params.chunk_size) as int } |
        map { [it, ref_fai, gaps] } |
        split_intervals |
        flatten() |
        map { [it, vcf, tbi, ref_fa, ref_fai] } |
        vcf_split_norm |
        flatten() |
        map { [it, ref_fa, ref_fai, vep_cache] } |
        vep |
        flatMap { [['vep', it[0]],  ['vep-modifier', it[1]], ['unannotated', it[2]]] } |
        collectFile(newLine: true, sort: { new File(it).toPath().fileName.toString() } ) {
            ["${it[0]}.files.txt", it[1].toString()] } |
        vcf_merge |
        filter { it[0] == 'vep' } |
        map { it[1] } |
        first() |
        combine(families) |
        vcf_family_subset |
        map { it[0..1] } |
        combine(ped_channel, by:0) |
        combine(list_channel, by:0) |
        combine(bam_channel, by:0) |
        cavalier
}