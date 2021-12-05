#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//imput params
params.mode = 'short'
params.id = ''
params.vcf = ''
params.ped = ''
params.bams = ''
params.lists = ''
params.chunk_size = params.mode == 'short' ? 200000 : 10000
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
// sv params
params.pop_sv = 'nstd166.GRCh38.variant_call.vcf.gz'
params.ref_gene = 'RefSeqGene.hg38.UCSC.txt'
params.sv_type_match = [DEL: ['DEL'], DUP: ['CNV', 'DUP']]
params.caller_id = 'SV'

include { path; read_tsv; get_families; date_ymd; checkMode } from './nf/functions'
include { check_inputs } from './nf/check_inputs'
include { prep_lists } from './nf/prep_lists'
include { vcf_sample_set } from './nf/vcf_sample_set'
include { vcf_count } from './nf/vcf_count'
include { split_intervals } from './nf/split_intervals'
include { vcf_split_norm } from './nf/vcf_split_norm'
include { vcf_split_sv_types } from './nf/vcf_split_sv_types'
include { vcf_split_sv_types as vcf_split_sv_types_pop } from './nf/vcf_split_sv_types'
include { vcf_stub } from './nf/vcf_stub'

include { vep; vep_svo } from './nf/vep'
include { vcf_merge } from './nf/vcf_merge'
include { vcf_merge as vcf_merge_ao } from './nf/vcf_merge' addParams(allow_overlap:true)
include { vcf_family_subset } from './nf/vcf_family_subset'
include { cavalier; cavalier_sv } from './nf/cavalier'
include { svpv } from './nf/svpv'

checkMode(params.mode)
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

if (params.mode == 'sv') {
    pop_sv = path(params.pop_sv)
    pop_sv_tbi = path(params.pop_sv + '.tbi')
    sv_type_match_rev = params.sv_type_match
        .collectMany { k, v -> v.collect { [it, k] }}
        .groupBy { it[0] }
        .collectEntries { k, v -> [(k) : v.collect{ it[1] }.unique()] }
    ref_gene = path(params.ref_gene)
}

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

    list_channel = prep_lists(lists)

    split_vcf =
        Channel.value([vcf, tbi]) |
        vcf_count |
        map { Math.ceil((it.toFile().text as int) / params.chunk_size) as int } |
        map { [it, ref_fai, gaps] } |
        split_intervals |
        flatten() |
        map { [it, vcf, tbi, ref_fa, ref_fai] } |
        vcf_split_norm |
        map { it.collect { it instanceof List ? it : [it] }} |
        flatMap { it.transpose() }

    if (params.mode == 'short') {

//        annot_vcf =
//            split_vcf |
//                map { [it, ref_fa, ref_fai, vep_cache] } |
//                vep |
//                flatMap { [['vep', it[0]],  ['vep-modifier', it[1]], ['unannotated', it[2]]] } |
//                collectFile(newLine: true, sort: { new File(it).toPath().fileName.toString() } ) {
//                    ["${it[0]}.files.txt", it[1].toString()] } |
//                map { [it.name.replaceAll('.files.txt', ''), it] } |
//                vcf_merge |
//                filter { it[0] == 'vep' } |
//                map { it[1..2] } |
//                first()

    } else if (params.mode == 'sv') {

        pop_sv_split = Channel.value([pop_sv, pop_sv_tbi]) |
            vcf_split_sv_types_pop |
            flatMap { it.transpose() } |
            map { [(it[0].name =~ /([A-Z]+)\.vcf\.gz$/)[0][1]] + it } |
            filter { sv_type_match_rev.keySet().contains(it[0]) } |
            flatMap { sv_type_match_rev[it[0]].collect {type -> [type] + it[1..2] } } |
            mix(vcf_stub(pop_sv) | map { ['STUB'] + it })
//
        annot_vcf =
            split_vcf |
                vcf_split_sv_types |
                flatMap { it.transpose() } |
                map { [(it[0].name =~ /([A-Z]+)\.vcf\.gz$/)[0][1]] + it } |
                map { params.sv_type_match.keySet().contains(it[0]) ?
                    it : ['STUB'] + it[1..2] } |
                combine(pop_sv_split, by:0) |
                groupTuple(by: 0..2) |
                map { it + [ref_fa, ref_fai, vep_cache] } |
                vep_svo |
                flatMap { [['vep', it[1]], ['unannotated', it[2]]] } |
                collectFile(newLine: true, sort: true ) { ["${it[0]}.files.txt", it[1].toString()] } |
                map { [it.name.replaceAll('.files.txt', ''), it] } |
                vcf_merge_ao |
                filter { it[0] == 'vep' } |
                map { it[1] } |
                first()
    }

    calalier_input = annot_vcf |
        combine(families) |
        vcf_family_subset |
        map { it[0..1] } |
        combine(ped_channel, by:0) |
        combine(list_channel, by:0) |
        combine(bam_channel, by:0)

    if (params.mode == 'sv') {
        calalier_input |
            cavalier_sv

        candidates = cavalier_sv.out |
            map { it[3] } |
            splitCsv(header: true)

        cavalier_sv.out |
            map { it[0, 2] } |
            combine(candidates.map{it.family}.unique(), by:0) |
            combine(bam_channel, by:0) |
            map { it +  [pop_sv, pop_sv_tbi, ref_gene] } |
            svpv

        candidates |
            first |
            map { (it.keySet() as List).join(',') } |
            concat(candidates.map { (it.values() as List).join(',') }) |
            collectFile(name: 'candidates.csv', storeDir: 'output',
                        newLine:true, sort: false, cache: false)
    }


}