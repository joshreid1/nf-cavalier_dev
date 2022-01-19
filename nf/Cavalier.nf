
include { pop_sv_channel; ref_gene_channel } from './functions'

workflow Cavalier {
    take:
        input  //set, fam, vcf, ped, lists, sam, bam, bai

    main:
    // run cavalier
    cavalier(input)
    // run svpv on candidate SVs
    cavalier.out |
        filter { it[0] == 'SV' } |
        filter { it[4].toFile().readLines().size() > 1 } |
        map { it[[1,3]] } | //fam, vcf
        combine( input |
                     filter { it[0] == 'SV' } |
                     map { it[[1,5,6,7]] }, // fam, sam, bam, bai,
                 by:0 ) |
        combine(pop_sv_channel()) |
        combine(ref_gene_channel()) |
        svpv

    candidates = cavalier.out |
        map { it[4] } |
        splitCsv(header: true) |
        branch { snp: it.set == 'SNP'; sv: true }


    candidates.snp |
        first |
        map { (it.keySet() as List).join(',') } |
        concat(
            candidates.snp |
                map { (it.values() as List).join(',') } |
                toSortedList |
                flatten
        ) |
        collectFile(name: 'SNP_candidates.csv', storeDir: 'output',
                    newLine:true, sort: false, cache: false)

    candidates.sv |
        first |
        map { (it.keySet() as List).join(',') } |
        concat(
            candidates.sv |
                map { (it.values() as List).join(',') } |
                toSortedList |
                flatten
        ) |
        collectFile(name: 'SV_candidates.csv', storeDir: 'output',
                    newLine:true, sort: false, cache: false)

//    emit: output // set, fam, pptx, vcf, csv
}

process cavalier {
    label 'C2M4T2'
    label 'Cavalier'
    container null
    module 'R/3.6.1'
    publishDir "output/cavalier", mode: 'copy', pattern: "*.pptx"
    tag { "$fam:$set" }

    input:
    tuple val(set), val(fam), path(vcf), path(ped), path(lists), val(sam), path(bam), path(bai)

    output:
    tuple val(set), val(fam), path("${pref}.pptx"), path("${pref}.candidates.vcf.gz"), path("${pref}.candidates.csv")

    script:
    pref = "$fam.$set"
    sam_bam = [sam, bam instanceof List ? bam: [bam]]
        .transpose().collect {it.join('=') }.join(' ')
    flags =(
        (set == 'SV' ? ['--sv ']: []) +
        (params.exclude_benign_missense ? ['--exclude-benign-missense']: [])
    ).join(' ')
    """
    cavalier_wrapper.R $vcf $ped $sam_bam $flags \\
        --out $pref \\
        --family $fam \\
        --genome ${params.ref_hg38 ? 'hg38' : 'hg19'} \\
        --gene-lists ${lists.join(',')} \\
        --maf-dom $params.maf_dom \\
        --maf-de-novo $params.maf_de_novo \\
        --maf-rec $params.maf_rec \\
        --maf-comp-het $params.maf_comp_het \\
        --max-cohort-af $params.max_cohort_af \\
        --max-cohort-ac $params.max_cohort_ac \\
        --min-impact $params.min_impact
    """
}

process svpv {
    label 'C2M4T2'
    container null
    conda '/stornext/Home/data/allstaff/m/munro.j/miniconda3/envs/numpy2'
    publishDir "output/svpv", mode: 'copy'
    tag { fam }

    input:
    tuple val(fam), path(vcf), val(sam), path(bam), path(bai), path(pop_sv), path(pop_sv_indx), path(ref_gene)

    output:
    tuple val(fam), path(fam)

    script:
    sam_bam = [sam, bam instanceof List ? bam : [bam]]
        .transpose().collect { it.join('=') }.join(' ')
    """
    SVPV \\
        -o $fam \\
        -samples ${sam.join(',')} \\
        -aln ${bam.join(',')} \\
        -vcf $vcf \\
        -ref_vcf gnomAD:$pop_sv \\
        -ref_gene $ref_gene
    """
}