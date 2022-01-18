
workflow Cavalier {
    take:
        input  //set, fam, vcf, ped, lists, sam, bam, bai

    main:
        output = cavalier(input)
//        candidates = cavalier_sv.out |
//            map { it[3] } |
//            splitCsv(header: true)
//
//        cavalier_sv.out |
//            map { it[0, 2] } |
//            combine(candidates.map{it.family}.unique(), by:0) |
//            combine(bam_channel, by:0) |
//            map { it +  [pop_sv, pop_sv_tbi, ref_gene] } |
//            svpv
//
//        candidates |
//            first |
//            map { (it.keySet() as List).join(',') } |
//            concat(candidates.map { (it.values() as List).join(',') }) |
//            collectFile(name: 'candidates.csv', storeDir: 'output',
//                        newLine:true, sort: false, cache: false)

    emit: output // set, fam, pptx, vcf, csv
}

process cavalier {
    label 'C2M4T2'
    label 'Cavalier'
    container null
    module 'R/3.6.1'
    publishDir "output/cavalier", mode: 'copy', pattern: "*.pptx"
    tag { "$set:$fam" }

    input:
    tuple val(set), val(fam), path(vcf), path(ped), path(lists), val(sam), path(bam), path(bai)

    output:
    tuple val(set), val(fam), path("${pref}.pptx"), path("${pref}.candidates.vcf.gz"), path("${pref}.candidates.csv")

    script:
    pref = "$set.$fam"
    sam_bam = [sam, bam instanceof List ? bam: [bam]]
        .transpose().collect {it.join('=') }.join(' ')
    flags =(
        (set == 'SV' ? ['--sv ']: []) +
        (params.exclude_benign_missense ? ['--exclude-benign-missense']: [])
    ).join(' ')
    """
    cavalier_wrapper.R $vcf $ped $sam_bam $flags \\
        --out $pref \\
        --genome ${params.ref_hg38 ? 'hg38' : 'hg19'} \\
        --gene-lists ${lists.join(',')} \\
        --maf-dom $params.maf_dom \\
        --maf-de-novo $params.maf_de_novo \\
        --maf-rec $params.maf_rec \\
        --maf-comp-het $params.maf_comp_het \\
        --max-cohort-af $params.max_cohort_af \\
        --min-impact $params.min_impact
    """
}
