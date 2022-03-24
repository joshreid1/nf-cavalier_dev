
include { pedigree_channel; bam_channel; pop_sv_channel; ref_gene_channel; make_path } from './functions'
include { Lists } from './Lists'

cavalier_cache_dir = make_path(params.cavalier_cache_dir)

workflow Report {
    take:
        input // set, vcf, tbi, fam, aff, unaff

    main:

   cavalier_input = input |
        map { it[[0,1,3,4,5]] } |
        family_subset |
        map { it[[1,0,2]] } | //fam, set, vcf
        combine(pedigree_channel(), by:0) |
        combine(bam_channel(), by:0) | //fam, set, vcf, ped, sam, bam, bai
        combine(Lists(), by:0) | //fam, set, vcf, ped, sam, bam, bai, lists
       map { it[[1,0,2,3,7,4,5,6]] + [cavalier_cache_dir] }   //set, fam, vcf, ped, lists, sam, bam, bai

    cavalier_input | cavalier

    // run svpv on candidate SVs
    cavalier.out |
        filter { it[0] == 'SV' } |
        filter { it[4].toFile().readLines().size() > 1 } |
        map { it[[1,3]] } | //fam, vcf
        combine( cavalier_input |
                     filter { it[0] == 'SV' } |
                     map { it[[1,5,6,7]] }, // fam, sam, bam, bai,
                 by:0 ) |
        combine(pop_sv_channel()) |
        combine(ref_gene_channel()) |
        svpv

    // report SNP and SV candidates
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
        collectFile(name: "${params.id}.SNP_candidates.csv", storeDir: 'output',
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
        collectFile(name: "${params.id}.SV_candidates.csv", storeDir: 'output',
                    newLine:true, sort: false, cache: false)

}

process family_subset {
    label 'C2M2T2'
    publishDir "output/family_subset", mode: 'copy'
    tag { "$fam:$set" }

    input:
    tuple val(set), path(vcf), val(fam), val(aff), val(unaff)

    output:
    tuple val(set), val(fam), path(out_vcf), path("${out_vcf}.tbi")

    script:
    out_vcf = "${set}.${fam}.subset.vcf.gz"
    """
    printf "${aff.join('\\n')}\\n" > aff
    bcftools view --no-update $vcf -Ou -s ${(aff + unaff).join(',')} |
        bcftools view  -i "GT[@aff]='alt'" -Oz -o $out_vcf
    bcftools index -t $out_vcf
    """
}

process cavalier {
    label 'C2M4T2'
//    container null
//    module 'R/4.1.2'
    publishDir "output/cavalier", mode: 'copy', pattern: "*.pptx"
    tag { "$fam:$set" }

    input:
    tuple val(set), val(fam), path(vcf), path(ped), path(lists), val(sam), path(bam), path(bai),
        path(cache_dir)

    output:
    tuple val(set), val(fam), path("${pref}.pptx"), path("${pref}.candidates.vcf.gz"), path("${pref}.candidates.csv")

    script:
    pref = "${params.id}.$fam.$set"
    sam_bam = [sam, bam instanceof List ? bam: [bam]]
        .transpose().collect {it.join('=') }.join(' ')
    flags =(
        (set == 'SV' ? ['--sv ']: []) +
            (params.exclude_benign_missense ? ['--exclude-benign-missense']: []) +
            (params.include_sv_csv ? ['--include-sv-csv']: [])
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
        --min-impact $params.min_impact \\
        --cache-dir $cache_dir
    """
}

process svpv {
    label 'C2M4T2'
    container = 'bahlolab/svpv:latest'
//    container null
//    conda '/stornext/Home/data/allstaff/m/munro.j/miniconda3/envs/numpy2'
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