
process cavalier {
    label 'C2M4T2'
    label 'cavalier'
//    container = null
//    module = 'R/3.6.1'
    publishDir "output/cavalier", mode: 'copy'
    tag { fam }

    input:
        tuple val(fam), path(vcf), path(ped), path(lists), val(sam), path(bam), path(bai)

    output:
        tuple val(fam), path("${fam}.pptx")

    script:
    sam_bam = [sam, bam instanceof List ? bam: [bam]]
        .transpose().collect {it.join('=') }.join(' ')
    """
    cavalier_wrapper.R $vcf $ped $sam_bam \\
        --out $fam \\
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

