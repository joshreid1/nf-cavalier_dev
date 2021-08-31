
process cavalier {
    label 'C2M4T2'
    container 'bahlolab/cavalier:latest'
//    container null
//    module 'R/3.6.1'
    publishDir "output/cavalier", mode: 'copy'
    tag { fam }

    input:
        tuple val(fam), path(vcf), path(ped), path(lists), val(sam), path(bam), path(bai), path(genemap2)

    output:
        tuple val(fam), path("${fam}.pptx")

    script:
    sam_bam = [sam, bam instanceof List ? bam: [bam]]
        .transpose().collect {it.join('=') }.join(' ')
    """
    cavalier_wrapper.R $vcf $ped $sam_bam \\
        --out $fam \\
        --gene-lists ${lists.join(',')} \\
        --omim-genemap2 $genemap2 \\
        --maf-dom $params.maf_dom \\
        --maf-de-novo $params.maf_de_novo \\
        --maf-rec $params.maf_rec \\
        --maf-comp-het $params.maf_comp_het \\
        --max-cohort-af $params.max_cohort_af \\
        --min-impact $params.min_impact
    """
}

