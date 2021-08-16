

process cavalier {
    label 'C2M4T2'
    container 'jemunro/cavalier:dev'
    publishDir "output/cavalier_singleton", mode: 'copy'
    tag { fam }

    input:
        tuple val(fam), path(vcf), path(ped), val(sam), path(bam)

    output:
        tuple val(fam), file(out)

    script:
    out = "$fam"
    """
    cavalier_singleton_panelapp.R $vcf $out $sample=$bam \\
        --gene-lists $lists \\
        --maf-dom $params.maf_dom \\
        --maf-rec $params.maf_rec \\
        --maf-comp-het $params.maf_comp_het \\
        --gtex-rpkm $params.gtex_rpkm \\
        --omim-genemap2 $params.omim_genemap2 \\
        --max-cohort-af $params.max_cohort_af
    """
}

