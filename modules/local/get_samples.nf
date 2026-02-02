
process GET_SAMPLES {
    label 'C2M2T2'
    label 'cavalier'

    input:
    path(bams)
    path(ped)
    path(short_vcf)
    path(struc_vcf)

    output:
    val  true                , emit: check
    path 'intersect_bams.tsv', emit: bams
    path 'intersect.ped'     , emit: ped
    path 'warnings.txt'      , emit: warnings
    
    script:
    """
    get_samples.R \\
        ${short_vcf.size() == 0 ? 'UNSET' : short_vcf} \\
        ${struc_vcf.size() == 0 ? 'UNSET' : struc_vcf} \\
        $bams \\
        intersect_bams.tsv \\
        ${ped.size() == 0 ? 'UNSET' : ped} \\
        intersect.ped 
    """
}

