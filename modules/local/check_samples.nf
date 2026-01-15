
process CHECK_SAMPLES {
    label 'C2M2T2'
    label 'cavalier'

    input:
    path(bams)
    path(ped)
    path(short_vcf)
    path(struc_vcf)

    output:
    val true
    
    script:
    """
    check_samples.R $bams \\
        ${ped.size() == 0 ? 'UNSET' : ped} \\
        ${short_vcf.size() == 0 ? 'UNSET' : short_vcf} \\
        ${struc_vcf.size() == 0 ? 'UNSET' : struc_vcf}
    """
}

