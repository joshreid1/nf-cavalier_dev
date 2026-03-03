
process GET_SAMPLES {
    label 'C2M2T2'
    label 'cavalier'

    input:
    path(alignments)
    path(ped)
    path(short_vcf)
    path(struc_vcf)

    output:
    val  true                , emit: check
    path 'intersect_alignments.tsv', emit: alignments
    path 'intersect.ped'     , emit: ped
    path 'warnings.txt'      , emit: warnings
    
    script:
    """
    get_samples.R \\
        ${short_vcf.size() == 0 ? 'UNSET' : short_vcf} \\
        ${struc_vcf.size() == 0 ? 'UNSET' : struc_vcf} \\
        $alignments \\
        intersect_alignments.tsv \\
        ${ped.size() == 0 ? 'UNSET' : ped} \\
        intersect.ped
        
    touch warnings.txt
    """
}

