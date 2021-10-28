
process vcf_count {
    label 'C1M1T1'
    publishDir "progress/vcf_sample_set", mode: 'symlink'

    input:
        tuple path(vcf), path(tbi)

    output:
        path('n.txt')

    script:
    """
    bcftools index --nrecords $vcf > n.txt
    """
}
