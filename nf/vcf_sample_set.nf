
process vcf_sample_set {
    label 'C1M1T1'
    publishDir "progress/vcf_sample_set", mode: 'symlink'

    input:
        path(vcf)

    output:
        path(out)

    script:
    out = "${params.id}.samples.txt"
    """
    bcftools query -l $vcf > $out
    """
}
