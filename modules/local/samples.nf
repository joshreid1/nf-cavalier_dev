
process SAMPLES {
    label 'C1M1T1'
    label 'bcftools'
    /*
        - return list of sample ids in VCF
    */

    input:
        tuple path(vcf), path(index)

    output:
        path(samples)

    script:
    samples = "samples.txt"
    """
    bcftools query -l $vcf > $samples
    """
}
