
process SAMPLES {
    label 'C2M2T2'
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
