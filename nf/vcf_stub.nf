
process vcf_stub {
    label 'C1T1M1'
    publishDir "progress/vcf_stub", mode: 'symlink'

    input:
        path(vcf)

    output:
        tuple path(stub), path("${stub}.tbi")

    script:
    stub = vcf.name.replaceAll(/(\.vcf\.gz)|(\.bcf)$/, '.stub.vcf.gz')
    """
    bcftools view --no-version -h $vcf | 
        bcftools view --no-version -Oz -o $stub
    bcftools index -t $stub
    """
}
