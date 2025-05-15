
process FILTER {
    label 'C2M2T2'
    label 'bcftools'
    tag { i }

    input:
    tuple val(i), path(vcf)
    val(include)

    output:
    tuple val(i), path(output)

    script:
    output = vcf.name.replace('.vcf.gz', '.flt.vcf.gz')
    """
    bcftools view $vcf --no-version --threads ${task.cpus} -i '$include' -Oz -o $output
    """
}