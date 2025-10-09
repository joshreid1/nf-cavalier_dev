
process FILTER {
    label 'C2M2T2'
    label 'bcftools'
    tag "$i"
    /*
        - Apply a filter to a VCF file
        - Used primarily to filter based on vcfanno annotation (e.g. population frequency)
    */

    input:
    tuple val(i), path(vcf)
    val(include_exp)

    output:
    tuple val(i), path(output)

    script:
    output = vcf.name.replace('.vcf.gz', '.flt.vcf.gz')
    """
    bcftools view $vcf --no-version --threads ${task.cpus} -i '$include_exp' -Oz -o $output
    """
}