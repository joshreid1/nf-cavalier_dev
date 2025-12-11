
process CLEAN {
    label 'C2M2T2'
    label 'bcftools'
    tag { i }
    /*
        Clean input VCF:
        - Filter for variants based on FILTER field
        - Drop unused annotations
        - Split multiallelics to single records and normalise representation
        - Fill AC/AF/AN if enabled
        - Drop non-variants (AF ==0) or variants with alt == '*' (GATK emits these)
    */

    input:
    tuple val(i), path(vcf)
    tuple path(ref), path(ref_fai)

    output:
    tuple val(i), path(output), path("${output}.tbi")

    script:
    output = vcf.name.replace('.vcf.gz', '.clean.vcf.gz')
    // decrease VCF size by dropping unused fields
    remove = 'QUAL,FILTER' +
        (params.short_info    ? ',^' + params.short_info  .collect{ "INF/$it" }.join(',') : '') +
        (params.short_format  ? ',^' + params.short_format.collect{ "FMT/$it" }.join(',') : '')
    """
    bcftools view $vcf --no-version -Ou ${params.short_vcf_filter ? "-f  '$params.short_vcf_filter'": ''} \\
        | bcftools annotate --no-version -Ou --remove '$remove' \\
        | bcftools norm --no-version -Ou -m-any -f $ref \\
        ${params.short_fill_tags ? '| bcftools +fill-tags --no-version -Ou -- -t AC,AF,AN' : ''} \\
        | bcftools view -e 'AF=0 || ALT="*"' --no-version -Oz -o $output --write-index=tbi
    """
}