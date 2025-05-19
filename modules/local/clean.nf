
process CLEAN {
    label 'C2M2T2'
    label 'bcftools'
    tag { i }

    input:
    tuple val(i), path(vcf)
    tuple path(ref), path(ref_fai)

    output:
    tuple val(i), path(output), path("${output}.tbi")

    script:
    output = vcf.name.replace('.vcf.gz', '.clean.vcf.gz')
    remove = 'QUAL,FILTER' +
        (params.snv_info    ? ',^' + params.snv_info  .collect{ "INF/$it" }.join(',') : '') +
        (params.snv_format  ? ',^' + params.snv_format.collect{ "FMT/$it" }.join(',') : '')
    """
    bcftools view $vcf --no-version -Ou -f "PASS,." \\
        | bcftools annotate --no-version -Ou --remove '$remove' \\
        | bcftools norm --no-version -Ou -m-any -f $ref \\
        ${params.fill_tags ? '| bcftools +fill-tags --no-version -Ou -- -t AC,AF,AN' : ''} \\
        | bcftools view -e 'AF=0 || ALT="*"' --no-version -Oz -o $output --write-index=tbi
    """
    // is ID necessary ? makes files larger
    //| bcftools annotate --no-version -Ou --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' \\
}