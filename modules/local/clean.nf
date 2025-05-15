
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
    """
    bcftools view $vcf --no-version -Ou -f "PASS,." \\
        | bcftools annotate --no-version -Ou --remove QUAL,FILTER,${params.snp_format_keep ? '^' + params.snp_format_keep :'FORMAT'} \\
        | bcftools annotate --no-version -Ou --remove ${params.snp_info_keep   ? '^' + params.snp_info_keep   : 'INFO'} \\
        | bcftools norm --no-version -Ou -m-any -f $ref \\
        ${params.fill_tags ? '| bcftools +fill-tags --no-version -Ou -- -t AC,AF,AN' : ''} \\
        | bcftools view -e 'AF=0 || ALT="*"' --no-version -Oz -o $output --write-index=tbi
    """
    // is ID necessary ? makes files larger
    //| bcftools annotate --no-version -Ou --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' \\
}