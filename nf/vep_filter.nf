

process vep_filter {
    label 'C2M8T4'
    publishDir "output/vep_filter", mode: 'symlink'

    input:
    tuple val(id), file(vcf), file(tbi)

    output:
    tuple val(id), file(out_vcf), file("${out_vcf}.tbi")

    script:
    out_vcf = "${id}.vep-filtered.vcf.gz"
    """
    filter_vep --input_file $vcf  \\
        --format vcf \\
        --only_matched \\
        --filter "MAX_AF < $params.max_af or not MAX_AF" \\
        --filter "IMPACT in ${params.vep_impact.join(',')}" |  \\
        bcftools view --no-version -Oz -o $out_vcf
    bcftools index -t $out_vcf
    """
}
