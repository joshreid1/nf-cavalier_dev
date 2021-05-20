
params.max_af = 0.10
params.vep_impact = ['MODERATE', 'HIGH']

process vep_filter {
    cpus 2
    memory '8 GB'
    time '4 h'
    container 'jemunro/nf-long-amplicon-typing:dev'
    publishDir "output/vep_filter", mode: 'copy'

    input:
    tuple val(id), file(vcf), file(tbi)

    output:
    tuple val(id), file(out_vcf), file("${out_vcf}.tbi")

    script:
    out_vcf = "${id}.vep-filtered.vcf.gz"
    """
    filter_vep --input_file $vcf  \\
        --format vcf \\
        --vcf_info_field ANN \\
        --only_matched \\
        --filter "MAX_AF < $params.max_af or not MAX_AF" \\
        --filter "IMPACT in ${params.vep_impact.join(',')}" |  \\
        bcftools view --no-version -Oz -o $out_vcf
    bcftools index -t $out_vcf
    """
}
