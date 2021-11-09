
params.naive = true

process vcf_merge {
    label 'C2M2T2'
    publishDir "output/vcf_merge", mode: 'copy'
    tag { set }

    input:
    tuple val(set), path(file_list)

    output:
    tuple val(set), path(out_vcf), path("${out_vcf}.csi")

    script:
    out_vcf = "${params.id}.${set}.bcf"
    """
    bcftools concat ${params.naive ? '--naive-force' : '-a' } \\
        --file-list $file_list -Ob -o $out_vcf
    bcftools index --threads 2 $out_vcf
    """
}
