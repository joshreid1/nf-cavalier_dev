
params.allow_overlap = false


process vcf_concat {
    label 'C2M2T2'
    publishDir "${params.outdir}/vcf_concat", mode: 'copy'
    tag { "$set:$type" }

    input:
    tuple val(set), val(type), path(file_list)

    output:
    tuple val(set), val(type), path(out_vcf), path("${out_vcf}.csi")

    script:
    out_vcf = "${set}.${type}.bcf"
    """
    bcftools concat ${params.allow_overlap ? '--allow-overlaps' : '--naive-force' } \\
        --file-list $file_list -Ob -o $out_vcf
    bcftools index --threads 2 $out_vcf
    """
}
