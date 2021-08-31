
process vcf_merge {
    label 'C2M2T2'
    publishDir "output/vcf_merge", mode: 'symlink'

    input:
    tuple file(vcfs), file(tbis)

    output:
    tuple file(out_vcf), file("${out_vcf}.tbi")

    script:
    out_vcf = "${params.id}.vep-filtered.merged.vcf.gz"
    """
    bcftools concat -a --no-version --threads $task.cpus -Oz ${vcfs.join(' ')} > $out_vcf
    bcftools index -t $out_vcf
    """
}
