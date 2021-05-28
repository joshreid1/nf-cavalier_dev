
process vcf_merge {
    cpus 2
    memory '4 GB'
    time '2 h'
    container 'jemunro/nf-long-amplicon-typing:dev'
    publishDir "output/vcf_merge", mode: 'copy'

    input:
    tuple val(id), file(vcfs), file(tbis)

    output:
    tuple val(id), file(out_vcf), file("${out_vcf}.tbi")

    script:
    out_vcf = "${id}.vep-filtered.merged.vcf.gz"
    """
    bcftools concat -a --no-version --threads 2 -Oz ${vcfs.join(' ')} > $out_vcf
    bcftools index -t $out_vcf
    """
}
