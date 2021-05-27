
process vcf_flatten_multi {
    cpus 1
    memory '1 GB'
    time '1 h'
    container 'jemunro/nf-long-amplicon-typing:dev'
    publishDir "output/vcf_flatten_multi", mode: 'copy'

    input:
    tuple val(id), file(vcf)

    output:
    tuple val(id), file(out_vcf)

    script:
    out_vcf = "${id}.flat_multi.vcf.gz"
    """
    bcftools norm -m-any --do-not-normalize -Oz $vcf -o $out_vcf
    """
}
