
process vcf_flatten_multi {
    label 'C2M2T2'
    publishDir "output/vcf_flatten_multi", mode: 'symlink'

    input:
        tuple val(id), path(vcf)

    output:
        tuple val(id), path(out_vcf)

    script:
        out_vcf = "${id}.flat_multi.vcf.gz"
        """
        bcftools norm -m-any --do-not-normalize -Oz $vcf -o $out_vcf
        """
}
