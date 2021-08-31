

process vcf_split {
    label 'C2M2T8'
    publishDir "output/vcf_split", mode: 'symlink'

    input:
        tuple path(vcf), path(tbi)

    output:
        path("$params.id-*.vcf.gz")

    script:
        """
        split_variants_chunked.py $vcf --n-chunk $params.n_split --out $params.id
        """
}
