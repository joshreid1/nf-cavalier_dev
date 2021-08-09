

process vcf_split {
    cpus 2
    memory '2 GB'
    time '10 h'
    container 'bahlolab/mps-geno:latest'
    publishDir "output/vcf_split", mode: 'symlink'

    input:
    tuple val(id), file(vcf), file(tbi)

    output:
    file("$id-*.vcf.gz")

    script:
    """
    split_variants_chunked.py $vcf --n-chunk $params.n_split --out $id
    """
}
