

process vcfsplit {
    cpus 1
    memory '2 GB'
    time '1 h'
    container 'bahlolab/mps-geno:latest'
    publishDir "output/vcfsplit", mode: 'copy'

    input:
    tuple val(i), val(id), file(vcf)

    output:
    tuple val(out_id), file(out_vcf)

    script:
    out_id = "$id-$i"
    out_vcf = "${out_id}.vcf.gz"
    """
    split_variants_parallel.py $vcf $i $params.n_split --mode u  | \\
        bcftools view --no-version -Oz -o $out_vcf
    """
}
