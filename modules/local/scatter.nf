
process SCATTER {
    label 'C4M4T1'
    label 'biopython'

    input:
    tuple path(vcf), path(tbi)
    val(n_shards)

    output:
    path("${prefix}*")

    script:
    prefix = vcf.name.replace('.vcf.gz', '').replace('.vcf.bgz', '') + '.shard'
    """
    scatter_vcf.py $vcf --n-shards $n_shards --output $prefix --threads $task.cpus
    """
}