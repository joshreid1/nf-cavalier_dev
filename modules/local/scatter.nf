
process SCATTER {
    label 'C4M4T2'
    label 'biopython'
    /*
        - Split VCF into pieces by using TBI index
        - requires BGZIPPED VCF
    */

    input:
    tuple path(vcf), path(tbi)
    val(n_shards)
    val(check)

    output:
    path("${prefix}*")

    script:
    prefix = vcf.name.replace('.vcf.gz', '').replace('.vcf.bgz', '') + '.shard'
    """
    scatter_vcf.py $vcf --n-shards $n_shards --output $prefix --threads $task.cpus
    """
}