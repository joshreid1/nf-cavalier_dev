
process GATHER {
    label 'C4M4T2'
    label 'bcftools'
    /*
        - Gather scattered annotated VCF filters
        - We can use a naive merge (much faster) since records are already
    */

    input:
    path(vcfs)

    output:
    tuple path(output), path("${output}.tbi")

    script:
    output = vcfs[0].name.replaceAll('shard\\.\\d+\\.', '')
    """
    bcftools concat $vcfs --threads $task.cpus --naive -Oz -o $output
    bcftools index -t --threads $task.cpus $output
    """
}