
process GATHER {
    label 'C4M4T2'
    label 'bcftools'

    input:
    path(vcfs)

    output:
    tuple path(output), path("${output}.tbi")

    script:
    // println vcfs[0]
    output = vcfs[0].name.replaceAll('shard\\.\\d+\\.', '')
    """
    bcftools concat $vcfs --threads $task.cpus --naive -Oz -o $output
    bcftools index -t  --threads $task.cpus $output
    """
}