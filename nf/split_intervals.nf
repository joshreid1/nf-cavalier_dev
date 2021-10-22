

process split_intervals {
    label 'C2M2T8'
    label 'cavalier'
    publishDir "output/split_intervals", mode: 'symlink'

    input:
        tuple val(n), path(ref_fai), path(gaps)

    output:
        path("split-*.intervals.tsv")

    script:
    """
    split_intervals.R $ref_fai $gaps $n split
    """
}
