
workflow CleanAndChunkVCF {
    take:
        vcfs_count // id, vcf, tbi, count
        ref_data // ref_fa, ref_fai, gaps, vep_cache

    main:
    counts = vcfs |
        proc |
        map { [it[0], it[1].toFile().text as int] }
    output = vcfs.combine(counts, by:0)

    emit: output // id, vcf, tbi, count
}

process split_intervals {
    label 'C2M2T8'
    label 'cavalier'
    publishDir "progress/split_intervals", mode: 'symlink'

    input:
    tuple val(n), path(ref_fai), path(gaps)

    output:
    path("split-*.intervals.tsv")

    script:
    """
    split_intervals.R $ref_fai $gaps $n split
    """
}