
workflow CountVCF {
    take: vcfs // id, vcf, tbi

    main:
    counts = vcfs |
        count_vcf |
        map { [it[0], it[1].toFile().text as int] }
    output = vcfs.combine(counts, by:0)

    emit: output // id, vcf, tbi, count
}

process count_vcf {
    label 'C1M1T1'
    // publishDir "progress/CountVCF", mode: 'symlink'
    tag { id }

    input:
        tuple val(id), path(vcf), path(tbi)

    output:
        tuple val(id), path(out)

    script:
    out = "${id}.nrecords"
    """
    bcftools index --nrecords $vcf > $out
    """
}
