
workflow GetSamples {
    take:
        vcfs // id, vcf, index

    main:
        output =
            proc(vcfs) |
            map { [it[0], it[1].toFile().readLines() as ArrayList] }

    emit:
        output // id, samples
}

process proc {
    label 'C1M1T1'
    publishDir "progress/GetSamples", mode: 'symlink'
    tag { id }

    input:
        tuple val(id), path(vcf), path(index)

    output:
        tuple val(id), path(out)

    script:
    out = "${params.id}.${id}.samples.txt"
    """
    bcftools query -l $vcf > $out
    """
}
