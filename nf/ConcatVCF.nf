
workflow ConcatVCF {
    take: input_vcfs // set, vcf, index

    main:
    vcfs = input_vcfs |
        groupTuple(by: 0) |
        branch {
            multi: it[1].size() > 1
            single: true }

    output = vcfs.multi |
        concat_vcf |
        mix(vcfs.single.map { [it[0], it[1][0], it[2][0]] })

    emit:
    output
}

process concat_vcf {
    label 'C2M2T2'
    // publishDir "progress/ConcatVCF", mode: 'symlink'
    tag { "$set" }

    input:
    tuple val(set), path(vcfs), path(indices)

    output:
    tuple val(set), path(out_vcf), path("${out_vcf}.csi")

    script:
    out_vcf = "${set}.bcf"
    """
    bcftools concat $vcfs --allow-overlaps -Ob -o $out_vcf
    bcftools index --threads 2 $out_vcf
    """
}