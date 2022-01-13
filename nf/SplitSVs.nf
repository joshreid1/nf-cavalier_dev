
workflow SplitSVs {
    take: vcfs // set, vcf, tbi

    main:
    if (params.sv_vcf) { // only run if sv vcf present
        output = vcfs |
            filter { it[0] == 'SV' } |
            map { it.drop(1) } |
            proc |
            flatMap { it.transpose() } |
            map { [(it[0].name =~ /([A-Z]+)\.vcf\.gz$/)[0][1]] + it } |
            filter { params.sv_types.contains(it[0]) } |
            mix(vcfs.filter { it[0] == 'SNP'} )

    } else {
        output = vcfs
    }

    emit: output // set, vcf, tbi
}

process proc {
    label 'C1M1T1'
    publishDir "progress/SplitSVs", mode: 'symlink'

    input:
        tuple path(vcf), path(index)

    output:
        tuple path("$pref*.vcf.gz"), path("$pref*.vcf.gz.tbi")

    script:
        pref = vcf.name.replaceAll(/(\.vcf\.gz)|(\.bcf)$/, '')
        """
        bcftools view --no-version $vcf |
             split_sv_types.py --input - \\
                --index \\
                --pref $pref
        """
}
