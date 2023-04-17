
workflow SplitSVs {
    take: vcfs // set, vcf, tbi

    main:
    all_sv_types = (params.sv_types.split(',') + params.sv_type_match.collectMany{ k,v -> v }).unique()

    if (params.sv_vcf) { // only run if sv vcf present
        output = vcfs |
            filter { it[0] == 'SV' } |
            map { it.drop(1) } |
            split_svs |
            flatMap { it[0] instanceof List ? it.transpose() : [it] } |
            map { [(it[0].name =~ /([A-Z]+)\.vcf\.gz$/)[0][1]] + it } |
            filter { all_sv_types.contains(it[0]) }

    } else {
        output = Channel.fromList([])
    }

    emit: output // set, vcf, tbi
}

process split_svs {
    label 'C1M1T1'
    // publishDir "progress/SplitSVs", mode: 'symlink'

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
