
process vcf_split_sv_types {
    label 'C1M1T1'
    publishDir "progress/vcf_split_sv_types", mode: 'symlink'

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
