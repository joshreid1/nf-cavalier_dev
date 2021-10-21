

process vcf_split_norm {
    label 'C2M2T8'
    publishDir "output/vcf_split_norm", mode: 'symlink'

    input:
        tuple path(bed), path(vcf), path(tbi), path(ref), path(ref_fai)

    output:
        path(out)

    script:
        out = bed.name.replaceAll('.bed', '.norm.vcf.gz')
        """
        bcftools view --no-version $vcf -R $bed -Ou |
            bcftools norm --no-version -m-any -f $ref -Ou |
            bcftools annotate --no-version --threads 2 -Oz \\
                --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' \\
                --output $out
        """
}
