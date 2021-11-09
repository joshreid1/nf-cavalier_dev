

process vcf_split_norm {
    label 'C2M2T2'
    publishDir "progress/vcf_split_norm", mode: 'symlink'

    input:
        tuple path(regions), path(vcf), path(tbi), path(ref), path(ref_fai)

    output:
        path("$pref*.bcf")

    script:
        pref = regions.name.replaceAll('.intervals.tsv', '.norm')
        if (params.mode == 'short')
            """
            bcftools view --no-version $vcf -R $regions -Ou |
                bcftools norm --no-version -m-any -f $ref -Ou |
                bcftools annotate --no-version --threads 2 -Ou \\
                    --remove INFO/CSQ \\
                    --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' | \\
                 split_variants.py --bcf \\
                    --chunk-size ${params.chunk_size} \\
                    --pref $pref
            """
        else if (params.mode == 'sv')
            """
            bcftools view --no-version $vcf -R $regions -Ou |
                bcftools annotate --no-version --threads 2 -Ou \\
                    --remove INFO/CSQ | \\
                 split_variants.py --bcf \\
                    --chunk-size ${params.chunk_size} \\
                    --pref $pref
            """
}
