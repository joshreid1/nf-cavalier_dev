
process vcf_family_subset {
    label 'C2M2T2'
    publishDir "output/vcf_family_subset", mode: 'symlink'
    tag { fam }

    input:
        tuple path(vcf), val(fam), val(aff), val(unaff)

    output:
        tuple val(fam), path(out_vcf), path("${out_vcf}.tbi")

    script:
        out_vcf = "${fam}.subset.vcf.gz"
        """
        printf "${aff.join('\\n')}\\n" > aff
        bcftools view --no-update $vcf -Ou -s ${(aff + unaff).join(',')} |
            bcftools view  -i "GT[@aff]='alt'" -Oz -o $out_vcf
        bcftools index -t $out_vcf
        """
}
