
process vcf_family_subset {
    cpus 2
    memory '1 GB'
    time '1 h'
    container 'jemunro/nf-long-amplicon-typing:dev'
    publishDir "output/vcf_family_subset", mode: 'copy'

    input:
    tuple file(vcf), val(samples)

    output:
    tuple val(sample), file(out_vcf)

    script:
    out_vcf = "${samples[0]}.subset.vcf.gz"
    sample = samples[0]
    """
    bcftools view --no-version $vcf -Ou -s ${samples.join(',')} | 
        bcftools view --no-version -i "GT[0]='alt'" -Oz -o $out_vcf
    """
}
