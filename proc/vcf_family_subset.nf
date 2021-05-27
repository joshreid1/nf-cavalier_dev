
params.max_cohort_ac = 1

process vcf_family_subset {
    cpus 2
    memory '1 GB'
    time '1 h'
    container 'jemunro/nf-long-amplicon-typing:dev'
    publishDir "output/vcf_family_subset", mode: 'copy'

    input:
    tuple file(vcf), val(samples)

    output:
    tuple val(samples), file(out_vcf), file("${out_vcf}.tbi")

    script:
    out_vcf = "${samples[0]}.subset.vcf.gz"
    """
    bcftools view --no-version $vcf -Ou -s ${samples.join(',')} | 
        bcftools view --no-version -i "GT[0]='alt' & AC<=$params.max_cohort_ac" -Oz -o $out_vcf
    bcftools index -t --threads 2 $out_vcf
    """
}
