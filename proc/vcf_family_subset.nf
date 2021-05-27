
params.max_cohort_ac = 1

process vcf_sample_subset {
    cpus 2
    memory '1 GB'
    time '1 h'
    container 'jemunro/nf-long-amplicon-typing:dev'
    publishDir "output/vcf_sample_subset", mode: 'copy'

    input:
    tuple file(vcf), val(sample)

    output:
    tuple val(sample), file(out_vcf), file("${out_vcf}.tbi")

    script:
    out_vcf = "${sample}.subset.vcf.gz"
    """
    bcftools view --no-version $vcf -Ou -s $sample | 
        bcftools view --no-version -i "GT='alt' & AC<=$params.max_cohort_ac" -Oz -o $out_vcf
    bcftools index -t --threads 2 $out_vcf
    """
}
