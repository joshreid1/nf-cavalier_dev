
params.vep_cache = '/stornext/Bioinf/data/lab_bahlo/ref_db/vep-cache'

process vep {
    cpus 2
    memory '8 GB'
    time '4 h'
    container 'jemunro/nf-long-amplicon-typing:dev'
    publishDir "output/vep", mode: 'copy'

    input:
    tuple val(id), file(vcf)

    output:
    tuple val(id), file(out_vcf), file("${out_vcf}.tbi")

    script:
    out_vcf = "${id}.vep.vcf.gz"
    """
    vep --input_file $vcf  \\
        --format vcf \\
        --vcf \\
        --cache \\
        --offline \\
        --everything \\
        --max_af \\
        --allele_number \\
        --variant_class \\
        --dont_skip \\
        --assembly GRCh37 \\
        --cache_version 101 \\
        --dir $params.vep_cache \\
        --allow_non_variant \\
        --flag_pick_allele_gene \\
        --output_file STDOUT | \\
        bcftools view --no-version -Oz -o $out_vcf
    bcftools index -t $out_vcf
    """
    // added --pick option to VEP, this will pick single best annotation per variant
    // this makes it simpler for cavalier
}
