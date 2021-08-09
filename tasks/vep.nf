
params.vep_cache = '/stornext/Bioinf/data/lab_bahlo/ref_db/vep-cache'
params.vep_cache_ver = '104'
params.vep_assembly = 'GRCh37'

process vep {
    cpus 2
    memory '8 GB'
    time '4 h'
    container 'jemunro/nf-long-amplicon-typing:dev'
    publishDir "output/vep", mode: 'symlink'

    input:
    tuple val(id), file(vcf)

    output:
    tuple val(id), file(out_vcf), file("${out_vcf}.tbi")

    script:
    out_vcf = "${id}.vep.vcf.gz"
    """
    vep --input_file $vcf \\
        --format vcf \\
        --vcf \\
        --cache \\
        --offline \\
        --everything \\
        --max_af \\
        --allele_number \\
        --variant_class \\
        --dont_skip \\
        --assembly $params.vep_assembly \\
        --cache_version $params.vep_cache_ver \\
        --dir $params.vep_cache \\
        --allow_non_variant \\
        --pick_allele_gene \\
        --output_file STDOUT | \\
        bcftools view --no-version -Oz -o $out_vcf
    bcftools index -t $out_vcf
    """
}
