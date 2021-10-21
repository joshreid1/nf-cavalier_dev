


process vep {
    label 'C2M8T4'
    publishDir "output/vep", mode: 'symlink'

    input:
        tuple val(id), path(vcf), path(fasta), path(fai), path(cache)

    output:
        tuple val(id), path(out_vcf), path("${out_vcf}.tbi")

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
        --hgvsg \\
        --fasta $fasta \\
        --assembly $params.vep_assembly \\
        --cache_version $params.vep_cache_ver \\
        --dir $cache \\
        --allow_non_variant \\
        --pick_allele_gene \\
        --output_file STDOUT | \\
        bcftools view --no-version -Oz -o $out_vcf
    bcftools index -t $out_vcf
    """
}
