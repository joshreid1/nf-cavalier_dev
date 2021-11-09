
process vep_plugin_svo {
    label 'C2M4T2'
    publishDir "progress/vep_plugin_svo", mode: 'symlink'
    tag { type }

    input:
        tuple val(type), path(vcf), path(tbi), path(pop_sv), path(pop_sv_tbi), path(fasta), path(fai), path(cache)

    output:
        tuple val(type), path(out_vcf), path("${out_vcf}.tbi")

    script:
    out_vcf = vcf.name.replaceAll('.vcf.gz', '.svo.vcf.gz')
    merge_cmd = pop_sv.size() > 1 ?
        "bcftools concat -a ${pop_sv.join(' ')} -Oz -o merged.vcf.gz && bcftools index -t merged.vcf.gz" :
        "mv ${pop_sv[0]} merged.vcf.gz && mv ${pop_sv[0]} merged.vcf.gz.tbi"
    """
    $merge_cmd
    bcftools view --no-version  $vcf |
        vep --input_file STDIN \\
            --format vcf \\
            --vcf \\
            --plugin StructuralVariantOverlap,reciprocal=1,file=meged.vcf.gz,cols=AF:SVTYPE:END,label=SVO \\
            --cache \\
            --offline \\
            --no_stats \\
            --fasta $fasta \\
            --assembly $params.vep_assembly \\
            --cache_version $params.vep_cache_ver \\
            --dir $cache \\
            --output_file STDOUT |
            bcftools view --no-version -Ob -o $out_vcf 
    bcftools index -t $out_vcf
    """
}
