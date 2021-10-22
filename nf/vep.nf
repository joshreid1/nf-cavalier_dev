
process vep {
    label 'C2M4T2'
    publishDir "output/vep", mode: 'symlink'

    input:
        tuple path(vcf), path(fasta), path(fai), path(cache)

    output:
        path(out_vcf)

    script:
    out_vcf = vcf.name.replaceAll('.bcf', '.vep.bcf')
    vep_output_opts = [
        '--sift b',
        '--polyphen b',
        '--ccds',
        '--hgvs',
        '--hgvsg',
        '--symbol',
        '--numbers',
        '--protein',
        '--af',
        '--af_1kg',
        '--af_gnomad',
        '--max_af',
        '--variant_class',
//        '--mane'
//        '--var_synonyms',
//        '--pubmed',
//        '--af_esp',
//        '--gene_phenotype',
//        '--appris',
//        '--tsl',
//        '--uniprot',
//        '--biotype',
//        '--canonical',
//        '--regulatory',
//        '--domains',
    ].join(' ')

    vep_filter_opts = [
        '--pick_allele_gene',
        '--no_intergenic'
//        '--allow_non_variant',
//        '--dont_skip',
    ].join(' ')
    """
    bcftools view --no-version  $vcf |
        vep --input_file STDIN \\
            $vep_output_opts \\
            $vep_filter_opts \\
            --fork 2 \\
            --format vcf \\
            --vcf \\
            --cache \\
            --offline \\
            --no_stats \\
            --fasta $fasta \\
            --assembly $params.vep_assembly \\
            --cache_version $params.vep_cache_ver \\
            --dir $cache \\
            --output_file STDOUT | \\
        bcftools view --threads 2 --no-version -Ob -o $out_vcf
    """
    //    bcftools index --threads 2 -t $out_vcf
}
