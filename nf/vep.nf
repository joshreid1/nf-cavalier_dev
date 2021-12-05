params.vep_output_opts =
    [
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

params.vep_filter_opts =
    [
        '--pick_allele_gene',
        '--no_intergenic'
//        '--allow_non_variant',
//        '--dont_skip',
    ].join(' ')


process vep {
    label 'C4M4T1'
    publishDir "progress/vep", mode: 'symlink'

    input:
        tuple path(vcf), path(fasta), path(fai), path(cache)

    output:
        tuple path(vep_vcf), path(mod_vcf), path(unann_vcf)

    script:
        vep_vcf = vcf.name.replaceAll('.bcf', '.vep.bcf')
        mod_vcf = vcf.name.replaceAll('.bcf', '.vep-modifier.bcf')
        unann_vcf = vcf.name.replaceAll('.bcf', '.unannotated.bcf')
    """
    mkfifo vep_out
    mkfifo filter_in
    bcftools view --no-version  $vcf |
        vep --input_file STDIN \\
            $params.vep_output_opts \\
            $params.vep_filter_opts \\
            --fork 4 \\
            --format vcf \\
            --vcf \\
            --cache \\
            --offline \\
            --no_stats \\
            --fasta $fasta \\
            --assembly $params.vep_assembly \\
            --cache_version $params.vep_cache_ver \\
            --dir $cache \\
            --output_file STDOUT |
            bcftools view --no-version -Ou |
            tee vep_out | \\
            bcftools view --no-version -i "INFO/CSQ ~ '.'" -Ov |
            tee filter_in |
            filter_vep \\
                --format vcf \\
                --only_matched \\
                --filter "IMPACT in LOW,MODERATE,HIGH" |
            bcftools view --no-version -Ob -o $vep_vcf &
    
    bcftools view --no-version vep_out \\
        -e "INFO/CSQ ~ '.'" \\
        -Ob -o $unann_vcf &

    cat filter_in |
        filter_vep \\
        --format vcf \\
        --only_matched \\
        --filter "IMPACT is MODIFIER" |
        bcftools view --no-version -Ob -o $mod_vcf

    wait && rm vep_out filter_in
    """
}


process vep_svo {
    label 'C4M4T2'
    publishDir "progress/vep_svo", mode: 'symlink'
    tag { type }

    input:
        tuple val(type), path(vcf), path(tbi), path(pop_sv), path(pop_sv_tbi), path(fasta), path(fai), path(cache)

    output:
        tuple val(type), path(vep_vcf), path(unann_vcf)

    script:
    vep_vcf = vcf.name.replaceAll(/(\.vcf\.gz)|(\.bcf)$/, '.vep.bcf')
    unann_vcf = vcf.name.replaceAll(/(\.vcf\.gz)|(\.bcf)$/, '.unannotated.bcf')
    merge_cmd = pop_sv.size() > 1 ?
        "bcftools concat -a ${pop_sv.join(' ')} -Oz -o popsv.vcf.gz && bcftools index -t popsv.vcf.gz" :
        "mv ${pop_sv[0]} popsv.vcf.gz && mv ${pop_sv[0]} popsv.vcf.gz.tbi"
    """
    $merge_cmd
    mkfifo vep_out
    bcftools view --no-version  $vcf |
        vep --input_file STDIN \\
            $params.vep_output_opts \\
            $params.vep_filter_opts \\
            --fork 4 \\
            --format vcf \\
            --vcf \\
            --plugin StructuralVariantOverlap,reciprocal=1,file=popsv.vcf.gz,cols=AF:SVTYPE:END,label=SVO \\
            --cache \\
            --offline \\
            --no_stats \\
            --fasta $fasta \\
            --assembly $params.vep_assembly \\
            --cache_version $params.vep_cache_ver \\
            --dir $cache \\
            --output_file STDOUT |
            bcftools view --no-version -Ou |
            tee vep_out | \\
            bcftools view --no-version -i "INFO/CSQ ~ '.'" -Ob -o $vep_vcf &&
            bcftools index $vep_vcf &
    
    bcftools view --no-version vep_out \\
        -e "INFO/CSQ ~ '.'" \\
        -Ob -o $unann_vcf &&
        bcftools index $unann_vcf &

    wait && rm vep_out
    """
}

