params.vep_output_opts = [
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
    '--variant_class'
].join(' ')

params.vep_filter_opts = [
    '--pick_allele_gene',
    '--no_intergenic'
].join(' ')


process vep {
    label 'C4M4T1'
    publishDir "progress/vep", mode: 'symlink'
    tag { "$set:$i:$j" }

    input:
        tuple val(set), val(i), val(j), path(vcf), path(fasta), path(fai), path(cache)

    output:
        tuple val(set), path(vep_vcf), path(mod_vcf), path(unann_vcf)

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


process vep_sv {
    label 'C4M4T2'
    publishDir "progress/vep_sv", mode: 'symlink'
    tag { "$set:$i:$j" }

    input:
        tuple val(set), val(i), val(j), path(vcf), path(pop_sv), path(pop_sv_tbi), path(fasta), path(fai), path(cache)

    output:
        tuple val(set), path(vep_vcf), path(unann_vcf)

    script:
    vep_vcf = vcf.name.replaceAll(/(\.vcf\.gz)|(\.bcf)$/, '.vep.bcf')
    unann_vcf = vcf.name.replaceAll(/(\.vcf\.gz)|(\.bcf)$/, '.unannotated.bcf')
    """
    mkfifo vep_out
    bcftools view --no-version  $vcf |
        vep --input_file STDIN \\
            $params.vep_output_opts \\
            $params.vep_filter_opts \\
            --fork 4 \\
            --format vcf \\
            --vcf \\
            --plugin StructuralVariantOverlap,reciprocal=1,file=$pop_sv,cols=AF:SVTYPE:END,label=SVO \\
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
            bcftools view --no-version -i "INFO/CSQ ~ '.'" -Ob -o $vep_vcf &
    
    bcftools view --no-version vep_out \\
        -e "INFO/CSQ ~ '.'" \\
        -Ob -o $unann_vcf

    wait && rm vep_out
    """
}

