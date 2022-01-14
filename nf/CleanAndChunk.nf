
workflow CleanAndChunk {
    take:
        vcf_counts // set, vcf, tbi, count
        ref_data // ref_fa, ref_fai, gaps, vep_cache

    main:

    output = vcf_counts |
        map { [it[0], Math.ceil(
            it[3] / (it[0] == 'SNP' ? params.chunk_size : params.sv_chunk_size )) as int ] } |
        combine(ref_data.map { it[1..2] } ) | // set, n_chunk, ref_fai, gaps
        split_intervals |
        flatMap { set, intvls ->
            intvls instanceof List ?
                intvls.collect{ [set, (it.name =~ /(?<=split-)([0-9]+)/)[0][1], it] } :
                [[set, (intvls.name =~ /(?<=split-)([0-9]+)/)[0][1], intvls]]
        } |
        combine(vcf_counts.map { it[0..2] }, by: 0) |
        combine(ref_data.map { it[0..1] }) | // set, i, intvls, vcf, tbi, ref_fa, ref_fai
        clean_and_chunk |
        flatMap { set, i,  vcf, index ->
            vcf instanceof List ?
                [vcf, index].transpose().collect { [set, i, (it[0].name =~ /(?<=clean-)([0-9]+)/)[0][1]] + it } :
                [[set, i, (vcf.name =~ /(?<=clean-)([0-9]+)/)[0][1], vcf, index]]
        } |
        map { [it[0], it[1] as int, it[2] as int, it[3], it[4] ] }

    emit: output // set, i, j, vcf, tbi
}

process split_intervals {
    label 'C2M2T8'
    label 'cavalier'
    publishDir "progress/split_intervals", mode: 'symlink'
    tag { set }

    input:
    tuple val(set), val(n), path(ref_fai), path(gaps)

    output:
    tuple val(set), path("$set-split-*.intervals.tsv")

    script:
    """
    split_intervals.R $ref_fai $gaps $n $set-split
    """
}

process clean_and_chunk {
    label 'C2M2T2'
    publishDir "progress/clean_and_chunk", mode: 'symlink'
    tag { "$set:$i" }

    input:
    tuple val(set), val(i), path(regions), path(vcf), path(tbi), path(ref), path(ref_fai)

    output:
    tuple val(set), val(i), path("$pref*.bcf"), path("$pref*.bcf.csi")

    script:
    pref = regions.name.replaceAll('.intervals.tsv', '.clean')
    if (set == 'SNP')
        """
        bcftools view --no-version $vcf -R $regions -Ou |
            bcftools norm --no-version -m-any -f $ref -Ou |
            bcftools annotate --no-version --threads 2 -Ou \\
                --remove INFO/CSQ \\
                --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' | \\
             split_variants.py --bcf --index \\
                --chunk-size ${params.chunk_size} \\
                --pref $pref
        """
    else
        """
        bcftools view --no-version $vcf -R $regions -Ou |
            bcftools annotate --no-version --threads 2 -Ou \\
                --remove INFO/CSQ | \\
             split_variants.py --bcf --index \\
                --chunk-size ${params.sv_chunk_size} \\
                --pref $pref
        """
}
