
process VEP {
    label 'C2M2T2'
    label 'vep'
    tag "$i"

    input:
        tuple val(i),  path(vcf_input)
        tuple path(fasta), path(fai)
        path(vep_cache)
        tuple path(spliceai_snv), path(spliceai_snv_index), path(spliceai_indel), path(spliceai_indel_index)
        tuple path(alphamiss), path(alphamiss_index)
        tuple path(revel), path(revel_index)
        // tuple path(dbnfsp), path(dbnfsp_index)

    output:
        tuple path(output)

    script:
    output = vcf_input.name.replaceAll('.vcf.gz', '.vep.vcf.gz')
    """
    zcat $vcf_input \\
        | vep \\
            --input_file STDIN \\
            --format vcf \\
            --vcf \\
            --offline \\
            --dir_plugins /usr/local/share/ensembl-vep-114.0-0 \\
            --cache \\
            --cache_version $params.vep_cache_ver \\
            --dir $vep_cache \\
            --assembly GRCh38 \\
            --fasta $fasta \\
            --pick_allele_gene \\
            --hgvs \\
            --hgvsg \\
            --symbol \\
            --numbers \\
            --protein \\
            --no_stats \\
            --no_intergenic \\
            --variant_class \\
            --dont_skip \\
            --sift b \\
            --polyphen b \\
            --fields "${params.vep_fields.join(',')}" \\
            --plugin SpliceAI,snv=$spliceai_snv,indel=$spliceai_indel \\
            --plugin AlphaMissense,file=$alphamiss,transcript_match=1 \\
            --plugin REVEL,file=$revel \\
            --output_file STDOUT \\
        | bgzip -c > $output
    """
     // --plugin dbNSFP,$dbnfsp,transcript_match=1,Ensembl_transcriptid,${params.vep_dbnsfp_fileds} \\

    // NIN=\$(bcftools view --threads $task.cpus -H $vcf_input | wc -l)
    // NOUT=\$(bcftools view --threads $task.cpus -H $vep_bcf | wc -l)
    // if [[ "\$NIN" != "\$NOUT" ]]; then
    //     echo "Error: Number of input variants (\$NIN) not equal to number of output variants (\$NOUT)"
    //     exit 1
    // fi
}

