
include { spliceai_enabled  } from '../../functions/vep_helpers'
include { alphamiss_enabled } from '../../functions/vep_helpers'
include { revel_enabled     } from '../../functions/vep_helpers'
include { utr_ann_enabled   } from '../../functions/vep_helpers'
include { get_vep_fields    } from '../../functions/vep_helpers'

process VEP {
    /*
        Run VEP with various plugins
    */
    label 'C2M4T8'
    label 'vep'
    tag "$i"

    input:
        tuple val(i),  path(vcf_input)
        tuple path(fasta), path(fai)
        path(vep_cache)
        tuple path(spliceai_snv), path(spliceai_snv_index), path(spliceai_indel), path(spliceai_indel_index)
        tuple path(alphamiss), path(alphamiss_index)
        tuple path(revel), path(revel_index)
        path(utr_annotator)

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
            --fork $task.cpus \\
            --offline \\
            --dir_plugins /usr/local/share/ensembl-vep-114.0-0 \\
            --cache \\
            --cache_version $params.vep_cache_ver \\
            --dir $vep_cache \\
            --assembly GRCh38 \\
            --fasta $fasta \\
            --pick_allele_gene \\
            --check_existing \\
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
            --fields "${get_vep_fields().join(',')}" \\
            ${spliceai_enabled()  ? "--plugin SpliceAI,snv=$spliceai_snv,indel=$spliceai_indel" : ''} \\
            ${alphamiss_enabled() ? "--plugin AlphaMissense,file=$alphamiss,transcript_match=1" : ''} \\
            ${revel_enabled()     ? "--plugin REVEL,file=$revel" : ''} \\
            ${utr_ann_enabled()   ? "--plugin UTRAnnotator,file=$utr_annotator" : ''} \\
            --output_file STDOUT \\
        | bgzip -c > $output
    """
    // --- FORMER BASH FOR CHECKING NO VARS DROPPED ---
    // NIN=\$(bcftools view --threads $task.cpus -H $vcf_input | wc -l)
    // NOUT=\$(bcftools view --threads $task.cpus -H $vep_bcf | wc -l)
    // if [[ "\$NIN" != "\$NOUT" ]]; then
    //     echo "Error: Number of input variants (\$NIN) not equal to number of output variants (\$NOUT)"
    //     exit 1
    // fi
}

