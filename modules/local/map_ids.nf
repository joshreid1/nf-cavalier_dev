
process MAP_IDS {
    label 'C1M1T1'
    label 'cavalier'
    tag "${input.name}"
    publishDir "${params.outdir}/gene_lists"
    /*
        - Map between various gene ids with cavalier
    */

    input: 
    path(input)
    path(cav_opts)

    output: path(output)

    script:
    output = input.name.replaceAll('\\.tsv$', '') + '.mapped.tsv'
    """
    cavalier_map_list_ids.R $input $output $cav_opts
    """
}