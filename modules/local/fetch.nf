
process FETCH {
    label 'C1M1T1'
    label 'cavalier'
    publishDir "${params.outdir}/gene_lists"
    tag "$id:$ver"

    input:
    tuple val(id), val(ver)
    path(cav_opts)

    output:
    tuple val(id), path(output)

    script:
    output = id.replaceAll(':', '_') + '_' + ver + '.tsv'
    """
    R --slave --vanilla -e "\\
        cavalier::set_options_from_json('$cav_opts');\\
        list <- cavalier::get_gene_list('$id', version='$ver');\\
        readr::write_tsv(list, '$output')
    "
    """
}