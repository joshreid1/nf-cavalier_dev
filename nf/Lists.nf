
include { list_channels; path; date_ymd; get_options_json } from './functions'

workflow Lists {

    main:
        
        list_channel = channel.fromList([])

        (web_lists, local_lists) = list_channels()
        
        if (web_lists != null) {
            list_channel =  web_lists |
                unique |
                collectFile(name: 'list_ids.txt', newLine: true) |
                combine([date_ymd()]) |
                update_versions |
                splitCsv(sep: '\t', skip: 1, strip: true) |
                pull_latest |
                map { it[1] }
        }
        if (local_lists != null) {
            list_channel = 
                local_lists | 
                map { path(it) } |
                map_list_ids |
                mix(list_channel)
        }

        list_channel = list_channel.collect()

    emit:
        list_channel
}

process map_list_ids {
    label 'C1M1T1'
    label 'cavalier'
    tag "${input.name}"
    publishDir "${params.outdir}/gene_lists"

    input: path(input)

    output: path(output)

    script:
    output = input.name.replaceAll('\\.tsv$', '') + '.mapped.tsv'
    """
    cavalier_map_list_ids.R $input $output '${get_options_json()}'
    """
}

process update_versions {
    label 'C1M1T1'
    label 'cavalier'
    tag { date }

    input:
    tuple path(id_file), val(date)

    output:
    path(output)

    script:
    output = "list_versions_${date}.tsv"
    cmd = "cavalier::get_gene_list_versions(readr::read_lines('$id_file'), '$output')"
    """
    R --slave --vanilla -e "\\
        cavalier::set_options_from_json('${get_options_json('\\\\\\"')}');\\
        cavalier::get_gene_list_versions(readr::read_lines('$id_file'), '$output')\\
    "
    """
}

process pull_latest {
    label 'C1M1T1'
    label 'cavalier'
    storeDir "${params.outdir}/gene_lists"
    tag "$id:$ver"

    input:
    tuple val(id), val(ver)

    output:
    tuple val(id), path(output)

    script:
    output = id.replaceAll(':', '_') + '_' + ver + '.tsv'
    """
    R --slave --vanilla -e "\\
        cavalier::set_options_from_json('${get_options_json('\\\\\\"')}');\\
        cavalier::get_gene_list('$id', version='$ver', save='$output')\\
    "
    """
}

