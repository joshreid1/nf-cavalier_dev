
include { read_lists; path; date_ymd; get_options_json } from './functions'

workflow Lists {

    main:
        lists = read_lists()

        if (lists.any { it.list  ==~ '^(HP|PA[AE]):.+'}) {
            lists = Channel.from(lists) |
                map { it.values() as ArrayList } |
                branch {
                    web: it[1]  ==~ '^(HP|PA.+|G4E):.+'
                    file: true
                }

            list_channel =
                lists.web |
                    map { it[1] } |
                    unique |
                    collectFile(name: 'list_ids.txt', newLine: true) |
                    combine([date_ymd()]) |
                    update_versions |
                    splitCsv(sep: '\t', skip: 1, strip: true) |
                    pull_latest |
                    combine(lists.web.map {it[[1,0]]}, by:0) |
                    map { [it[2], it[1]] } |
                    mix(lists.file.map { [it[0], path(it[1])] }) |
                    groupTuple(by: 0)
        } else {
            list_channel = Channel.from(lists) |
                map { it.values() as ArrayList } |
                map { [it[0], path(it[1])] } |
                groupTuple(by: 0)
        }
    emit:
        list_channel // fam, lists
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

