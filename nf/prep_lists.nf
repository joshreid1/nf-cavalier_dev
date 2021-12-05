
include { path; date_ymd } from './functions'

workflow prep_lists {
    take:
        lists

    main:
        if (lists.any { it.list  ==~ '^(HP|PA[AE]):.+'}) {
            lists = Channel.from(lists) |
                map { it.values() as ArrayList } |
                branch {
                    web: it[1]  ==~ '^(HP|PA[AE]):.+'
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
        list_channel
}

process pull_latest {
    label 'C1M1T1'
    label 'cavalier'
    publishDir "output/pull_latest_list", mode: 'copy'
    tag { "$id:v$ver" }

    input:
    tuple val(id), val(ver)

    output:
    tuple val(id), path(output)

    script:
    output = id.replaceAll(':', '_') + '_v' + ver + '.tsv'
    cmd = "cavalier::get_web_list('$id', version='$ver', save='$output', secure=FALSE)"
    """
    R --slave --vanilla -e "$cmd"
    """
}

process update_versions {
    label 'C1M1T1'
    label 'cavalier'
    publishDir "progress/update_list_versions", mode: 'symlink'
    tag { date }

    input:
    tuple path(id_file), val(date)

    output:
    path(output)

    script:
    output = "list_versions_${date}.tsv"
    cmd = "cavalier::get_web_list_version(readr::read_lines('$id_file'), '$output')"
    """
    R --slave --vanilla -e "$cmd"
    """
}