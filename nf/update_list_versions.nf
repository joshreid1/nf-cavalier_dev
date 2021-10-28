
process update_list_versions {
    label 'C1M1T1'
    label 'cavalier'
//    container = null
//    module = 'R/3.6.1'
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

