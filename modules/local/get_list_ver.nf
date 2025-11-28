
process GET_LIST_VER {
    label 'C1M1T1'
    label 'cavalier'
    tag { date }
    /*
        - Get latest version of gene lists specified
        - runs at most once per day (cached based on date)
    */  

    input:
    tuple path(id_file), val(date)
    path(cav_opts)
    path(cache_dir)

    output:
    path(output)

    script:
    output = "list_versions_${date}.tsv"
    cmd = "cavalier::get_gene_list_versions(readr::read_lines('$id_file'), '$output')"
    """
    R --slave --vanilla -e "\\
        cavalier::set_options_from_json('$cav_opts');\\
        cavalier::get_gene_list_versions(readr::read_lines('$id_file'), '$output')\\
    "
    """
}