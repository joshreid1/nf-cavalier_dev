
process INIT_CACHE {
    label 'C1M1T1'
    // label 'cavalier'
    tag "$date_ymd"
    module 'R/flexiblas/4.5.2'
    publishDir "${params.outdir}/gene_lists"
    /*
        - Initialise cavalier cache, checking for latest version of HGNC, HPO, MI_OMIM
    */

    input: 
    val(date_ymd)
    path(cav_opts)
    path(cache_dir)

    output: 
    path(output)

    script:
    output = cav_opts.name.replaceAll('\\.json$', '') + '.cache.json'
    """
    cavalier_init_cache.R $cav_opts $output
    """
}