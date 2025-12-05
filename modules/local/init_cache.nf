
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
    val(cav_opts)
    path(cache_dir)

    output: 
    path('cavalier_options.cache.json')

    script:
"""
cat > cavalier_options.json <<< '${cav_opts}'

cavalier_init_cache.R cavalier_options.json cavalier_options.cache.json
"""
}