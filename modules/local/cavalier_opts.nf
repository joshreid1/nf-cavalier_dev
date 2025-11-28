
process CAVALIER_OPTS {
    /*
    create a json file of parameters to pass to cavalier R package
    */
    container null
    label 'C1M1T1'

    output:
    path 'cavalier_opts.json'

    script:
    def opts = (params.cavalier_options ?: [:]) + [cache_dir: params.cavalier_cache_dir]
    def cav_opts_txt = groovy.json.JsonOutput.prettyPrint(
          groovy.json.JsonOutput.toJson(opts)
    )
"""
cat << 'EOF' > cavalier_opts.json
${cav_opts_txt}
EOF
"""
}
