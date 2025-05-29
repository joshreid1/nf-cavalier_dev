
process CAVALIER_OPTS {
    /*
    create a json file of parameters to pass to cavalier R package
    */
    container null
    label 'C1M1T1'

    output:
    path 'cavalier_opts.json'

    script:
    def cav_opts_txt = groovy.json.JsonOutput.prettyPrint(
          groovy.json.JsonOutput.toJson(params.cavalier_options ?: [:])
    )
"""
cat << 'EOF' > cavalier_opts.json
${cav_opts_txt}
EOF
"""
}
