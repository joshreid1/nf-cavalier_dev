
process REPORT_CONF {
    /*
    create a json file of parameters to pass to cavalier R package
    */
    container null

    label 'C1M1T1'

    output:
    path 'report_conf.json'

    script:
    def report_conf = [
        snv_freq_filters: params.snv_freq_filters,
        snv_mutate      : params.snv_mutate,
        snv_subsets     : params.snv_subsets,
        snv_report      : params.snv_report,
        snv_min_depth   : params.snv_min_depth
    ]
    def report_conf_txt = groovy.json.JsonOutput.prettyPrint(
          groovy.json.JsonOutput.toJson(report_conf)
    )
"""
cat << 'EOF' > report_conf.json
${report_conf_txt}
EOF
"""
}
