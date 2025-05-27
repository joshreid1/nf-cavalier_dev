
process REPORT_CONF {
    /*
    create a json file of parameters to pass to cavalier R package
    */
    container null

    label 'C1M1T1'

    output:
    path 'report_conf.json'

    script:
    def snv_conf = params
        .keySet()
        .findAll{ it.startsWith('snv_report_') }
        .collectEntries { key -> [(key.replace('snv_report_', '')): params[key]] }
    
    def report_conf = [snv: snv_conf]

    def report_conf_txt = groovy.json.JsonOutput.prettyPrint(
          groovy.json.JsonOutput.toJson(report_conf)
    )
"""
cat << 'EOF' > report_conf.json
${report_conf_txt}
EOF
"""
}
