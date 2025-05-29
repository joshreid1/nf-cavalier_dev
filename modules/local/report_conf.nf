
process REPORT_CONF {
    container null
    label 'C1M1T1'
    /*
    create a JSON file of parameters to pass to cavalier R package
    */
    
    input:
    val(report_conf)

    output:
    path 'report_conf.json'

    script:
    def report_conf_txt = groovy.json.JsonOutput.prettyPrint(
          groovy.json.JsonOutput.toJson(report_conf)
    )
"""
cat << 'EOF' > report_conf.json
${report_conf_txt}
EOF
"""
}
