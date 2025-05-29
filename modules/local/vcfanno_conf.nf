
process VCFANNO_CONF {
    container null
    label 'C1M1T1'
    /*
        - Write VCFanno configuration
    */  

    input:
        val blocks

    output:
        path 'vcfanno.conf'

    script:
"""
cat << 'VCFANNO_CONF_EOF' > vcfanno.conf
${blocks.join('\n')}
VCFANNO_CONF_EOF
"""
}
