workflow VcfAnno {
    take: input // tuple val(set), val(i), val(j), path(vcf)

    main:
    
    if (params.vcfanno_config) {
        def configFile = params.vcfanno_config

        processed = input.map { tuple -> tuple + [configFile] } | run_vcfanno
    } else {
        // No extra file, just pass through the input
        processed = input
    }
    
    emit:
    processed
}

process run_vcfanno {
    label 'C4M4T1'
    
    tag { "$set" }

    input:
      tuple val(set), val(i), val(j), path(vcf), val(configFile)

    output:
      tuple val(set), val(i), val(j), path(vcf)

    script:
    vcf_path = vcf
    formatted_vcf = vcf_path.name.replaceAll(/(\.vcf\.gz)|(\.bcf)$/, '.vcf')
    zipped_vcf = vcf_path.name.replaceAll(/(\.vcf\.gz)|(\.bcf)$/, '.vcf.gz')
    annotated_vcf = vcf_path.name.replaceAll(/(\.vcf\.gz)|(\.bcf)$/, '.gnomad.vcf')

    """
    bcftools view ${vcf_path} > ${formatted_vcf}
    bgzip -f ${formatted_vcf} > ${zipped_vcf}
    ${params.vcfanno} ${configFile} ${zipped_vcf} > ${annotated_vcf}
    bcftools view ${annotated_vcf} -o ${vcf_path}
    """
}
