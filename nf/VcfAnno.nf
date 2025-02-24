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
    label 'C1M20T1'
    
    tag { "$set" }

    input:
      tuple val(set), val(i), val(j), path(vcf), val(configFile)

    output:
      tuple val(set), val(i), val(j), path(annotated_vcf_path)

    script:
    vcf_path = vcf
    formatted_vcf = vcf_path.name.replaceAll(/(\.vcf\.gz)|(\.bcf)$/, '.vcf')
    zipped_vcf = vcf_path.name.replaceAll(/(\.vcf\.gz)|(\.bcf)$/, '.vcf.gz')
    annotated_vcf = vcf_path.name.replaceAll(/(\.vcf\.gz)|(\.bcf)$/, '.gnomad.vcf')

    // Define annotated_vcf based on the original file's suffix
    if (vcf_path.name.endsWith('.vcf.gz')) {
        annotated_vcf_path = vcf_path.name.replaceAll(/\.vcf\.gz$/, '.gnomad.vcf.gz')
    } else if (vcf_path.name.endsWith('.bcf')) {
        annotated_vcf_path = vcf_path.name.replaceAll(/\.bcf$/, '.gnomad.bcf')
    }

    """
    bcftools view ${vcf_path} > ${formatted_vcf}
    bgzip -f ${formatted_vcf} > ${zipped_vcf}
    ${params.vcfanno} ${configFile} ${zipped_vcf} > ${annotated_vcf}
    bcftools view ${annotated_vcf} -o ${annotated_vcf_path}
    """
}