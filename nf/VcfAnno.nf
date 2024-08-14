workflow VcfAnno {
    take: input // tuple val(set), val(i), val(j), path(vcf), path(fasta), path(fai), path(cache)

    main:
    
    if (params.vcfanno_config) {
        def configFile = params.vcfanno_config

        // Process the input with the pass_through_process, including the additional file
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
    
    tag { "vcfanno:$set:$i:$j" }

    input:
      tuple val(set), val(i), val(j), path(vcf), path(fasta), path(fai), path(cache), val(configFile)

    output:
      tuple val(set), val(i), val(j), path(vcf), path(fasta), path(fai), path(cache)

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
