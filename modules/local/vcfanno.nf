
process VCFANNO {
    label 'C2M2T2'
    label 'vcfanno'
    tag "$i"
    /*
        - Run user specified VCFanno annotations
    */ 

    input:
        tuple val(i), path(in_vcf), path(idx)
        path(conf)
        path(ann_files)

    output:
        tuple val(i), path(output)

    script:
    output  = in_vcf.name.replace('.vcf.gz', ".vcfanno.vcf.gz")
    """
    vcfanno -p ${task.cpus} $conf $in_vcf | gzip -c > ${output}
    """
}
