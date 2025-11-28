
process VCFANNO {
    label 'C2M2T2'
    label 'bcftools'
    tag "$i"
    /*
        - Run specified VCFanno annotations
        - optionally filter output
    */ 

    input:
        tuple val(i), path(in_vcf), path(idx)
        path(conf)
        path(ann_files)
        path(vcfanno)

    output:
        tuple val(i), path(output)

    script:
    output  = in_vcf.name.replace('.vcf.gz', ".vcfanno.vcf.gz")
    filter = params.snv_vcfanno_filter ? "-i '$params.snv_vcfanno_filter'" : ''

    """
    [ -x "$vcfanno" ] || chmod +x "$vcfanno"
    
    ./$vcfanno -p ${task.cpus} $conf $in_vcf \\
        | bcftools view --no-version --threads ${task.cpus} $filter -Oz -o $output
    """
}
