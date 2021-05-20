
params.VCFANNO_EXEC='/stornext/Bioinf/data/lab_bahlo/software/apps/vcfanno/vcfanno_v0.3.1/vcfanno'
params.VCFANNO_TOML='/stornext/Bioinf/data/lab_bahlo/software/apps/vcfanno/config/conf_gnomADv2.1.toml'

process vcfanno {
    cpus 8
    memory '4 GB'
    time '6 h'
    module 'htslib/1.12'
    publishDir "output/vcfanno", mode: 'copy'
    maxForks 10

    input:
    tuple val(id), file(vcf)

    output:
    tuple val(id), file(out_vcf), file("${out_vcf}.tbi")

    script:
    out_vcf = "${id}.vcfanno.vcf.gz"
    """
    $params.VCFANNO_EXEC -p ${task.cpus} $params.VCFANNO_TOML $vcf | 
        bgzip > $out_vcf
    tabix -p vcf $out_vcf
    """
}
