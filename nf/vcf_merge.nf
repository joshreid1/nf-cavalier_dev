
process vcf_merge {
    label 'C2M2T2'
    publishDir "output/vcf_merge", mode: 'copy'

    input:
    path(file_list)

    output:
    tuple val(set), path(out_vcf), path("${out_vcf}.csi")

    script:
    set = file_list.name.replaceAll('.files.txt', '')
    out_vcf = "${params.id}.${set}.bcf"
    """
    bcftools concat --naive-force --file-list $file_list -Ob -o $out_vcf
    bcftools index --threads 2 $out_vcf
    """
}
