
process IGV_REPORT {
    label 'C2M2T2'
    label 'igvreports'
    tag "$fam"
    publishDir "${params.outdir}/report/$fam", mode: 'copy', pattern: '*.igv_report.html'
    /*
        - Generate html igv-reports for candidate variants identified by cavalier
        - Firstly for all samples combined with VCF 
        - Secondly individually to extract PNGs to put into slides
    */

    input:
    tuple val(fam), path(sites), path(ped), path(vcf), path(tbi), val(ids), path(bams), path(bais)

    output:
    tuple val(fam), path(output)                    , emit: combined
    tuple val(fam), path("${fam}.igv_report.*.html"), emit: individual

    script:
    output = "${fam}.igv_report.html"
    """
    create_report $sites \\
        --genome hg38 \\
        --flanking 250 \\
        --tracks $vcf $bams \\
        --output $output

    IDS=(${ids.join(' ')})
    BAM=($bams)
    for i in "\${!IDS[@]}"; do
    create_report $sites \\
        --standalone \\
        --genome hg38 \\
        --flanking 100 \\
        --tracks \${BAM[\$i]} \\
        --output "${fam}.igv_report.\${IDS[\$i]}.html"
    done
    """
}