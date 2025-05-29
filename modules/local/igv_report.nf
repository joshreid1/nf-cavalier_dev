
process IGV_REPORT {
    label 'C2M2T2'
    label 'igvreports'
    tag "$fam"
    publishDir "${params.outdir}/report/$fam", mode: 'copy'
    /*
        - Generate html igv-reports for candidate variants identified by cavalier
    */

    input:
    tuple val(fam), path(sites), path(ped), path(bams), path(bais)

    output:
    path(output)

    script:
    output = "${fam}.igv_report.html"
    """
    create_report $sites \\
        --genome hg38 \\
        --flanking 100 \\
        --tracks $bams \\
        --output $output
    """
}