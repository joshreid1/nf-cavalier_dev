
process PPT_TO_PDF {
    container 'linuxserver/libreoffice:7.6.7'
    label 'C2M4T4'
    maxRetries 3
    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }    
    tag { "$fam" }
    publishDir "${params.outdir}/report/$fam", mode: 'copy', pattern: "*.pdf"
    /*
        - Convert powerpoints to PDFs
        - Error set to ignore because sometimes randomly fails
    */

    input:
    tuple  val(fam), path(pptx)

    output:
    tuple val(fam), path(pdf)

    script:
    pdf = pptx.name.replaceAll('.pptx', '.pdf')
    """
    HOME=\$PWD soffice --headless --convert-to pdf $pptx || echo 'done'
    """
}