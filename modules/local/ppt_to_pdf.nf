
process PPT_TO_PDF {
    container 'linuxserver/libreoffice:7.6.7'
    errorStrategy 'ignore'
    memory '4G'
    tag { "$fam" }
    publishDir "${params.outdir}/cavalier", mode: 'copy', pattern: "*.pdf"
    
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