
process IGV_TO_PNG {
    label 'C2M2T2'
    label 'puppeteer'
    maxForks 100
    tag "$fam"
    publishDir "${params.outdir}/report/$fam/snapshot", mode: 'copy'
    beforeScript "mkdir -p home"
    containerOptions "--overlay \"\$PWD/home\":/home/pptruser"
    /*
        - Generate PNGs for igv-report
    */

      input:
    tuple val(fam), path(reports)

    output:
    tuple val(fam), path("*.png")

    script:
    def width = 600
    def height = 800
    def scale = 3

    cmds = reports.collect{
      def id = (it.name =~ /igv_report\.(.*?)\.html/)[0][1]
      "export_igv_png.js $it ${id} $width $height $scale"
    }.join('\n')

"""
export NODE_PATH=/home/pptruser/node_modules
export HOME=/home/pptruser/

$cmds
"""
}