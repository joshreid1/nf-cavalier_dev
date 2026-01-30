
process MAKE_SLIDES {
    label 'C2M8T2'
    label 'cavalier'
    publishDir "${params.outdir}/report/$fam", mode: 'copy'
    tag "$fam"
    /*
        - Create slides for candidate variants
        - Loads data directly from FILTER output
    */  

    input:
    tuple val(fam), path(ped), path(short_var), path(igv), path(struc_var), path(svpv), path(samplot)
    path(lists)
    val(slide_info)
    path(cav_opts)
    path(cache_dir)

    output:
    tuple val(fam), path("${fam}.pptx")

    
    script:
"""
cat > slide_info.json <<< '$slide_info'

make_slides.R $ped ${lists.join(',')} slide_options.json \\
    --short-var ${short_var.size() > 0 ? "$short_var"           : 'NONE' } \\
    --struc-var ${struc_var.size() > 0 ? "$struc_var"           : 'NONE' } \\
    --igv       ${igv.size()       > 0 ? "${igv.join(',')}"     : 'NONE' } \\
    --svpv      ${svpv.size()      > 0 ? "${svpv.join(',')}"    : 'NONE' } \\
    --samplot   ${samplot.size()   > 0 ? "${samplot.join(',')}" : 'NONE' } \\
    --slide-info slide_info.json \\
    --cav-opts $cav_opts \\
    --output $fam
"""
}