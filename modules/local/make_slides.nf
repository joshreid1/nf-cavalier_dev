
process MAKE_SLIDES {
    label 'C2M16T4'
    label 'cavalier'
    publishDir "${params.outdir}/by_family/$fam", mode: 'copy', pattern: "${fam}.pptx"
    tag "$fam"
    /*
        - Create slides for candidate variants
        - Loads data directly from FILTER output
    */

    input:
    tuple val(fam), path(ped),
      path(short_var), path(short_flt_plot), path(igv),
      path(struc_var), path(struc_flt_plot), path(svpv), path(samplot)
    path(lists)
    val(slide_info)
    path(cav_opts)
    path(cache_dir)

    output:
    tuple val(fam), path("*.pptx")

    
    script:
"""
cat > slide_info.json <<< '$slide_info'

make_slides.R $ped ${lists.join(',')} slide_options.json \\
    --short-var ${short_var.size() > 0 ? "$short_var"           : 'NONE' } \\
    --short-flt-plot $short_flt_plot \\
    --struc-var ${struc_var.size() > 0 ? "$struc_var"           : 'NONE' } \\
    --struc-flt-plot $struc_flt_plot \\
    --igv       ${igv.size()       > 0 ? "${igv.join(',')}"     : 'NONE' } \\
    --svpv      ${svpv.size()      > 0 ? "${svpv.join(',')}"    : 'NONE' } \\
    --samplot   ${samplot.size()   > 0 ? "${samplot.join(',')}" : 'NONE' } \\
    --slide-info slide_info.json \\
    --cav-opts $cav_opts \\
    --output $fam \\
    --max-short-per-deck $params.max_short_per_deck \\
    --max-struc-per-deck $params.max_struc_per_deck
"""
}