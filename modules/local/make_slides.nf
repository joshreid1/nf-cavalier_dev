
process MAKE_SLIDES {
    label 'C2M4T4'
    label 'cavalier'
    publishDir "${params.outdir}/report/$fam", mode: 'copy'
    tag "$fam"
    /*
        - Create slides for candidate variants
        - Loads data directly from FILTER output
    */  

    input:
    tuple val(fam), path(short_var), path(ped), path(igv_png)
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
    --short-var $short_var \\
    --igv ${igv_png.join(',')} \\
    --slide-info slide_info.json \\
    --cav-opts $cav_opts \\
    --output $fam
"""
}