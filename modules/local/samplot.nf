
process SAMPLOT {
    label 'C2M4T4'
    label 'samplot'
    tag "$fam"
    publishDir "${params.outdir}/by_family/$fam/samplot", mode: 'copy'

    /*
        - Generate samplot for structural variants
    */

    input:
    tuple val(fam), path(sites), val(ids), path(bams), path(bais)

    output:
    tuple val(fam), path("*.png")

    script:
    """
    export MPLCONFIGDIR=\$PWD/mpl_tmp

    while IFS=\$'\\t' read -r NAME CHROM START END TYPE; do
    samplot plot \\
        -n ${ids.join(' ')} \\
        -b ${bams.join(' ')} \\
        -o samplot_${fam}_\$NAME.png \\
        -c \$CHROM \\
        -s \$START \\
        -e \$END \\
        -t \$TYPE
    done < ${sites}
    """
}