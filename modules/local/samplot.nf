
process SAMPLOT {
    label 'C2M4T4'
    label 'samplot'
    tag "$fam"
    publishDir "${params.outdir}/by_family/$fam/samplot", mode: 'copy'

    /*
        - Generate samplot for structural variants
    */

    input:
    tuple val(fam), val(sites), val(ids), path(bams), path(bais)

    output:
    tuple val(fam), path("*.png")

    script:
"""
export MPLCONFIGDIR=\$PWD/mpl_tmp
cat > sites <<< '${sites}'
(
    while IFS=\$'\\t' read -r NAME CHROM START END TYPE; do

    WINDOW_FLAG=""
    if [[ "\$TYPE" == "INS" && "\$START" == "\$END" ]]; then
        WINDOW_FLAG="-w 1000"
    fi

    echo "samplot plot \\
        -n ${ids.join(' ')} \\
        -b ${bams.join(' ')} \\
        -o samplot_${fam}_\$NAME.png \\
        -c \$CHROM \\
        -s \$START \\
        -e \$END \\
        -t \$TYPE \\
        \$WINDOW_FLAG"
    done < sites
) | xargs -n 1 -P ${task.cpus} -I {} sh -c '{}'
"""
}