
process JIGV {
    label 'C2M2T2'
    container null
    tag "$fam"
    publishDir "${params.outdir}/cavalier/$fam", mode: 'copy'


    input:
    tuple val(fam), val(proband), path(sites), path(ped), path(bams), path(bais)
    tuple path(ref), path(ref_fai)
    path(jigv_binary)

    output:
    path(output)

    script:
    // println vcfs[0]
    output = "${fam}.variants.jigv.html"
    """
    [ -x $jigv_binary ] || chmod +x $jigv_binary
    
    ./$jigv_binary $bams \\
        --sample $proband \\
        --ped $ped \\
        --sites $sites \\
        --genome-build=hg38 \\
        --fasta $ref \\
        > $output
    """
    // Plu_PK86442.jigv.bed.gz
    

        // --annotation hg38.refGene.bed.gz \ # see: https://github.com/brentp/jigv/wiki/bed12
        // --annotation LCR-hs38.bed.gz \     # specify as many of these as needed.
}