
process REPORT {
    label 'C2M4T2'
    label 'cavalier'
    publishDir "${params.outdir}/cavalier", mode: 'copy'
    tag "$fam"

    input:
    tuple val(fam), path(tsv), path(ped), val(sam), path(bam), path(bai)
    path(config)
    path(cav_opts)
    path(lists)
    path(cache_dir)
    path(func_source)
    
    output:
    tuple val(fam), path("${fam}.snv.pptx"), path("${fam}.snv_candidates.csv"), path("${fam}.snv_filter_stats.csv"),  path("${fam}.snv_filter_reason.csv.gz"), emit: cands
    tuple val(fam), path("${fam}.igv.bed.gz"), emit: igv
    
    script:
    sam_bam = [sam, bam instanceof List ? bam: [bam]]
        .transpose()
        .collect { it.join('=') }
        .join(',')

    """
    report.R $tsv $ped $sam_bam ${lists.join(',')} $config \\
        --out $fam \\
        --cav-opts $cav_opts \\
        --func-source ${func_source.join(',')} \\
        ${params.no_slides ? '--no-slides': ''}
    """
}