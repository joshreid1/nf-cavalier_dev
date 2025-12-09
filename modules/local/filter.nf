
process FILTER {
    label 'C2M4T4'
    label 'cavalier'
    publishDir "${params.outdir}/report/$fam", mode: 'copy', pattern: "*.csv*"
    tag "$fam"
    /*
        - Read TSV formatted variants for a given family
        - Identify and report candidate variants with cavalier script 
        - Custom filtering and manipulation of variants by modifying functions in ./misc/scripts/filtering_logic.R
        - Currently on SNV/Indel implemented, but will extend to SV
    */  

    input:
    tuple val(fam), path(short_var), path(ped)
    path(lists)
    val(filter_opts)
    path(cav_opts)
    path(cache_dir)

    output:
    tuple val(fam), path("${fam}*.short.filtered_variants.rds") , emit: short_rds
    tuple val(fam), path("${fam}*.short.filtered_variants.csv") , emit: short_csv
    tuple val(fam), path("${fam}*.short.igv.bed.gz")            , emit: short_igv
    tuple val(fam), path("${fam}*.short.count")                 , emit: short_count
    tuple val(fam), path("${fam}*.short.reason_filtered.csv.gz"), emit: short_reason

    tuple val(fam), path("${fam}*.struc.filtered_variants.rds") , emit: struc_rds
    tuple val(fam), path("${fam}*.struc.filtered_variants.csv") , emit: struc_csv
    tuple val(fam), path("${fam}*.struc.count")                 , emit: struc_count
    tuple val(fam), path("${fam}*.struc.reason_filtered.csv.gz"), emit: struc_reason
    
    script:
"""
cat > filter_options.json <<< '${filter_opts}'

filter.R $ped ${lists.join(',')} filter_options.json \\
    --short-var $short_var \\
    --output $fam \\
    --cav-opts $cav_opts
"""
}