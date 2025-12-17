
process CADDSV {
    label 'C8M8T16'
    label 'caddsv'
    beforeScript 'mkdir -p output && mkdir -p beds'
    tag "$i"
    containerOptions { 
        "--writable-tmpfs " + 
        "--bind $annotations:/CADD-SV/annotations " + 
        "--bind $bed:/CADD-SV/input/$bed " + 
        "--bind beds:/CADD-SV/beds " + 
        "--bind output:/CADD-SV/output" 
    }
    /*
        Not ready for production
    */

    input: 
    tuple val(i), path(bed)
    path(annotations)

    script:
    """
    cd /CADD-SV
    printf "dataset:\\n  - shard_${i}\\n" > input.yaml
    snakemake  --use-conda --configfile input.yaml -j $task.cpus
    """
}