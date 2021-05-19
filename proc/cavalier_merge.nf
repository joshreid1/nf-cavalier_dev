
params.GTEx_median_rpkm_file = '/stornext/Bioinf/data/lab_bahlo/public_datasets/GTEx/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz'
params.OMIM_genemap2_file = '/wehisan/bioinf/lab_bahlo/ref_db/human/OMIM/OMIM_2019-05-04/genemap2.txt'

process cavalier_merge {
    cpus 1
    memory '4 GB'
    time '4 h'
    container 'jemunro/cavalier:latest'
    publishDir "output/cavalier_merge", mode: 'copy'

    input:
    tuple val(id), file(inputs)

    output:
    tuple val(id), file(out)

    script:
    out = "${id}_qualvars_filtered_exonic.rds"
    """
    cavalier_merge.R $out $inputs
    """
}

