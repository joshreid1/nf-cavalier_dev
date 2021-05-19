
params.GTEx_median_rpkm_file = '/stornext/Bioinf/data/lab_bahlo/public_datasets/GTEx/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz'
params.OMIM_genemap2_file = '/wehisan/bioinf/lab_bahlo/ref_db/human/OMIM/OMIM_2019-05-04/genemap2.txt'

process cavalier_prep {
    cpus 1
    memory '4 GB'
    time '4 h'
    container 'jemunro/cavalier:latest'
    publishDir "output/cavalier_prep", mode: 'copy'

    input:
    tuple val(id), file(vcf), file(tbi)

    output:
    file(out)

    script:
    out = "${id}_qualvars_filtered_exonic.rds"
    """
    cavalier_prep.R $vcf ${id}
    """
}

