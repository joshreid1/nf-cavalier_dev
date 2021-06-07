
params.gtex_rpkm = '/stornext/Bioinf/data/lab_bahlo/public_datasets/GTEx/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz'
params.omim_genemap2 = '/stornext/Bioinf/data/lab_bahlo/ref_db/human/OMIM/OMIM_2020-04-29/genemap2.txt'
params.maf_dom = 0.0001
params.maf_rec = 0.01
params.maf_comp_het = 0.01

process cavalier_singleton {
    cpus 1
    memory '4 GB'
    time '4 h'
    container 'jemunro/cavalier:dev'
    publishDir "output/cavalier_singleton", mode: 'copy'
    tag { sample }

    input:
    tuple val(sample), file(vcf), file(bam), file(bai), val(lists)

    output:
    tuple val(sample), file(out)

    script:
    out = "$sample"
    """
    cavalier_singleton_panelapp.R $vcf $out $sample=$bam \\
        --gene-lists $lists \\
        --maf-dom $params.maf_dom \\
        --maf-rec $params.maf_rec \\
        --maf-comp-het $params.maf_comp_het \\
        --gtex-rpkm $params.gtex_rpkm \\
        --omim-genemap2 $params.omim_genemap2
    """
}
//process cavalier_singleton {
//    cpus 1
//    memory '4 GB'
//    time '4 h'
//    container 'jemunro/cavalier:dev'
//    publishDir "output/cavalier_singleton", mode: 'copy'
//
//    input:
//    tuple val(sample), file(vcf), file(bam), file(bai), file(lists)
//
//    output:
//    tuple val(sample), file(out)
//
//    script:
//    out = "$sample"
//    """
//    cavalier_singleton.R $vcf $out $sample=$bam \\
//        --gene-lists ${lists.join(',')} \\
//        --maf-dom $params.maf_dom \\
//        --maf-rec $params.maf_rec \\
//        --maf-comp-het $params.maf_comp_het \\
//        --gtex-rpkm $params.gtex_rpkm \\
//        --omim-genemap2 $params.omim_genemap2
//    """
//}

