#!/usr/bin/env nextflow

include { validate_params } from './functions/validate'

include { SETUP    } from './workflows/setup'
include { ANNOTATE } from './workflows/annotate'
include { CAVALIER } from './workflows/cavalier'

workflow {

    validate_params()

    SETUP()

    ANNOTATE(
        SETUP.out.vcfanno_binary
    )

    if (!params.annotate_only) {
        CAVALIER(
            SETUP.out.gene_set,
            SETUP.out.lists,
            SETUP.out.cavalier_opts,
            ANNOTATE.out.short_vcf,
            ANNOTATE.out.struc_vcf
        )
    }

}
