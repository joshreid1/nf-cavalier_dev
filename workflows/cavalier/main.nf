
/* ----------- funtions ----------------*/
include { get_report_conf     } from '../../functions/helpers'
include { cache_dir_channel   } from '../../functions/channels'
include { vcf_channel         } from '../../functions/channels'
include { pedigree_channel    } from '../../functions/channels'
include { bam_channel         } from '../../functions/channels'
include { func_source_channel } from '../../functions/channels'

/* ----------- subworkflows ----------------*/
include { SNV           } from '../../subworkflows/local/snv'
include { LISTS         } from '../../subworkflows/local/lists'
include { CHECK         } from '../../subworkflows/local/check'

/* ----------- processes ----------------*/
include { FAM_VARS      } from '../../modules/local/fam_vars'
include { CAVALIER_OPTS } from '../../modules/local/cavalier_opts'
include { REPORT_CONF   } from '../../modules/local/report_conf'
include { REPORT        } from '../../modules/local/report'
include { PPT_TO_PDF    } from '../../modules/local/ppt_to_pdf'
include { IGV_REPORT    } from '../../modules/local/igv_report'

workflow CAVALIER {
    /*
        - Run nf-cavalier main workflow
    */
    if (params.snv_vcf_annotated) {
        println "INFO: Skipping SNV annotation, using annotated VCF - $params.snv_vcf_annotated"
        
        CHECK(
            vcf_channel(params.snv_vcf_annotated),
            'snv_vcf_annotated'
        )
        snv = CHECK.out.vcf

    } else if (params.snv_vcf) {
        
        CHECK(
            vcf_channel(params.snv_vcf),
            'snv_vcf'
        )
        SNV(
            CHECK.out.vcf
        )
        snv = SNV.out
    }
    
    FAM_VARS(
        snv,
        CHECK.out.families
    )

    CAVALIER_OPTS()

    LISTS(CAVALIER_OPTS.out)

    REPORT_CONF(
        get_report_conf()
    )

    report_input = FAM_VARS.out.tsv // fam, tsv
        .combine(pedigree_channel(), by: 0) //fam, tsv, ped
        .combine(bam_channel()     , by: 0) // fam, tsv, bed, [sam], [bam], [bai]
    
    REPORT(
        report_input,
        REPORT_CONF.out,
        CAVALIER_OPTS.out,
        LISTS.out,
        cache_dir_channel(),
        func_source_channel()
    )

    PPT_TO_PDF(
        REPORT.out.cands.map { [it[0], it[1]] }
    )

    IGV_REPORT(
        REPORT.out.igv
            .combine(pedigree_channel(), by:0)
            .combine(bam_channel().map { [it[0], it[2], it[3]]}, by:0)
    )
}

