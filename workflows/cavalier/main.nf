
/* ----------- funtions ----------------*/
include { path              } from '../../functions/helpers'
include { pedigree_channel  } from '../../functions/helpers'
include { bam_channel       } from '../../functions/helpers'
include { cache_dir_channel } from '../../functions/helpers'
include { ref_fa_channel    } from '../../functions/helpers'
include { get_func_sources  } from '../../functions/helpers'


/* ----------- workflows ----------------*/
include { SNV     } from '../../subworkflows/local/snv'
include { LISTS   } from '../../subworkflows/local/lists'

/* ----------- processes ----------------*/

include { CAVALIER_OPTS } from '../../modules/local/cavalier_opts'
include { REPORT_CONF   } from '../../modules/local/report_conf'
include { REPORT        } from '../../modules/local/report'
include { PPT_TO_PDF    } from '../../modules/local/ppt_to_pdf'
include { IGV_REPORT    } from '../../modules/local/igv_report'

workflow CAVALIER {

    if (params.snv_vcf) {
        SNV()
    }
    // TODO: SVs
    // if (params.sv_vcf) {
    //     SV()
    // }

    CAVALIER_OPTS()

    LISTS(CAVALIER_OPTS.out)

    REPORT_CONF()

    report_input = SNV.out // fam, tsv
        .combine(pedigree_channel(), by:0) //fam, tsv, ped
        .combine(bam_channel(), by:0) // fam, tsv, bed, [sam], [bam], [bai]
    
    REPORT(
        report_input,
        REPORT_CONF.out,
        CAVALIER_OPTS.out,
        LISTS.out,
        cache_dir_channel(),
        get_func_sources()
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

