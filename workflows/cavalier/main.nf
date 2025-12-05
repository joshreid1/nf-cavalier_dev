
/* ----------- funtions ----------------*/
include { get_filter_opts     } from '../../functions/helpers'
include { collect_csv         } from '../../functions/helpers'
include { cache_dir_channel   } from '../../functions/channels'
include { pedigree_channel    } from '../../functions/channels'
include { bam_channel         } from '../../functions/channels'
include { func_source_channel } from '../../functions/channels'


/* ----------- subworkflows ----------------*/
include { CHECK_VCF } from '../../subworkflows/local/check_vcf'


/* ----------- processes ----------------*/
include { SPLIT_VEP   } from '../../modules/local/split_vep'
include { FILTER      } from '../../modules/local/filter.nf'
include { PPT_TO_PDF  } from '../../modules/local/ppt_to_pdf'
include { IGV_REPORT  } from '../../modules/local/igv_report'
include { IGV_TO_PNG } from '../../modules/local/igv_to_png.nf'

workflow CAVALIER {
    take:
    lists
    cavalier_opts
    short_vcf
    struc_vcf // TODO

    main:
    /*
        - Filter and report variants
    */    
    CHECK_VCF(
        short_vcf,
        'short_vcf'
    )

    SPLIT_VEP(
        short_vcf,
        CHECK_VCF.out.families
    )

    filter_input = SPLIT_VEP.out.tsv // fam, tsv
        .combine(pedigree_channel(), by: 0) //fam, tsv, ped
        // .combine(bam_channel()     , by: 0) // fam, tsv, bed, [sam], [bam], [bai]
    
    FILTER(
        filter_input,
        get_filter_opts(),
        cavalier_opts,
        lists,
        cache_dir_channel()
    )

    short_count = FILTER.out.short_count.map { [it[0], it[1].text as int] }

    IGV_REPORT(
        FILTER.out.short_igv
            .combine(short_count.filter{it[1] > 0}.map{it[0]}, by:0)
            .combine(pedigree_channel(), by:0)
            .combine(SPLIT_VEP.out.vcf, by:0)
            .combine(bam_channel(), by:0)
            
    )
 
    IGV_TO_PNG(
        IGV_REPORT.out.individual
    )

    // PPT_TO_PDF(
    //     REPORT.out.cands.map { [it[0], it[1]] }
    // )

    // collect_csv(
    //     REPORT.out.cands.map { it[2] },
    //     "snp_candidates.csv"
    // )
}

