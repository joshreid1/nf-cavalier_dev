
/* ----------- funtions ----------------*/
include { path                } from '../../functions/helpers'
include { get_filter_opts     } from '../../functions/helpers'
include { collect_csv         } from '../../functions/helpers'
include { cache_dir_channel   } from '../../functions/channels'
include { pedigree_channel    } from '../../functions/channels'
include { bam_channel         } from '../../functions/channels'
include { func_source_channel } from '../../functions/channels'
include { get_slide_info      } from '../../functions/helpers.nf'
include { get_inf             } from '../../functions/helpers.nf'
include { get_fmt             } from '../../functions/helpers.nf'
include { get_fam_aff_un      } from '../../functions/helpers.nf'
include { get_struc_inf       } from '../../functions/helpers.nf'
include { get_struc_fmt       } from '../../functions/helpers.nf'

/* ----------- processes ----------------*/
include { SPLIT_VEP   } from '../../modules/local/split_vep'
include { FILTER      } from '../../modules/local/filter.nf'
include { PPT_TO_PDF  } from '../../modules/local/ppt_to_pdf'
include { IGV_REPORT  } from '../../modules/local/igv_report'
include { IGV_TO_PNG  } from '../../modules/local/igv_to_png.nf'
include { MAKE_SLIDES } from '../../modules/local/make_slides.nf'
include { SAMPLOT     } from '../../modules/local/samplot.nf'

workflow CAVALIER {
    take:
    gene_set
    lists
    cavalier_opts
    short_vcf
    struc_vcf
    check

    main:
    /*
        - Filter and report variants
    */
    SPLIT_VEP(
        short_vcf.map { ['SHORT'] + it }
            .mix(struc_vcf.map { ['STRUC'] + it })
            .map { it + [get_inf(it[0]), get_fmt(it[0])] }
            .combine(get_fam_aff_un()),
        check
    )

    pedigree = pedigree_channel()

    filter_input = SPLIT_VEP.out.tsv
        .filter {it[0] == 'SHORT' }
        .map { it[[1,2]] }
        .join(
            SPLIT_VEP.out.tsv
                .filter {it[0] == 'STRUC' }
                .map { it[[1,2]] },
            remainder: true
        )
        .map { [it[0], it[1] ?: [], it[2] ?: [] ] }
        .combine(pedigree, by: 0) //fam, short, struc, ped

    FILTER(
        filter_input,
        gene_set,
        get_filter_opts(),
        cavalier_opts,
        cache_dir_channel()
    )

    /* ----- Visualise short variants ----- */

    short_count = FILTER.out.short_count.map { [it[0], it[1].text as int] }

    IGV_REPORT(
        FILTER.out.short_igv
            .combine(short_count.filter{it[1] > 0}.map{it[0]}, by:0)
            .combine(pedigree, by:0)
            .combine(SPLIT_VEP.out.vcf.filter {it[0] == 'SHORT' }.map { it[[1,2,3]] }, by:0)
            .combine(bam_channel(), by:0)
    )
 
    IGV_TO_PNG(
        IGV_REPORT.out.individual
    )

    /* ----- Visualise struct variants ----- */

    struc_count = FILTER.out.struc_count.map { [it[0], it[1].text as int] }

    SAMPLOT(
        FILTER.out.struc_bamplot
            .combine(struc_count.filter{it[1] > 0}.map{it[0]}, by:0)
            .combine(bam_channel(), by:0)
    )

    /* ----- Create Slides ----- */

    MAKE_SLIDES(    
        FILTER.out.short_rds
            .combine(pedigree, by:0)
            .combine(IGV_TO_PNG.out, by:0), // fam, vars, ped, igv_imgs
            lists,
            get_slide_info(),
            cavalier_opts,
            cache_dir_channel()
    )

    PPT_TO_PDF(
        MAKE_SLIDES.out
    )

    /* ----- Save CSV results ----- */

    collect_csv(
        FILTER.out.short_csv.map { it[1] },
        'short_candidates.csv'
    )

    collect_csv(
        FILTER.out.struc_csv.map { it[1] },
        'struc_candidates.csv'
    )
}

