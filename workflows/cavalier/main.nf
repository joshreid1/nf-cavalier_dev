
/* ----------- funtions ----------------*/
include { path                } from '../../functions/helpers'
include { get_filter_opts     } from '../../functions/helpers'
include { collect_csv         } from '../../functions/helpers'
include { cache_dir_channel   } from '../../functions/channels'
include { get_slide_info      } from '../../functions/helpers.nf'
include { get_inf             } from '../../functions/helpers.nf'
include { get_fmt             } from '../../functions/helpers.nf'
include { get_fam_aff_un      } from '../../functions/helpers.nf'

/* ----------- processes ----------------*/
include { SPLIT_VEP   } from '../../modules/local/split_vep'
include { FILTER      } from '../../modules/local/filter.nf'
include { IGV_REPORT  } from '../../modules/local/igv_report'
include { IGV_TO_PNG  } from '../../modules/local/igv_to_png.nf'
include { SAMPLOT     } from '../../modules/local/samplot.nf'
include { SVPV        } from '../../modules/local/svpv.nf'
include { SVPV_TO_PNG } from '../../modules/local/svpv_to_png.nf'
include { MAKE_SLIDES } from '../../modules/local/make_slides.nf'
include { PPT_TO_PDF  } from '../../modules/local/ppt_to_pdf'
include { PDF_UNITE   } from '../../modules/local/pdf_unite.nf'
include { PDF_COPY    } from '../../modules/local/pdf_copy.nf'
include { PDF_SPLIT   } from '../../modules/local/pdf_split.nf'

workflow CAVALIER {
    take:
    gene_set
    lists
    cavalier_opts
    short_vcf
    struc_vcf
    pedigree_channel
    bam_channel
    check

    main:
    /*
        - Filter and report variants
    */
    SPLIT_VEP(
        short_vcf.map { ['SHORT'] + it }
            .mix(struc_vcf.map { ['STRUC'] + it })
            .map { it + [get_inf(it[0]), get_fmt(it[0])] }
            .combine(get_fam_aff_un(pedigree_channel, bam_channel)),
        check
    )

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
        .combine(pedigree_channel, by: 0) //fam, short, struc, ped

    FILTER(
        filter_input,
        gene_set,
        get_filter_opts(),
        cavalier_opts,
        cache_dir_channel()
    )

    /* ----- Visualise short variants ----- */
    short_count = FILTER.out.short_count
        .map    { [it[0], it[1].text as int] }
    
    samples_short = short_count
        .filter { it[1] > 0 }
        .map    { it[0]     }
    
    short_count
        .filter { it[1] > params.max_short_per_deck }
        .map    { it[0] }
        .toList()
        .filter { it.size() > 0 }
        .map { 
            log.warn("${it.size()} samples with more than ${params.max_short_per_deck} short variants - slides will be truncated")
        }

    IGV_REPORT(
        FILTER.out.short_igv
            .combine(samples_short, by:0)
            .combine(pedigree_channel, by:0)
            .combine(SPLIT_VEP.out.vcf.filter {it[0] == 'SHORT' }.map { it[[1,2,3]] }, by:0)
            .combine(bam_channel, by:0)
    )
 
    IGV_TO_PNG(
        IGV_REPORT.out.individual
    )

    /* ----- Visualise struct variants ----- */
    struc_count = FILTER.out.struc_count
        .map    { [it[0], it[1].text as int] }

    samples_struc = struc_count
        .filter { it[1] > 0 }
        .map    { it[0]     } 

    struc_count
        .filter { it[1] > params.max_struc_per_deck }
        .map    { it[0] }
        .toList()
        .filter { it.size() > 0 }
        .map { 
            log.warn("${it.size()} samples with more than ${params.max_short_per_deck} struc variants - slides will be truncated")
        }

    SVPV(
        SPLIT_VEP.out.vcf.filter { it[0] == 'STRUC' }.map { it[[1,2]] } // fam, vcf
            .combine(samples_struc, by:0)
            .combine(FILTER.out.struc_csv, by:0) // fam, vcf, csv
            .combine(bam_channel, by:0), // fam, vcf, csv, ids, bams, bais
        path(params.ref_gene)
    )

    SVPV_TO_PNG(
        SVPV.out
    )

    SAMPLOT(
        FILTER.out.struc_samplot
            .combine(samples_struc, by:0)
            .combine(bam_channel, by:0)
    )

    if (params.make_slides) {
         /* ----- Create Slides ----- */
        MAKE_SLIDES(   
            pedigree_channel // fam, ped
                .join(FILTER.out.short_rds.combine(samples_short, by:0), by: 0, remainder: true) // fam, ped, short_rds
                .join(FILTER.out.short_flt_plot, by: 0, remainder: true)
                .join(IGV_TO_PNG.out           , by: 0, remainder: true) // fam, ped, short_rds, short_igv
                .join(FILTER.out.struc_rds.combine(samples_struc, by:0), by: 0, remainder: true) // fam, ped, short_rds, short_igv, struc_rds
                .join(FILTER.out.struc_flt_plot, by: 0, remainder: true)
                .join(SVPV_TO_PNG.out          , by: 0, remainder: true) // fam, ped, short_rds, short_igv, struc_rds, svpv
                .join(SAMPLOT.out              , by: 0, remainder: true) // fam, ped, short_rds, short_igv, struc_rds, svpv, samplot
                .map    { row -> row.collect { x -> x ?: [] }   },  // replace missing/null with [] to avoid errors
            lists,
            get_slide_info(),
            cavalier_opts,
            cache_dir_channel()
        )

        PPT_TO_PDF(
            MAKE_SLIDES.out
        )

        /* ----- Create PDFs by gene ----- */
        PDF_SPLIT(
            PPT_TO_PDF.out
        )

        by_gene = 
            PDF_SPLIT.out.flatten()
            .map  { [(it.name =~ /([^.]+)\.pdf/)[0][1], it] }
            .groupTuple()

        PDF_UNITE(
            by_gene.filter { it[1].size() >  1 }
        )

        PDF_COPY(
            by_gene.filter { it[1].size() == 1 }
        )

        by_gene
            .map { [it[0], it[1].size()] }
            .toSortedList { it[0] }
            .map { "SYMBOL,n_samples\n" + it.collect { it.join(',') }.join('\n') }
            .collectFile(name: 'by_gene_counts.csv',storeDir: "${params.outdir}")
    }

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
