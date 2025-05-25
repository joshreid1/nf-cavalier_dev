
/* ----------- funtions ----------------*/
include { path                } from '../../functions/helpers'
include { ref_fa_channel      } from '../../functions/helpers'
include { get_vcfanno_conf    } from '../../functions/snv_helpers'
include { get_vcfanno_files   } from '../../functions/snv_helpers'
include { get_spliceai_files  } from '../../functions/vep_helpers'
include { get_alphamiss_files } from '../../functions/vep_helpers'
include { get_revel_files     } from '../../functions/vep_helpers'
include { get_utr_ann_files   } from '../../functions/vep_helpers'

// include { get_dbnsfp_files   } from '../../functions/snv_helpers'

/* ----------- workflows ----------------*/
include { CHECK        } from './check'

/* ----------- processes ----------------*/
include { SCATTER      } from '../../modules/local/scatter'
include { CLEAN        } from '../../modules/local/clean'
include { VCFANNO_CONF } from '../../modules/local/vcfanno_conf'
include { VCFANNO      } from '../../modules/local/vcfanno'
include { FILTER       } from '../../modules/local/filter'
include { VEP          } from '../../modules/local/vep'
include { GATHER       } from '../../modules/local/gather'
include { FAM_VARS     } from '../../modules/local/fam_vars'

workflow SNV {

    vcf_channel = Channel.value([path(params.snv_vcf), path(params.snv_vcf + '.tbi')])
    
    if (params.snv_vcfanno) {
        VCFANNO_CONF(
            get_vcfanno_conf()
        )
    }

    CHECK(vcf_channel, 'snv_vcf')

    SCATTER(
        CHECK.out.vcf,
        params.snv_n_shards
    )

    vcf_shards = SCATTER.out.flatMap().map{ [((it.name =~ /(?<=\.shard\.)([0-9]+)/)[0][1]), it] }
    
    CLEAN(
        vcf_shards,
        ref_fa_channel()
    )

    if (params.snv_vcfanno) {

        VCFANNO(
            CLEAN.out,
            VCFANNO_CONF.out,
            get_vcfanno_files()
        )

        if (params.snv_vcfanno_filter) {
            FILTER(
                VCFANNO.out,
                params.snv_vcfanno_filter
            )
            vep_input = FILTER.out
        } else {
            vep_input = VCFANNO.out
        }

    } else {
        vep_input = CLEAN.out.map{ [it[0], it[1]] }
    }

    VEP(
        vep_input,
        ref_fa_channel(),
        path(params.vep_cache),
        get_spliceai_files(),
        get_alphamiss_files(),
        get_revel_files(),
        get_utr_ann_files()
    )

    GATHER(
        VEP.out.toSortedList { a, b -> a.name <=> b.name }
    )

    FAM_VARS (
        GATHER.out,
        CHECK.out.families
    )

    emit:
    FAM_VARS.out.tsv

}

