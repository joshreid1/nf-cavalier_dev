
/* ----------- funtions ----------------*/
include { path                } from '../../functions/helpers'
include { get_ref_fa_fai      } from '../../functions/helpers'
include { get_vcfanno_conf    } from '../../functions/snp_helpers'
include { get_vcfanno_files   } from '../../functions/snp_helpers'
include { get_spliceai_files  } from '../../functions/snp_helpers'
include { get_alphamiss_files } from '../../functions/snp_helpers'
include { get_revel_files     } from '../../functions/snp_helpers'
// include { get_dbnsfp_files   } from '../../functions/snp_helpers'


/* ----------- workflows ----------------*/
include { CHECK   } from './check'

/* ----------- processes ----------------*/
include { SCATTER      } from '../../modules/local/scatter'
include { CLEAN        } from '../../modules/local/clean'
include { DROP         } from '../../modules/local/drop'
include { VCFANNO_CONF } from '../../modules/local/vcfanno_conf'
include { VCFANNO      } from '../../modules/local/vcfanno'
include { FILTER       } from '../../modules/local/filter'
include { VEP          } from '../../modules/local/vep'
include { GATHER       } from '../../modules/local/gather'

workflow SNP {

    vcf_channel = Channel.value([path(params.snp_vcf), path(params.snp_vcf + '.tbi')])
    CHECK(vcf_channel, 'snp_vcf')

    SCATTER(
        CHECK.out,
        params.snp_n_shards
    )

    vcf_shards = SCATTER.out.flatMap().map{ [((it.name =~ /(?<=\.shard\.)([0-9]+)/)[0][1]), it] }
    
    CLEAN(
        vcf_shards,
        get_ref_fa_fai()
    )

    if (params.snp_vcfanno) {

        VCFANNO_CONF(
            get_vcfanno_conf()
        )

        VCFANNO(
            CLEAN.out,
            VCFANNO_CONF.out,
            get_vcfanno_files()
        )

        if (params.snp_vcfanno_filter) {
            FILTER(
                VCFANNO.out,
                params.snp_vcfanno_filter
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
        get_ref_fa_fai(),
        path(params.vep_cache),
        get_spliceai_files(),
        get_alphamiss_files(),
        get_revel_files()
    )

    

    GATHER(
        VEP.out.toSortedList { a, b -> a.name <=> b.name }
    )



//     if (params.snp_ann_sources) {
//         DROP(CLEAN.out)
//         ann_sources = Channel.fromList(params.snp_ann_sources).map { [it[0], it[1], path(it[2]),  path(it[2] + it[3])] }
// //   path("${it[2]}.*")]
//         BCFTOOLS_ANN(
//             DROP.out,
//             ann_sources
//         )
//     }

    


}

