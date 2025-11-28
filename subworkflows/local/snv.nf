/* ----------- funtions ----------------*/
include { ref_fasta_channel   } from '../../functions/channels'
include { vcf_channel         } from '../../functions/channels'
include { get_vcfanno_conf    } from '../../functions/vcfanno_helpers'
include { get_vcfanno_files   } from '../../functions/vcfanno_helpers'
include { get_spliceai_files  } from '../../functions/vep_helpers'
include { get_alphamiss_files } from '../../functions/vep_helpers'
include { get_revel_files     } from '../../functions/vep_helpers'
include { get_utr_ann_files   } from '../../functions/vep_helpers'
include { get_vep_fields      } from '../../functions/vep_helpers'
include { get_vep_cache       } from '../../functions/vep_helpers'

/* ----------- processes ----------------*/
include { VCFANNO_CONF } from '../../modules/local/vcfanno_conf'
include { SCATTER      } from '../../modules/local/scatter'
include { CLEAN        } from '../../modules/local/clean'
include { VCFANNO      } from '../../modules/local/vcfanno'
include { VEP          } from '../../modules/local/vep'
include { GATHER       } from '../../modules/local/gather'

workflow SNV {
    take: 
    vcf

    main:
    /*
        - Preprocess and annotate SNV/INDEL variants
    */    
    if (params.snv_vcfanno) {
        VCFANNO_CONF(
            get_vcfanno_conf()
        )
    }

    SCATTER(
        vcf,
        params.snv_n_shards
    )

    vcf_shards = SCATTER.out.flatMap().map{ [((it.name =~ /(?<=\.shard\.)([0-9]+)/)[0][1]), it] }
    
    CLEAN(
        vcf_shards,
        ref_fasta_channel()
    )

    if (params.snv_vcfanno) {

        VCFANNO(
            CLEAN.out,
            VCFANNO_CONF.out,
            get_vcfanno_files(),
            file(params.vcfanno_binary)
        )

        vep_input = VCFANNO.out

    } else {
        vep_input = CLEAN.out.map{ [it[0], it[1]] }
    }

    VEP(
        vep_input,
        ref_fasta_channel(),
        get_vep_cache(),
        get_spliceai_files(),
        get_alphamiss_files(),
        get_revel_files(),
        get_utr_ann_files(),
        get_vep_fields()
    )

    GATHER(
        VEP.out.toSortedList { a, b -> a.name <=> b.name }
    )

    emit:
    GATHER.out

}

