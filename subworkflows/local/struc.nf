/* ----------- funtions ----------------*/
include { ref_fasta_channel   } from '../../functions/channels'
include { get_vep_cache       } from '../../functions/vep_helpers'
include { get_vep_fields      } from '../../functions/vep_helpers'
include { get_svafdb          } from '../../functions/vep_helpers.nf'


/* ----------- processes ----------------*/
include { SCATTER              } from '../../modules/local/scatter'
include { CLEAN_STRUC as CLEAN } from '../../modules/local/clean'
// include { SV_TO_BED            } from '../../modules/local/sv_to_bed.nf'
// include { CADDSV               } from '../../modules/local/caddsv.nf'
include { SVAFOTATE            } from '../../modules/local/svafotate.nf'
include { VEP_STRUC as VEP     } from '../../modules/local/vep'
include { GATHER               } from '../../modules/local/gather'



workflow STRUC {
    take: 
    vcf

    main:
    /*
        - Preprocess and annotate structural variants variants
    */

    SCATTER(
        vcf,
        params.struc_n_shards
    )

    vcf_shards = SCATTER.out.flatMap().map{ [((it.name =~ /(?<=\.shard\.)([0-9]+)/)[0][1]), it] }
    
    CLEAN(
        vcf_shards,
    )

    SVAFOTATE(
        CLEAN.out,
        get_svafdb()
    )

    

    // SV_TO_BED(
    //     CLEAN.out
    // )

    // CADDSV(
    //     SV_TO_BED.out,
    //     file('/vast/projects/munro_data/cache/cadd-sv/annotations')
    // )

    VEP(
        SVAFOTATE.out,
        ref_fasta_channel(),
        get_vep_cache(),
        get_vep_fields(false)
    )

    GATHER(
        VEP.out.toSortedList { a, b -> a.name <=> b.name }
    )

    emit:
    GATHER.out
}

