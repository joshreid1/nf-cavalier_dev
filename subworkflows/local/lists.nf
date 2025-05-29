
/* ----------- funtions ----------------*/
include { list_channels    } from '../../functions/helpers'
include { path             } from '../../functions/helpers'
include { date_ymd         } from '../../functions/helpers'

/* ----------- processes ----------------*/
include { MAP_IDS    } from '../../modules/local/map_ids'
include { UPDATE_VER } from '../../modules/local/update_ver'
include { FETCH      } from '../../modules/local/fetch'

workflow LISTS {
    /*
        - Pull down latest versions of gene lists
        - Normalised gene_ids (i.e. convert to ensemble)
    */
    take:
    cavalier_opts

    main:
    list_channel = channel.fromList([])

    (web_lists, local_lists) = list_channels()
    
    if (web_lists != null) {
        update_input = web_lists
            .unique()
            .collectFile(name: 'list_ids.txt', newLine: true)
            .combine([date_ymd()])
        
        UPDATE_VER(
            update_input,
            cavalier_opts
        )

        FETCH(
            UPDATE_VER.out.splitCsv(sep: '\t', skip: 1, strip: true),
            cavalier_opts
        )

        list_channel = FETCH.out.map { it[1] }
    }

    if (local_lists != null) {

        MAP_IDS(
            local_lists.map { path(it) },
            cavalier_opts
        )
        list_channel = MAP_IDS.out.mix(list_channel)
    }

    list_channel = list_channel.toSortedList { a, b -> a.name <=> b.name }.flatMap().collect()

    emit:
    list_channel
}