
/* ----------- funtions ----------------*/
include { list_channels    } from '../../functions/helpers'
include { path             } from '../../functions/helpers'
include { date_ymd         } from '../../functions/helpers'
include { cache_dir_channel } from '../../functions/channels.nf'

/* ----------- processes ----------------*/
include { CAVALIER_OPTS } from '../../modules/local/cavalier_opts'
include { INIT_CACHE    } from '../../modules/local/init_cache'
include { MAP_LIST_IDS  } from '../../modules/local/map_list_ids'
include { GET_LIST_VER  } from '../../modules/local/get_list_ver'
include { FETCH_LIST    } from '../../modules/local/fetch_list'


workflow SETUP {
    /*
        - Initialise cavalier
        - Pull down latest versions of gene lists
        - Normalised gene_ids (i.e. convert to ensemble)
    */
    main:
    CAVALIER_OPTS()

    INIT_CACHE(
        date_ymd(),
        CAVALIER_OPTS.out,
        cache_dir_channel()
    )

    cavalier_opts  = INIT_CACHE.out

    list_channel = channel.fromList([])

    (web_lists, local_lists) = list_channels()
    
    if (web_lists != null) {
        update_input = web_lists
            .unique()
            .collectFile(name: 'list_ids.txt', newLine: true)
            .combine([date_ymd()])
        
        GET_LIST_VER(
            update_input,
            cavalier_opts,
            cache_dir_channel()
        )

        FETCH_LIST(
            GET_LIST_VER.out.splitCsv(sep: '\t', skip: 1, strip: true),
            cavalier_opts,
            cache_dir_channel()
        )

        list_channel = FETCH_LIST.out.map { it[1] }
    }

    if (local_lists != null) {

        MAP_LIST_IDS(
            local_lists.map { path(it) },
            cavalier_opts,
            cache_dir_channel()
        )
        list_channel = MAP_LIST_IDS.out.mix(list_channel)
    }

    list_channel = list_channel.toSortedList { a, b -> a.name <=> b.name }.flatMap().collect()

    emit:
    lists  = list_channel
    cavalier_opts = cavalier_opts
}