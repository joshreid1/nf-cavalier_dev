
/* ----------- funtions ----------------*/
include { path             } from '../../functions/helpers'
include { get_ref_fa_fai   } from '../../functions/helpers'

/* ----------- workflows ----------------*/
include { CHECK   } from './check'

/* ----------- processes ----------------*/
include { SCATTER      } from '../../modules/local/scatter'
include { CLEAN        } from '../../modules/local/clean'
include { DROP         } from '../../modules/local/drop'
include { BCFTOOLS_ANN } from '../../modules/local/bcftools_ann'
include { VCFANNO      } from '../../modules/local/vcfanno'

def get_vcfanno_conf() {
    Channel
        .fromList(params.snp_vcfanno)
        .map { ann ->
            // choose file key and toml keyword
            def file   = file(ann.vcf ?: ann.tsv).name
            def is_vcf = ann.containsKey('vcf')
            def key    = is_vcf ? 'fields' : 'columns'
            // build the three arrays
            def values = ann.fields.values().collect { it instanceof String ? "\"${it}\"" : "${it}" }.join(', ')
            def names  = ann.fields.keySet().collect { "\"$it\"" }.join(', ')
            def ops    = ann.fields.collect { '"self"' }.join(', ')
            // one TOML block per annotation
            def block =  
                "[[annotation]]\n" +
                "file    = \"${file}\"\n" +
                "${key}  = [${values}]\n" +
                "names   = [${names}]\n" +
                "ops     = [${ops}]"
            // println block
            return block
        }
        .collectFile(name: 'vcfanno.conf', sort: false, newLine: true)
        .first()
}

def get_vcfanno_files() {
    Channel
        .fromList(params.snp_vcfanno)
        .flatMap { [path(it.vcf ?: it.tsv), path("${it.vcf ?: it.tsv}.${it.index}")] }
        .collect()
}



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
        VCFANNO(
            CLEAN.out,
            get_vcfanno_conf(),
            get_vcfanno_files()
        )
    }

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

