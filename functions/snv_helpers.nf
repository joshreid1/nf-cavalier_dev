
include { path } from './helpers'

def get_vcfanno_conf() {
    Channel
        .fromList(params.snv_vcfanno)
        .map { ann ->
            // choose file key and toml keyword
            def file   = file(ann.vcf ?: ann.tsv).name
            def is_vcf = ann.containsKey('vcf')
            def key    = is_vcf ? 'fields' : 'columns'
            def op     = ann.op ?: 'self'
            // build the three arrays
            def values = ann.fields.values().collect { it instanceof String ? "\"${it}\"" : "${it}" }.join(', ')
            def names  = ann.fields.keySet().collect { "\"$it\"" }.join(', ')
            def ops    = ann.fields.collect { "\"$op\"" }.join(', ')
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
        .collect()
}

def get_vcfanno_files() {
    Channel
        .fromList(params.snv_vcfanno)
        .flatMap { [path(it.vcf ?: it.tsv), path("${it.vcf ?: it.tsv}.${it.csi ? 'csi' : 'tbi'}")] }
        .toSortedList()
}

def get_spliceai_files() {
    if (params.vep_spliceai_snv && params.vep_spliceai_indel) {
        snv       = path(params.vep_spliceai_snv)
        snv_index = path(params.vep_spliceai_snv + '.tbi')
        ind       = path(params.vep_spliceai_indel)
        ind_index = path(params.vep_spliceai_indel + '.tbi')
    } else {
        snv       = path("$projectDir/misc/dummy/spliceai_snv")
        snv_index = path("$projectDir/misc/dummy/spliceai_snv_index")
        ind       = path("$projectDir/misc/dummy/spliceai_indel")
        ind_index = path("$projectDir/misc/dummy/spliceai_index")
    }
    Channel.value([snv, snv_index, ind, ind_index])
}

def get_alphamiss_files() {
    if (params.vep_alphamissense) {
        amiss  = path(params.vep_alphamissense)
        index  = path(params.vep_alphamissense + '.tbi')
    } else {
        amiss  = path("$projectDir/misc/dummy/alphamissense")
        index  = path("$projectDir/misc/dummy/alphamissense_index")
    }
    Channel.value([amiss, index])
}

def get_revel_files() {
    if (params.vep_revel) {
        revel  = path(params.vep_revel)
        index  = path(params.vep_revel + '.tbi')
    } else {
        revel  = path("$projectDir/misc/dummy/revel")
        index  = path("$projectDir/misc/dummy/revel_index")
    }
    Channel.value([revel, index])
}

// def get_dbnsfp_files() {
//     if (params.vep_dbnsfp) {
//         dbnsfp       = path(params.vep_dbnsfp)
//         dbnsfp_index = path(params.vep_dbnsfp + '.tbi')
//     } else {
//         dbnsfp       = path("$projectDir/misc/dummy/dbnfsp")
//         dbnsfp_index = path("$projectDir/misc/dummy/dbnfsp_index")
//     }
//     Channel.value([dbnsfp, dbnsfp_index])
// }