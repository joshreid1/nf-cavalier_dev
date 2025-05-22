
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
