
Path path(filename) {
    file(filename, checkIfExists: true).toAbsolutePath()
}

Path make_path(filename) {
    if (! new File(filename).isDirectory()) {
        new File(filename).mkdir()
    }
    path(filename)
}

ArrayList<Map> read_tsv(Path path, List<String> names ) {
    path.toFile().readLines().with { lines ->
        lines.each { assert it.split('\t').size() == names.size() }
        lines.collect {
            [names, it.split('\t')].transpose().collectEntries { k, v -> [(k): v] }
        }
    }
}

ArrayList<Map> read_csv(Path path, List<String> names ) {
    path.toFile().readLines().with { lines ->
        lines.each {assert it.split(',').size() == names.size() }
        lines.collect {
            [names, it.split(',')].transpose().collectEntries { k, v -> [(k): v] }
        }
    }
}
/* read_csv but with flexible column names and ordering */
ArrayList<Map> read_csv2(Path path, List<String> req_names, List<String> opt_names ) {
    path.toFile().readLines().with { lines ->
        def names = lines[0].split(',') as ArrayList<String>
        assert req_names.all { names.contains(it) }
        // TODO - check opt_names
        lines.each {assert it.split(',').size() == names.size() }
        lines.collect {
            [names, it.split(',')].transpose().collectEntries { k, v -> [(k): v] }
        }
    }
}

String date_ymd() {
    def date = new Date()
    def sdf = new java.text.SimpleDateFormat("yyyy-MM-dd")
    sdf.format(date)
}

def checkMode(mode) {
    if (! mode instanceof String) {
        throw new Exception("ERROR: Mode must be a string")
    }
    def valid_modes = ['short', 'sv']
    if (! valid_modes.contains(mode)) {
        throw new Exception("ERROR: Mode must be on of: '${valid_modes.join("', '")}'")
    }
}

def read_ped() {
    def ped
    if (params.ped) {
        ped = read_tsv(path(params.ped), ['fid', 'iid', 'pid', 'mid', 'sex', 'phe'])
    } else {
        ped = read_bams().collect { it -> [fid: it.iid, iid: it.iid, pid: '0', mid: '0', sex: '0', phe: '2'] }
    }
    
    def ped_families = ped.collect { it.fid }.unique()
    
    def fam_w_aff = ped
        .groupBy { it.fid }
        .collect { k, v -> [
            k, v.findAll {it.phe == '2'}.collect {it.iid} ] 
        }
        .findAll { it[1].size() > 0 }
        .collect { it[0] }

    if (fam_w_aff.size() == 0) {
        throw new Exception("No affected individuals in pedigree")
    }
    
    def fam_wo_aff = ped_families - fam_w_aff
    
    if (fam_wo_aff.size() > 0) {
        def n = fam_wo_aff.size()
        def fams =  n > 5 ? fam_wo_aff[0..4] + ['...'] : fam_wo_aff
        println "WARNING: $n famil${n > 1 ? 'ies':'y'} with no affected members will be excluded: ${fams.join(', ')}"
    }

    ped.findAll { fam_w_aff.contains(it.fid) }
}

def read_bams() {
    read_tsv(path(params.bams), ['iid', 'bam'])
}

def get_external_lists() {
    def lists_list = params.lists.split(',') as ArrayList
    def reg1 = /^(HP|PA[A-Z]+|HGNC|G4E|chr):.*/
    def reg2 = /^[^:]+:[0-9]+-[0-9]+$/
    lists_list.findAll { it ==~ reg1 ||  it ==~ reg2 }
}

def get_local_lists() {
    def lists_list = params.lists.split(',') as ArrayList
    def reg1 = /^(HP|PA[A-Z]+|HGNC|G4E|chr):.*/
    def reg2 = /^[^:]+:[0-9]+-[0-9]+$/
    def loc = lists_list.findAll { !(it ==~ reg1) }.findAll { !(it ==~ reg2) }
    if (loc) {
        Channel.value(loc.collect { path(it) })
    } else {
        Channel.value(path("$projectDir/misc/dummy/local_list"))
    }
}


def list_channels() {

    if (params.lists == null) {
        error("params.lists must not be null")
    }
    def lists_list = params.lists.split(',') as ArrayList
    def regex = /(HP|PA[A-Z]+|HGNC|G4E):.*/
    def web = lists_list.findAll { it ==~ regex }
    def local =  lists_list.findAll { !(it ==~ regex) }

    [ 
       (web ? Channel.fromList(web) : null),
       (local ? Channel.fromList(local) : null) 
    ]
}

def ref_fa_channel() {
    def ref_fa = path(params.ref_fasta)
    def ref_fai = path(params.ref_fasta + '.fai')
    Channel.value([ref_fa, ref_fai])
}

def ref_data_channel() {
    def ref_fa = path(params.ref_fasta)
    def ref_fai = path(params.ref_fasta + '.fai')
    def gaps = params.ref_hg38 ?
        path("${workflow.projectDir}/data/hg38.gaps.bed.gz") :
        path("${workflow.projectDir}/data/hg19.gaps.bed.gz")
    def vep_cache = path(params.vep_cache)
    Channel.value([ref_fa, ref_fai, gaps, vep_cache])
}

def get_variants_override() {
    params.variants_override ? path(params.variants_override) : path("$projectDir/data/dummy/VARIANTS_OVERRIDE")
}

def pop_sv_channel() {
    Channel.value([path(params.pop_sv), path(params.pop_sv + '.tbi')])
}

def ref_gene_channel() {
    Channel.value([path(params.ref_gene)])
}

def families_aff_un() {

    def fam_af_un = read_ped()
        .groupBy { it.fid }
        .collect { k, v -> [
            k,
            v.findAll {it.phe == '2'}.collect {it.iid},
            v.findAll {it.phe == '1'}.collect {it.iid}
        ] }

    Channel.fromList(fam_af_un) // fam, aff, unaff
    
}

def get_cavalier_opts() {
    def cav_opts = (params.cavalier_options ?: [:]) + [cache_dir: params.cavalier_cache_dir]
    groovy.json.JsonOutput.prettyPrint(
          groovy.json.JsonOutput.toJson(cav_opts)
    )
}

def get_filter_opts() {
    def filter_opts = params.findAll { k, v -> k.startsWith('FILTER_') && v != null }
    groovy.json.JsonOutput.prettyPrint(
          groovy.json.JsonOutput.toJson(filter_opts)
    )
}

def get_slide_info() {
    def info = [
        SHORT: params.SLIDE_INFO_SHORT,
        STRUC: params.SLIDE_INFO_STRUC
    ]
    groovy.json.JsonOutput.prettyPrint(
          groovy.json.JsonOutput.toJson(info)
    )
}

def collect_csv(csv_channel, filename) {
    
    def split_csv = csv_channel.splitCsv(header: true)

    split_csv
        .first()
        .map { (it.keySet() as List).join(',') }
        .concat(
            split_csv
                .map { (it.values() as List).join(',') }
                .toSortedList()
                .flatten()
        )
        .collectFile(
            name: filename, 
            storeDir: params.outdir,
            newLine: true,
            sort: false,
            cache: false
        )
}

def short_enabled() {
    params.short_vcf ? true : params.short_vcf_annotated ? true : false
}

def sv_enabled() {
    params.sv_vcf ? true : params.sv_vcf_annotated ? true : false
}



def get_short_fmt() {
    (params.short_format ?: []).toList().unique()
}

def get_short_inf() {
    (
        (params.short_info ?: []) + 
        (params.short_vcfanno ? params.short_vcfanno.collectMany { it.fields.keySet() } : [])
    ).unique()
}

def get_fmt(type) {
    if (type == 'SHORT') {
        return(get_short_fmt())
    }
    if (type == 'STRUC') {
        return(get_struc_fmt())
    }
}

def get_struc_fmt() {
    (params.struc_format ?: []).toList().unique()
}

def get_struc_inf() {
    (
        (params.struc_info ?: []) + 
        (['Max_AF', 'Max_Het', 'Max_HomAlt', 'Max_PopMax_AF', 'ThousG_Count', 'gnomAD_Count', 'CCDG_Count', 'TOPMed_Count']) // from SVAFotate
    ).unique()
}

def get_inf(type) {
    if (type == 'SHORT') {
        return(get_short_inf())
    }
    if (type == 'STRUC') {
        return(get_struc_inf())
    }
}
