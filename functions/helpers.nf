import java.text.SimpleDateFormat
import groovy.json.JsonOutput

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
        names = lines[0].split(',') as ArrayList<String>
        assert req_names.all { names.contains(it) }
        // TODO - check opt_names
        lines.each {assert it.split(',').size() == names.size() }
        lines.collect {
            [names, it.split(',')].transpose().collectEntries { k, v -> [(k): v] }
        }
    }
}

String date_ymd() {
    date = new Date()
    sdf = new SimpleDateFormat("yyyy-MM-dd")
    sdf.format(date)
}

void checkMode(mode) {
    if (! mode instanceof String) {
        throw new Exception("ERROR: Mode must be a string")
    }
    valid_modes = ['short', 'sv']
    if (! valid_modes.contains(mode)) {
        throw new Exception("ERROR: Mode must be on of: '${valid_modes.join("', '")}'")
    }
}

def read_ped() {
    ped = read_tsv(path(params.ped), ['fid', 'iid', 'pid', 'mid', 'sex', 'phe'])
    
    ped_families = ped.collect { it.fid }.unique()
    
    fam_w_aff = ped
        .groupBy { it.fid }
        .collect { k, v -> [
            k, v.findAll {it.phe == '2'}.collect {it.iid} ] 
        }
        .findAll { it[1].size() > 0 }
        .collect { it[0] }

    if (fam_w_aff.size() == 0) {
        throw new Exception("No affected individuals in pedigree")
    }
    
    fam_wo_aff = ped_families - fam_w_aff
    
    if (fam_wo_aff.size() > 0) {
        n = fam_wo_aff.size()
        fams =  n > 5 ? fam_wo_aff[0..4] + ['...'] : fam_wo_aff
        println "WARNING: $n famil${n > 1 ? 'ies':'y'} with no affected members will be excluded: ${fams.join(', ')}"
    }

    ped.findAll { fam_w_aff.contains(it.fid) }
}

def read_bams() {
    read_tsv(path(params.bams), ['iid', 'bam'])
}

def list_channels() {

    if (params.lists == null) {
        error("params.lists must not be null")
    }
    lists_list = params.lists.split(',') as ArrayList
    regex = /(HP|PA[A-Z]+|HGNC|G4E):.*/
    web = lists_list.findAll { it ==~ regex }
    local =  lists_list.findAll { !(it ==~ regex) }

    [ 
       (web ? Channel.fromList(web) : null),
       (local ? Channel.fromList(local) : null) 
    ]
}

def get_ref_fa_fai() {
    ref_fa = path(params.ref_fasta)
    ref_fai = path(params.ref_fasta + '.fai')
    Channel.value([ref_fa, ref_fai])
}

def ref_data_channel() {
    ref_fa = path(params.ref_fasta)
    ref_fai = path(params.ref_fasta + '.fai')
    gaps = params.ref_hg38 ?
        path("${workflow.projectDir}/data/hg38.gaps.bed.gz") :
        path("${workflow.projectDir}/data/hg19.gaps.bed.gz")
    vep_cache = path(params.vep_cache)
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

    fam_af_un = read_ped()
        .groupBy { it.fid }
        .collect { k, v -> [
            k,
            v.findAll {it.phe == '2'}.collect {it.iid},
            v.findAll {it.phe == '1'}.collect {it.iid}
        ] }

    Channel.fromList(fam_af_un) // fam, aff, unaff
    
}

def pedigree_channel() {

    Channel.fromList(read_ped()) |
        unique |
        map { it.values() as ArrayList } |
        collectFile(newLine:true) {
            [ "${it[0]}.ped", it.join('\t')]
        } |
        map { [it.name.replaceAll('.ped', ''), it] }
    // fam, ped
}

def bam_channel() {

    Channel.fromList(read_bams()) |
        unique |
        map { [it.iid, path(it.bam), path(it.bam + '.bai')] } |
        combine(read_ped().collect { [it.iid, it.fid] }, by: 0) |
        map { it[[3,0,1,2]] } |
        groupTuple(by: 0)
    // fam, iid, bam, bai
}

def get_options_json(escape='"') {
    options = (params.cavalier_options ?: [:]) +
    [ 
        database_mode: params.database_mode,
        cache_dir : params.cache_dir,
        ref_genome: params.ref_hg38 ? 'hg38' : 'hg19'
    ]
    JsonOutput.toJson(options).replace('"', escape)
}
