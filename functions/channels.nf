include { read_ped  } from './helpers'
include { read_bams } from './helpers'
include { make_path } from './helpers'

// helper functions to create input channels for processes

def file_channel(filename) {
    Channel.value(file(filename, checkIfExists: true).toAbsolutePath())
}
// return various HTS files with extensions
def index_channel(filename, extension) {
    Channel.value(
        [
            file(filename, checkIfExists: true).toAbsolutePath(), 
            extension instanceof List ? 
                (extension.collect { ext -> 
                    file(filename + ".$ext", checkIfExists: true).toAbsolutePath()
                }) :
                file(filename + ".$extension", checkIfExists: true).toAbsolutePath()
        ]
    )
}

def tabix_channel(filename) {
    index_channel(filename, 'tbi')
}

def vcf_channel(filename, extension = 'tbi') {
   index_channel(filename, extension)
}

def fasta_channel(filename) {
   index_channel(filename, 'fai')
}

def ref_fasta_channel() {
    fasta_channel(params.ref_fasta)
}

def func_source_channel() {
    def src_list = ["$projectDir/bin/snv_functions.R"] +
        (
            params.report_func_source ? 
            (params.report_func_source.split(',')) :
            []
        )

    Channel.value(
        src_list.collect { file(it, checkIfExists: true).toAbsolutePath() }
    )
}

def pedigree_channel() {

    Channel.fromList(read_ped())
        .unique()
        .map { it.values() as ArrayList }
        .collectFile(newLine: true) {
            [ "${it[0]}.ped", it.join('\t')]
        }
        .map { [it.name.replaceAll('.ped', ''), it] }
    // fam, ped
}

def bam_channel() {

    Channel.fromList(read_bams())
        .unique()
        .map { 
            [
                it.iid,
                file(it.bam, checkIfExists: true).toAbsolutePath(),
                file(it.bam + '.bai', checkIfExists: true).toAbsolutePath()
            ] 
        }
        .combine(read_ped().collect { [it.iid, it.fid] }, by: 0)
        .map { it[[3,0,1,2]] }
        .groupTuple(by: 0)  
    // fam, iid, bam, bai
}

def cache_dir_channel() {
    if (params.cavalier_options.cache_dir) {
        Channel.value(make_path(params.cavalier_options.cache_dir))
    } else {
        throw new Exception("Please define params.cavalier_options.cache_dir")
    }
}





