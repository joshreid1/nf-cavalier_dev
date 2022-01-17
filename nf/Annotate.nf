include { path; get_ref_data } from './functions'
include { vep; vep_svo } from './vep'
include { vcf_concat } from './vcf_concat'
include { SplitSVs as SplitPopSVs } from './SplitSVs'

workflow Annotate {
    take:
        vcf_chunks // // set, i, j, vcf, index

    main:
    ref_data = get_ref_data()

    // SNP/Indel VCFs
    if (params.snp_vcf) {
        snp_ann = vcf_chunks |
            filter { it[0] == 'SNP' } |
            map { it[0..3] } |
            combine(ref_data.map { it[[0, 1, 3]] }) |
            vep |
            flatMap {
                [[it[0], 'vep', it[1]],
                 [it[0], 'vep-modifier', it[2]],
                 [it[0], 'unannotated', it[3]]]
            } |
            collectFile(newLine: true, sort: { new File(it).toPath().fileName.toString() }) {
                ["${it[0]}.${it[1]}.files.txt", it[2].toString()]
            } |
            map { (it.name =~ /([^.]+)\.([^.]+)\.files\.txt/)[0][1..2] + [it] } |
            vcf_concat |
            filter { it[1] == 'vep' } |
            map { it[0, 2, 3] }
    } else {
        snp_ann = Channel.fromList([])
    }

//    // SV VCFs
    if (params.sv_vcf) {

        sv_type_match_rev = params.sv_type_match
            .collectMany { k, v -> v.collect { [it, k] }}
            .groupBy { it[0] }
            .collectEntries { k, v -> [(k) : v.collect{ it[1] }.unique()] }

        pop_sv_split = Channel.value(['SV', path(params.pop_sv), path(params.pop_sv + '.tbi')]) |
            SplitPopSVs |
            filter { sv_type_match_rev.keySet().contains(it[0]) } |
            flatMap { sv_type_match_rev[it[0]].collect {type -> [type] + it[1..2] } } |
            view
//            mix(vcf_stub(pop_sv) | map { ['STUB'] + it })

        sv_ann = vcf_chunks |
            filter { it[0] != 'SNP' }
    } else {
        sv_ann = Channel.fromList([])
    }
//    output = mix(snp_ann, sv_ann)
//
//
//    emit: output
}
