include { vep; vep_svo } from './vep'
include { vcf_concat } from './vcf_concat'

workflow Annotate {
    take:
        vcf_chunks // // set, i, j, vcf, index
        ref_data // ref_fa, ref_fai, gaps, vep_cache
        sv_ref_data // pop_sv, pop_sv_tbi, ref_gene

    main:
    // SNP/Indel VCFs
    if (params.snp_vcf) {
        snp_ann = vcf_chunks |
            filter { it[0] == 'SNP' } |
            map { it[0..3] } |
            combine(ref_data.map { it[[0, 1, 3]] }) |
            vep |
            flatMap {
                [
                    [it[0], 'vep', it[1]],
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
//    if (params.sv_vcf) {
//        sv_ann = vcf_chunks |
//            filter { it[0] != 'SNP' }
//    } else {
//        sv_ann = Channel.fromList([])
//    }
//    output = mix(snp_ann, sv_ann)
//
//
//    emit: output
}
