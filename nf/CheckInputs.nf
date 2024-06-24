include { read_ped; read_bams; list_channels } from './functions'

workflow CheckInputs {
    take: vcf_samples

    main:
    ped = read_ped()
    bams = read_bams()

    // check families
    ped_families = ped.collect { it.fid }.unique()
    ped_w_aff =
        ped.groupBy { it.fid }
        .collect { k, v -> [
            k,  v.findAll {it.phe == '2'}.collect {it.iid} ] }
        .findAll { it[1].size() > 0 }
        .collect { it[0] }

    [['with no affected members', ped_families - ped_w_aff]]
        .findAll { it[1].size() > 0}
        .forEach { warn, fams ->
            n = fams.size()
            fams = n > 5 ? fams[0..4] + ['...'] : fams
            println "WARNING: $n famil${n > 1 ? 'ies':'y'} $warn: ${fams.join(', ')}"
        }

    // check samples
    ped_samples = ped.collect { it.iid }.unique()
    bam_samples = bams.collect { it.iid }.unique()

    vcf_samples |
        map { set, samples -> ["in \"bams\" but not in \"$set\"", bam_samples - samples] } |
        mix(vcf_samples
            .flatMap { set, samples ->
                ne = samples - (bam_samples + ped_samples)
                [["in $set VCF and \"bams\" but not in \"ped\"", samples - ne - ped_samples],
                 ["in $set VCF and \"ped\" but not in \"bams\"", samples - ne - bam_samples],
                 ["in $set VCF but not in \"ped\" or \"bams\"", ne],]}) |
        filter { it[1].size() > 0 } |
        map { warn, sm ->
            n = sm.size()
            sm = n > 5 ? sm[0..4] + ['...'] : sm
            println "WARNING: $n sample${n > 1 ? 's':''} $warn: ${sm.join(', ')}"
        }

    //     check there is any work to be done
    sets =
        vcf_samples |
        map { set, samples ->
            complete = samples.intersect(bam_samples)
            if (complete.size() == 0) {
                throw new Exception("ERROR: No samples to process in $set VCF")
            }
            return set
        }

    emit: sets
}