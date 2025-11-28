

include { read_ped        } from '../../functions/helpers'
include { read_bams       } from '../../functions/helpers'
include { families_aff_un } from '../../functions/helpers'


include { SAMPLES } from '../../modules/local/samples'


workflow CHECK {
    /*
        - CHECK that VCF samples line up with pedigree and bam input samples
        - Error if there is no work to be done
    */
    take:
    vcf_channel
    name
    
    main:
    def ped = read_ped()
    def bams = read_bams()

    // check samples
    def ped_samples = ped.collect { it.iid }.unique()
    def bam_samples = bams.collect { it.iid }.unique()

    vcf_samples = SAMPLES(vcf_channel).splitText().map{ it.trim() }.toSortedList()

    vcf_samples \
        | flatMap { samples -> 
                [
                    ["in \"${name}\" and \"bams\" but not in \"ped\"", samples.intersect(bam_samples) - ped_samples],
                    ["in \"${name}\" and \"ped\" but not in \"bams\"", samples.intersect(ped_samples) - bam_samples],
                    ["in \"${name}\" but not in \"ped\" or \"bams\"",  samples - ped_samples.intersect(bam_samples)],
                    ["in \"bams\" but not in \"${name}\"",  bam_samples - samples],
                ]
          } \
        | filter { it[1].size() > 0 } \
        | map { warn, sm ->
            def n = sm.size()
            sm = n > 5 ? sm[0..4] + ['...'] : sm
            println "WARNING: $n sample${n > 1 ? 's':''} $warn: ${sm.join(', ')}"
        }

    // //     check there is any work to be done
    vcf_output = vcf_channel
        .combine(vcf_samples.map {[it]})
        .map { vcf, tbi, samples ->
            def complete = samples.intersect(bam_samples).intersect(ped_samples)
            if (complete.size() == 0) {
                throw new Exception("ERROR: No samples to process in $set VCF")
            }
            return [vcf, tbi]
        }
        .first()

    fam_output = vcf_samples
        .map { [it] }
        .combine(families_aff_un())
        .map { sam, fam, af, un -> [fam, af.intersect(sam), un.intersect(sam)] }
        .filter { it[1].size() > 0 }

    emit: 
    vcf      = vcf_output
    families = fam_output
}