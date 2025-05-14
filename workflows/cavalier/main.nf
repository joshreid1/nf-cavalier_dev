
// include { vcf_channel; families_channel } from './nf/functions'
// include { GetSamples } from './nf/GetSamples'
// include { CheckInputs } from './nf/CheckInputs'
// include { CleanAndChunk } from './nf/CleanAndChunk'
// include { Annotate } from './nf/Annotate'
// include { Report } from './nf/Report'

include { SNP             } from '../../subworkflows/local/snp'
// include { SV } from '../../subworkflows/local/snp'

workflow CAVALIER {


    if (params.snp_vcf) {
        SNP()
    }
    // if (params.sv_vcf) {
    //     SV()
    // }
    

    /* ---------------OLD CODE ---------------*/
    // vcfs = vcf_channel()

    // vcf_samples = GetSamples(vcfs)

    // ann_vcf = vcf_samples |
    //     CheckInputs |
    //     combine(vcfs, by: 0) |
    //     CleanAndChunk |
    //     Annotate |
    //     combine(families_channel(vcf_samples), by: 0) |
    //     Report
}

