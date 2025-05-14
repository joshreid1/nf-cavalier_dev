
void validate_params() {
    if (!params.snp_vcf & !params.sv_vcf){
        throw new Exception("ERROR: Must specify at least one of 'params.snp_vcf' or 'params.sv_vcf'")
    }
}