
def validate_params() {
    /*
        TODO: more validation of params
    */
    if (!params.snv_vcf & !params.sv_vcf){
        throw new Exception("ERROR: Must specify at least one of 'params.snv_vcf' or 'params.sv_vcf'")
    }
}