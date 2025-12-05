
def validate_params() {
    /*
        TODO: more validation of params
    */
    if (!params.short_vcf & !params.sv_vcf){
        throw new Exception("ERROR: Must specify at least one of 'params.short_vcf' or 'params.sv_vcf'")
    }
}