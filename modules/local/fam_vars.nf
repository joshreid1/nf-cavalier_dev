
process FAM_VARS {
    /*
        - Restrict VCF to select samples
        - Filter for only variants carried by an affected individual
        - Create TSV output in addition to VCF for simpler parsing in R using bcftools +split-vep
    */
    label 'C2M2T2'
    label 'bcftools'
    tag "$fam"

    input:
    tuple path(vcf), path(index)
    tuple val(fam), val(aff), val(unaff)

    output:
    tuple val(fam), path(out_vcf), path("${out_vcf}.tbi"), emit: vcf
    tuple val(fam), path(out_tsv)                        , emit: tsv

    script:
    out_vcf = vcf.name.replace('.vcf.gz', ".family_${fam}.vcf.gz")
    out_tsv = vcf.name.replace('.vcf.gz', ".family_${fam}.tsv.gz")
    // extract INFO fields from VCF
    def inf = (
        (params.snv_info ?: []) + 
        (params.snv_vcfanno ? params.snv_vcfanno.collectMany { it.fields.keySet() } : [])
    ).unique().sort()
    def inf_hdr = inf.join('\\t')
    def inf_qry = inf.collect {"%$it"}.join('\\t')
    // extract sample FORMAT fields from VCF
    def samples = (aff + unaff).unique().sort()
    def format = params.snv_format.toList().unique().sort()
    def fmt_hdr = format
        .collect { fmt -> samples.collect { sam -> "FMT_${fmt}_${sam}" } }
        .flatten()
        .join('\t')
    def fmt_query = format.collect{"[%$it\\t]"}.join('')
    """
    printf "${aff.join('\\n')}\\n" > aff
    
    bcftools view --threads ${task.cpus} --no-version --no-update $vcf -Ou -s ${samples.join(',')} |
        bcftools view --no-version --threads ${task.cpus} -i "GT[@aff]='alt'" \\
            --write-index=tbi -Oz -o $out_vcf

    VEP_HDR=\$( bcftools +split-vep $vcf -l | cut -f2- | paste -sd '\\t' - )
    (
        echo -e "CHROM\\tPOS\\tREF\\tALT\\t${inf_hdr}\\t${fmt_hdr}\\t\${VEP_HDR}"
        bcftools +split-vep ${out_vcf} -d -A tab \\
            -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t${inf_qry}\\t${fmt_query}%CSQ\\n' \\
    ) | bgzip > $out_tsv
    """
}