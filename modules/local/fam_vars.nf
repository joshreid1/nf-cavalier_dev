
process FAM_VARS {
    /*
        Restrict VCF to select samples
        Filter for any variants carried by any affected
        Create TSV output in addition to VCF for simpler parsing in R
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
    inf = (
        (params.snv_info ?: []) + 
        (params.snv_vcfanno ? params.snv_vcfanno.collectMany { it.fields.keySet() } : [])
    ).unique().sort()
    inf_hdr = inf.join('\\t')
    inf_qry = inf.collect {"%$it"}.join('\\t')
    // extract sample FORMAT fields from VCF
    samples = (aff + unaff).unique().sort()
    format = params.snv_format.unique().sort()
    fmt_hdr = format
        .collect { fmt -> samples.collect { sam -> "${fmt}_${sam}" } }
        .flatten()
        .join('\t')
    fmt_query = format.collect{"[%$it\\t]"}.join('')
    """
    printf "${aff.join('\\n')}\\n" > aff
    
    bcftools view --threads ${task.cpus} --no-version --no-update $vcf -Ou -s ${samples.join(',')} |
        bcftools view --no-version --threads ${task.cpus} -i "GT[@aff]='alt'" -Oz -o $out_vcf
    bcftools index -t $out_vcf

    VEP_HDR=\$( bcftools +split-vep $vcf -l | cut -f2- | paste -sd '\\t' - )
    (
        echo -e "CHROM\\tPOS\\tREF\\tALT\\t${inf_hdr}\\t${fmt_hdr}\\t\${VEP_HDR}"
        bcftools +split-vep ${out_vcf} -d -A tab \\
            -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t${inf_qry}\\t${fmt_query}%CSQ\\n' \\
    ) | bgzip > $out_tsv
    """
}