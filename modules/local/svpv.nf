
process SVPV {
    label 'C2M4T4'
    label 'svpv'
    tag "$fam"
    /*
        - Generate samplot for structural variants
    */

    input:
    tuple val(fam), path(vcf), path(flt_svs), val(ids), path(bams), path(bais)
    path(ref_gene)
    path(pop_sv)

    output:
    tuple val(fam), path(output)

    script:
    output = "$fam"
    """
    LINES=\$(awk -F',' 'NR==1{for(i=1;i<=NF;i++) if(\$i=="LINE_ID") c=i; if(!c) exit 1} NR>1{print \$c}' $flt_svs)
    
    zcat $vcf \\
        | awk -v lines="\$LINES" '
            BEGIN { split(lines, lst, " "); for(i in lst) keep[lst[i]] = 1 }
            {
                if(/^##/ || /^#CHROM/) { print; next }   # print VCF header
                dataLine++                               # count only data lines
                if(dataLine in keep) print
            }' \\
        | gzip -c > filtered.vcf.gz

    SVPV \\
        -o $output \\
        -samples ${ids.join(',')} \\
        -aln ${bams.join(',')} \\
        -vcf filtered.vcf.gz \\
        -ref_vcf gnomAD_v4:$pop_sv \\
        -ref_gene $ref_gene
    """
}