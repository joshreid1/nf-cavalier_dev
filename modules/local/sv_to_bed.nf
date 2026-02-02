
process SV_TO_BED {
    label 'C2M2T2'
    label 'bcftools'
    tag { i }
    /*
        Output DELS/DUPs as BED file
    */

    input:
    tuple val(i), path(vcf), path(tbi)

    output:
    tuple val(i), path(output)

    script:
    output = "id_shard_${i}.bed"
    """
    bcftools view $vcf --no-version -Ou -i "SVTYPE='DEL' || SVTYPE='DUP'" \\
        | bcftools query -f "%CHROM\\t%POS\\t%INFO/END\\t%INFO/SVTYPE\\n" \\
        | sort -k1,1 -k2,2n -k3,3n -k4,4 \\
        > tmp.bed

    # CADDSV failing on CNVnator bounds extending past chrom length
    grep -w "\$(head -n 1 tmp.bed | cut -f 1)" tmp.bed > $output
    """
}
