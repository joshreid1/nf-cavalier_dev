
process SVPV {
    label 'C2M4T4'
    label 'svpv'
    tag "$fam"
    publishDir "${params.outdir}/by_family/$fam/svpv", mode: 'copy'

    /*
        - Generate SVPV for structural variants
    */

    input:
    tuple val(fam), path(vcf), val(lines), val(ids), path(bams), path(bais)
    path(ref_gene)

    output:
    tuple val(fam), path("$fam/**.pdf")

    script:
    """
    LINES="${lines.split('\n').join(' ')}"
    
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
        -o $fam \\
        -samples ${ids.join(',')} \\
        -aln ${bams.join(',')} \\
        -vcf filtered.vcf.gz \\
        -ref_gene $ref_gene
    """
}