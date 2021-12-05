
process svpv {
    label 'C2M4T2'
    container null
    conda '/stornext/Home/data/allstaff/m/munro.j/miniconda3/envs/numpy2'
    publishDir "output/svpv", mode: 'copy'
    tag { fam }

    input:
        tuple val(fam), path(vcf), val(sam), path(bam), path(bai), path(pop_sv), path(pop_sv_indx), path(ref_gene)

    output:
        tuple val(fam), path(fam)

    script:
        sam_bam = [sam, bam instanceof List ? bam: [bam]]
            .transpose().collect {it.join('=') }.join(' ')
        """
        SVPV \\
            -o $fam \\
            -samples ${sam.join(',')} \\
            -aln ${bam.join(',')} \\
            -vcf ${params.caller_id}:$vcf \\
            -ref_vcf gnomAD:$pop_sv \\
            -ref_gene $ref_gene
        """
}
