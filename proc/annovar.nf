
params.ANNOVAR_DIR = '/stornext/Bioinf/data/lab_bahlo/software/apps/annovar/annovar_2018-04-16'
params.ANNOVAR_SCRIPT = "$params.ANNOVAR_DIR/table_annovar.pl"
params.ANNOVAR_HUMANDB_DIR ="$params.ANNOVAR_DIR/humandb/"

process annovar {
    cpus 1
    memory '2 GB'
    time '1 h'
    module 'perl/5.22.2:htslib/1.12'
    publishDir "output/annovar", mode: 'copy'

    input:
    tuple val(id), file(vcf), file(tbi)

    output:
    tuple val(id), file(out_vcf), file("${out_vcf}.tbi")

    script:
    out_base = "${id}.annovar"
    out_vcf = "${out_base}.hg19_multianno.vcf.gz"
    """
    perl $params.ANNOVAR_SCRIPT $vcf \\
        $params.ANNOVAR_HUMANDB_DIR -buildver hg19 \\
        -vcfinput -out $out_base -remove \\
        -protocol refGene,genomicSuperDups,exac03,gnomad_exome,gnomad_genome,avsnp150,dbnsfp35a,clinvar_20180603 -operation g,r,f,f,f,f,f,f -nastring . \\
        -arg '-exonicsplicing -splicing 8',,,,,,,
    sed 's:avsnp150:avsnp147:g' ${out_base}.hg19_multianno.vcf | \\
        bgzip > $out_vcf
    tabix -p vcf $out_vcf
    gzip ${out_base}.avinput
    gzip ${out_base}.hg19_multianno.txt
    """
}

