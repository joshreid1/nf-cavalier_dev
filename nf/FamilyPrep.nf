
include { pedigree_channel; bam_channel } from './functions'
include { Lists } from './Lists'


workflow FamilyPrep {
    take:
        input // set, vcf, tbi, fam, aff, unaff

    main:
    output = input |
        map { it[[0,1,3,4,5]] } |
        subset |
        map { it[[1,0,2]] } | //fam, set, vcf
        combine(pedigree_channel(), by:0) |
        combine(bam_channel(), by:0) | //fam, set, vcf, ped, sam, bam, bai
        combine(Lists(), by:0) | //fam, set, vcf, ped, sam, bam, bai, lists
        map { it[[1,0,2,3,7,4,5,6]] } //set, fam, vcf, ped, lists, sam, bam, bai

    emit: output //set, fam, vcf, ped, lists, sam, bam, bai

}

process subset {
    label 'C2M2T2'
    publishDir "output/family_subset", mode: 'copy'
    tag { "$set:$fam" }

    input:
    tuple val(set), path(vcf), val(fam), val(aff), val(unaff)

    output:
    tuple val(set), val(fam), path(out_vcf), path("${out_vcf}.tbi")

    script:
    out_vcf = "${params.id}.${set}.${fam}.subset.vcf.gz"
    """
        printf "${aff.join('\\n')}\\n" > aff
        bcftools view --no-update $vcf -Ou -s ${(aff + unaff).join(',')} |
            bcftools view  -i "GT[@aff]='alt'" -Oz -o $out_vcf
        bcftools index -t $out_vcf
        """
}