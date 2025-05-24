
include { path } from './helpers'

def spliceai_enabled() {
    (params.vep_spliceai_snv   ? true : false) && 
    (params.vep_spliceai_indel ? true : false)
}

def get_spliceai_files() {
    if (spliceai_enabled()) {
        snv       = path(params.vep_spliceai_snv)
        snv_index = path(params.vep_spliceai_snv + '.tbi')
        ind       = path(params.vep_spliceai_indel)
        ind_index = path(params.vep_spliceai_indel + '.tbi')
    } else {
        snv       = path("$projectDir/misc/dummy/spliceai_snv")
        snv_index = path("$projectDir/misc/dummy/spliceai_snv_index")
        ind       = path("$projectDir/misc/dummy/spliceai_indel")
        ind_index = path("$projectDir/misc/dummy/spliceai_indel_index")
    }
    Channel.value([snv, snv_index, ind, ind_index])
}

def alphamiss_enabled() {
    params.vep_alphamissense ? true : false
}

def get_alphamiss_files() {
    if (alphamiss_enabled()) {
        amiss  = path(params.vep_alphamissense)
        index  = path(params.vep_alphamissense + '.tbi')
    } else {
        amiss  = path("$projectDir/misc/dummy/alphamissense")
        index  = path("$projectDir/misc/dummy/alphamissense_index")
    }
    Channel.value([amiss, index])
}

def revel_enabled() {
    params.vep_revel ? true : false
}

def get_revel_files() {
    if (revel_enabled()) {
        revel  = path(params.vep_revel)
        index  = path(params.vep_revel + '.tbi')
    } else {
        revel  = path("$projectDir/misc/dummy/revel")
        index  = path("$projectDir/misc/dummy/revel_index")
    }
    Channel.value([revel, index])
}

def utr_ann_enabled() {
    params.vep_utr_annotator ? true : false
}

def get_utr_ann_files() {
    if (utr_ann_enabled()) {
        utrann  = path(params.vep_utr_annotator)
    } else {
        utrann  = path("$projectDir/misc/dummy/utr_annotator")
    }
    Channel.value(utrann)
}

def get_vep_fields() {
    def vep_fields = [
          'SYMBOL', 'Gene', 'VARIANT_CLASS', 'Consequence', 'IMPACT', 'Feature_type', 'Feature',
          'BIOTYPE', 'EXON', 'INTRON', 'HGVSc', 'HGVSp', 'HGVSg', 'Amino_acids',
          'HGNC_ID', 'MANE', 'MANE_SELECT', 'MANE_PLUS_CLINICAL', 'CCDS', 'ENSP', 'SIFT', 'PolyPhen', 'CLIN_SIG',
        ] +
        (alphamiss_enabled() ? ['am_class', 'am_pathogenicity'] : [] ) +
        (revel_enabled() ? ['REVEL'] : [] ) +
        (spliceai_enabled() ? [
            'SpliceAI_pred_DP_AG', 'SpliceAI_pred_DP_AL', 'SpliceAI_pred_DP_DG', 'SpliceAI_pred_DP_DL', 
            'SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL', 
            'SpliceAI_pred_SYMBOL'] :
            []
        ) +
        (utr_ann_enabled() ? [
            '5UTR_annotation', '5UTR_consequence', 'Existing_InFrame_oORFs', 'Existing_OutOfFrame_oORFs',
            'Existing_uORFs' ] : 
            []
        )
    vep_fields
}
