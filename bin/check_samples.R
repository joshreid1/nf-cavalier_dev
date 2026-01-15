#!/usr/bin/env Rscript

stopifnot(
  require(tidyverse)
)

# set command line arguments
args      <- commandArgs(trailingOnly = TRUE)
bams      <- args[1]
ped       <- args[2]
short_vcf <- args[3]
struc_vcf <- args[4]
Sys.setenv("VROOM_CONNECTION_SIZE" = 1048576)

bam_ids <- read_tsv(bams, col_names = c("sample_id", "bam_path"), show_col_types = FALSE)$sample_id

if (ped != "UNSET") {
    ped_ids <- read_tsv(ped, col_names = c("fid", "sid", "pid", "mid", "sex", "phe"), show_col_types = FALSE)$sid

    sample_set <- intersect(bam_ids, ped_ids)

    message("\n----- Checking sample set consistency between input files -----\n")

    if (setequal(bam_ids, ped_ids)) {
        message("✔ BAM and PED sample IDs match\n")
    } else if (length(sample_set) > 0) {
        warning(
            "ℹ BAM and PED sample IDs do not match:\n\t",
            "|union| = ", length(union(bam_ids, ped_ids)), "\n\t",
            "|intersection| = ", length(sample_set), "\n",
            "Using union ..."
        )
    } else {
        stop("✖ No matching sample IDs between BAM and PED files\n")
    }
} else {
    sample_set <- bam_ids
    message("\nℹ No PED file provided, using BAM sample IDs only\n")
}

# check short_vcf_ids match sample_set
if (short_vcf != "UNSET") {

    short_vcf_ids <- colnames(read_tsv(short_vcf, comment = "##", n_max = 0, show_col_types = FALSE))[-c(1:9)]
    if (setequal(short_vcf_ids, sample_set)) {
        message("✔ Short VCF and PED/BAM sample IDs match\n")
    } else if (length(intersect(short_vcf_ids, sample_set)) > 0) {
        warning(
            "ℹ Short VCF and PED/BAM sample IDs do not match:\n\t",
            "|union| = ", length(union(short_vcf_ids, sample_set)), "\n\t",
            "|intersection| = ", length(intersect(short_vcf_ids, sample_set)), "\n"
        )
    } else {
        stop("✖ No matching sample IDs between Short VCF and PED/BAM files\n")
    }
} else {
    short_vcf_ids <- NULL
    message("\nℹ No Short VCF file provided, skipping Short VCF sample ID checks\n")
}

# check struc_vcf_ids match sample_set
if (struc_vcf != "UNSET") {
    struc_vcf_ids <- colnames(read_tsv(struc_vcf, comment = "##", n_max = 0, show_col_types = FALSE))[-c(1:9)]

    if (setequal(struc_vcf_ids, sample_set)) {
        message("✔ Struc VCF and PED/BAM sample IDs match\n")
    } else if (length(intersect(struc_vcf_ids, sample_set)) > 0) {
        warning(
            "ℹ Struc VCF and PED/BAM sample IDs do not match:\n\t",
            "|union| = ", length(union(struc_vcf_ids, sample_set)), "\n\t",
            "|intersection| = ", length(intersect(struc_vcf_ids, sample_set)), "\n"
        )
    } else {
        stop("✖ No matching sample IDs between Struc VCF and PED/BAM files\n")
    }

    
} else {
    struc_vcf_ids <- NULL
    message("\nℹ No Struc VCF file provided, skipping Struc VCF sample ID checks\n")
}

if (!is.null(short_vcf_ids) & !is.null(struc_vcf_ids)) {
    # check short_vcf_ids match struc_vcf_ids exactly
    if (setequal(short_vcf_ids, struc_vcf_ids)) {
        message("✔ Short VCF and Struc VCF sample IDs match\n")
    } else {
        stop("✖ Short VCF and Struc VCF sample IDs do not match\n")
    }
}