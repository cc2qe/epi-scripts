#!/bin/bash

CHROM=$1

mkdir -p chroms/${CHROM} &&

# load bgen file
plink2 \
    --threads 1 \
    --memory 19500 \
    --bgen imp_genome_files/ukb22828_c${CHROM}_b0_v3.bgen ref-first \
    --sample imp_genome_files/ukb22828_c${CHROM}_b0_v3_s487202.sample \
    --keep samples_unrel_caucasian.txt \
    --snps-only \
    --mac 1 \
    --extract range exons_protein_coding.GRCh37.bed.gz \
    --chr ${CHROM} \
    --make-pgen \
    --out chroms/${CHROM}/chr${CHROM} &&

# get sample counts
plink2 \
    --threads 1 \
    --memory 10000 \
    --pfile chroms/${CHROM}/chr${CHROM} \
    --sample-counts \
    --out chroms/${CHROM}/chr${CHROM} &&

# recode variant names (chrom_pos_ref_alt rather than rsid)
PREFIX=chroms/${CHROM}/chr${CHROM} &&
cp ${PREFIX}.pvar ${PREFIX}.pvar.bak &&
cat ${PREFIX}.pvar \
    | awk 'NR==1 { print; next } {print $1,$2,$1"_"$2"_"$4"_"$5,$4,$5}' OFS="\t" \
    > ${PREFIX}.pvar.recode &&
mv ${PREFIX}.pvar.recode ${PREFIX}.pvar


