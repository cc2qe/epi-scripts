#!/bin/bash

CHROM=$1

VCF_PREFIX=chroms/${CHROM}/chr${CHROM}.af.cadd.snpEff.polarized
PLINK_PREFIX=chroms/${CHROM}/chr${CHROM}.exonic.pol

# A copy of the VCF without the genotypes, so that we can polarize `--a2-allele`
zcat ${VCF_PREFIX}.vcf.gz | cut -f -8 > ${VCF_PREFIX}.noGT.vcf &&

# load VCF into plink
plink1.9 \
    --threads 1 \
    --memory 8000 \
    --vcf ${VCF_PREFIX}.vcf.gz \
    --biallelic-only \
    --a2-allele ${VCF_PREFIX}.noGT.vcf 4 3 '#' \
    --make-bed \
    --out ${PLINK_PREFIX} &&

plink1.9 \
    --threads 4 \
    --memory 16000 \
    --bfile ${PLINK_PREFIX} \
    --keep-allele-order \
    --r gz d 'with-freqs' \
    --ld-window 100000 \
    --ld-window-kb 1000 \
    --ld-window-r2 0 \
    --out ${PLINK_PREFIX}
