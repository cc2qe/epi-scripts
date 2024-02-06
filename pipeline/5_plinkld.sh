#!/bin/bash

CHROM=$1

VCF_PREFIX=chroms/${CHROM}/chr${CHROM}.af.cadd.snpEff.polarized
PLINK_PREFIX=chroms/${CHROM}/chr${CHROM}.exonic.pol
OUT_PREFIX=chroms/${CHROM}

# A copy of the VCF without the genotypes, so that we can polarize `--a2-allele`
zcat ${VCF_PREFIX}.vcf.gz | cut -f -8 > ${VCF_PREFIX}.noGT.vcf &&

# The VCF file contains multi-allelic sites (on separate lines).
# These should be removed because the multi-allelic LD
# structure that will confound results. I should have addressed this
# earlier in the pipeline but I overlooked it. (Note: plink's `--max-alleles 2`
# does not work because these variants are on multiple lines.
#
# As a workaround, I will follow instructions from
# https://groups.google.com/g/plink2-users/c/O5UdjS2AAW4 to exclude
# the multiallelics based on position by variant id.
cat ${VCF_PREFIX}.noGT.vcf \
    | vawk '{ print $2,$3 }' \
    | zapdups -k1 -v | cut -f 2 \
    > ${VCF_PREFIX}.multiallelic.txt &&

# make impact file, tab delimited file of functional impact by variant
cat ${VCF_PREFIX}.noGT.vcf \
    | vawk --header 'I$AF>0' \
    | vawk '{
        split(I$ANN,x,",");
        for(idx in x) { split(x[idx],y,"|");
            print $1,$2,$3,y[4],y[2],I$CADDPHRED }
        }' \
    | zapdups -u \
    | bgzip -c \
    > ${OUT_PREFIX}/impact.txt.gz &&

# load VCF into plink
plink1.9 \
    --threads 4 \
    --memory 16000 \
    --vcf ${VCF_PREFIX}.vcf.gz \
    --exclude ${VCF_PREFIX}.multiallelic.txt \
    --biallelic-only \
    --a2-allele ${VCF_PREFIX}.noGT.vcf 4 3 '#' \
    --make-bed \
    --out ${PLINK_PREFIX} &&

# compute LD from plink bed file
plink1.9 \
    --threads 4 \
    --memory 16000 \
    --bfile ${PLINK_PREFIX} \
    --keep-allele-order \
    --r gz d 'with-freqs' \
    --ld-window 100000 \
    --ld-window-kb 1000 \
    --ld-window-r2 0 \
    --out ${PLINK_PREFIX} &&

# From the LD file, extract only variant pairs affecting the same gene
zcat ${PLINK_PREFIX}.ld.gz | tr -s ' ' '\t' | sed 's/^\t//g' \
    | zjoin -p "CHR" -a stdin -b <(zcat ${OUT_PREFIX}/impact.txt.gz) -1 3 -2 3 \
    | cut -f -10,15 \
    | zjoin -p "CHR" -a stdin -b <(zcat ${OUT_PREFIX}/impact.txt.gz) -1 7 -2 3 \
    | cut -f -11,15 \
    | awk 'NR==1 { OFS="\t"; print $0,"GENE_A","GENE_B"; next} { if ($(NF-1)==$NF) print }' \
    | zapdups \
    | bgzip -c \
    > ${PLINK_PREFIX}.ld.genematch.gz
