#!/bin/bash

CHROM=$1
VCF_PREFIX=chroms/${CHROM}/chr${CHROM}.af.cadd.snpEff
OUT_PREFIX=chroms/${CHROM}

zcat ${VCF_PREFIX}.vcf.gz \
    | vawk --header 'I$AF>0' \
    | vawk '{
        split(I$ANN,x,",");
        for(idx in x) { split(x[idx],y,"|");
            print $1,$2,$3,y[4],y[2],I$CADDPHRED }
        }' \
    | zapdups -u \
    | bgzip -c \
    > ${OUT_PREFIX}/impact.txt.gz



