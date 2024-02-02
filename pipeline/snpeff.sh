#!/bin/bash

CHROM=$1
VCF_PREFIX=chroms/${CHROM}/chr${CHROM}.af.cadd

java -Xmx16g -jar /net/home/colby/src/snpEff/snpEff.jar \
    -verbose \
    -noStats \
    GRCh37.p13 \
    ${VCF_PREFIX}.vcf.gz \
    | bcftools view -O z \
    > ${VCF_PREFIX}.snpEff.vcf.gz &&

tabix -p vcf -f ${VCF_PREFIX}.snpEff.vcf.gz
